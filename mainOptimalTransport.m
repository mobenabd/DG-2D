clear; clc; close all;
addpath './src_DG/' './optimal_transport/' './obstacle_problem'
%%%%%%%%%%%%% Solve L2MK optimal transport problem %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%profile on
global k x0 xN hxy
n     =  8;          %% grid size in 1D direction
x0 = 0; xN = 1;      %% domain [x0,xN]^2
hxy = (xN-x0)/n;     %% grid size step
sigma =  10;         %% DG penalization parameter
eps   =  1;          %% DG symmetrization parameter
k    =  4;           %% polynomial degree 
gdim = n^2*(k+1)^2;  %% global dimention (degrees of freedom vector length)
c    = 1e4;          %% c parameter of the fixed point algorithm
alpha = 0.5;%1e-2;   %% gradient step size.
itmax = 200;         %% maxi number of gradient iterations
b1 = 0.3; b2 = 0.75; %% Armijo-Goldstein rule constants
a1 = 0.8; a2 = 1.25;



beta     = @(x,y) 0.5*((x-0.5).^2 + (y-0.5).^2);  %% -> 0.5 |xi|^2
beta_xi = computeDirectProjection(n,k,beta);      %% L2 projection in the discrete space

%%%%%%% densities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TC1 %%%%%%%%%%%%
rho_0 = @(x,y) 1 + 0.8*sin(2*pi*2*x)*sin(2*pi*2*y);
rho_1 = @(x,y) 1 + 0.4*cos(2*pi*4*x)*cos(2*pi*4*y);
%%% TC0 %%%%%%%%%%%%
% nu0 = 0.4; nu1 = 0.6;
% ss = 1/100;
% rho_0 = @(x,y)  1/(2*pi*ss) * exp(-0.5*(ss^-1)*((x-nu0)^2 + (y-nu0)^2));
% rho_1 = @(x,y)  1/(2*pi*ss) * exp(-0.5*(ss^-1)*((x-nu1)^2 + (y-nu1)^2)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


rho0_proj = computeDirectProjection(n,k,rho_0);  %% L2 projection in the discrete space
rho1_proj = computeDirectProjection(n,k,rho_1);
%OT_plot_sol(n,k, rho0_proj);


%%%%%% initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Phi_xi_0 = zeros(gdim,1);

Phi_xi = Phi_xi_0;
Alpha_xi = beta_xi + Phi_xi;   %-> the convex function ψ = 0.5 |xi|^2 + Phi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

res=1;      %% L2 relative residual of the displacement ∇(Phi_xi)
tol=1e-4;   %% stoping criterion
it=1;       %% iteration increment
E=0; dE=0;  %% H1-norm (Energy functional)
GRAD=zeros(gdim,2);
while (res > tol)
    %%%%%%%  compute Tilde rho0_xi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %compute det(I + \nabla^2_xi \Phi(\xi))
    det_xi = computeDet_ddl(n,k,Alpha_xi);    
    
    %plot_sol2(n,k,det_xi)
    fprintf('min(det)=%5e\n', min(det_xi));
    
    if(~isempty(find(det_xi<0, 1)))
        fprintf('Warning: det<0\n'); %break;
        det_xi(det_xi<0) = 0;
    end
    
    rho0_tilde_xi = OT2_compute_rho0Tilde_xi(n,k, rho_1, Phi_xi).*det_xi;
    
    % normalizing density
    vol =  integrate(rho0_tilde_xi,n);
    rho0_tilde_xi = 1/vol * rho0_tilde_xi ;
    OT_plot_sol(n,k, rho0_tilde_xi);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%% solve poisson equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % source term
    f_src = @(xi1,xi2) -compute_sol(xi1, xi2, rho0_tilde_xi, n, k) + rho_0(xi1, xi2);
    % Neumann b.c 
    g_bc = @(xi1, xi2) 0;
    if (it == 1)
        [phi_xi, Ldecomp, Udecomp, Mass] = solve_poisson(f_src, g_bc, n, sigma, eps);
    else
        phi_xi = solve_poisson2(f_src, g_bc, n, sigma, eps, Ldecomp, Udecomp, Mass);
    end   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%% update Phi_xi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dPhi_xi = -alpha*phi_xi;
    Phi_xi = Phi_xi + dPhi_xi;
    Alpha_xi = beta_xi + Phi_xi;
    %%% convexify ψ (-concavization)
    % Dirichlet b.c 
    g_bc = @(xi1,xi2) - compute_sol(xi1,xi2, Alpha_xi, n, k) ;
    acc = solve_concav(g_bc, -Alpha_xi, n, sigma, eps, c); 
    
    Alpha_xi = -acc; %%convexification = -acc
    Phi_xi = Alpha_xi - beta_xi;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%% check energy minimization and update gradient step %%%%%%%%%%%%
    Eold = E;
    dEold = dE;
    [E, dE] = OT2_compute_energy(phi_xi, rho_1, Phi_xi, n);
    if (it>1)
        fprintf('E=%5e | dE=%5e | alpha=%5e| (Eold-E)/(alpha*dE)=%f\n', E, dE, alpha, (Eold-E)/(alpha*dEold));
        %%% update step size
        if ( E-Eold > -b1*alpha*dEold )
            alpha = a1*alpha;
        elseif ( E-Eold < -b2*alpha*dEold )
            alpha = a2*alpha;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%% compute cost (L2 wasserstein distance) and residual %%%%%%%%%%%
    cost =  OT_compute_cost(rho_0, Phi_xi, n);
    GRADold = GRAD;
    % compute the displacement ∇(Phi_xi)
    GRAD    = computeGrad_ddl(n,k,Phi_xi);
    res     = sqrt(sum((GRAD(:,1)-GRADold(:,1)).^2+(GRAD(:,2)-GRADold(:,2)).^2))/sqrt(sum(GRAD(:,1).^2+GRAD(:,2).^2));
    % L2 error between rho0 and Tile rho0
    l2error = integrate((rho0_proj-rho0_tilde_xi).^2,n)^0.5;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('volume=%5e, cost =%5e, res_grad=%5e, l2-error=%5e, counter=%i\n', vol, cost, res, l2error, it); 
    it = it +1;
    
    
end

OT_plot_displ(n,k,rho0_proj,Phi_xi)


