clear; clc; close all;
%%%%%%%%%%%%% Solve optimal transport problem %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% test with initializing \delta\Phi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%profile on
global k
n     =  8;
sigma =  10; eps   =  1;
k    =  4; gdim = n^2*(k+1)^2;
global x0 xN hxy
x0 = -2; xN = 2; hxy = (xN-x0)/n;
h=hxy;


c    = 1e4;
alpha = 0.5;%1e-2;
itmax = 200;
b1 = 0.3; b2 = 0.75;
a1 = 0.8; a2 = 1.25;


beta     = @(x,y) 0.5*((x-0.5).^2 + (y-0.5).^2);
beta_old = @(x,y) 0.5*((x-0.5).^2 + (y-0.5).^2);
beta_xi = computeDirectProjection(n,k,beta_old);

%%%%%%% densities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TC1 %%%%%%%%%%%%
% rho_0 = @(x,y) 1 + 0.8*sin(2*pi*2*x)*sin(2*pi*2*y);
% rho_1 = @(x,y) 1 + 0.4*cos(2*pi*4*x)*cos(2*pi*4*y);
%%% TC0 %%%%%%%%%%%%
nu0 = -0.4; nu1 = 0.4;
ss = 1/8;
rho0_0 = @(x,y)  1/(2*pi*ss) * exp(-0.5*(ss^-1)*((x-nu0)^2 + (y-nu0)^2)) ;
rho0_1 = @(x,y)  1/(2*pi*ss) * exp(-0.5*(ss^-1)*((x-nu1)^2 + (y-nu1)^2)) ;

%rho0_0 = @(x,y) 3.97887*exp(-12.5*(x^2 +(y +0.1)^2 ));
%rho0_1 = @(x,y) 4.59441*exp(-16.66666*(x^2 +x*(y-0.1)+(y-0.1)^2));


rho0_proj = computeDirectProjection(n,k,rho0_0);
gamma = 0.1;
eps_k = gamma*max(rho0_proj)/(1-gamma);

rho_0 = @(x,y) rho0_0(x,y) + eps_k;
rho_1 = @(x,y) rho0_1(x,y) + eps_k;


vol0 = integrate2(rho_1,n);
rho_0 = @(x,y)  rho_0(x,y);
rho_1 = @(x,y)  rho_1(x,y);

rho0_proj = computeDirectProjection(n,k,rho_0);
rho1_proj = computeDirectProjection(n,k,rho_1);



% OT_plot_sol(n,k, rho0_proj./rho1_proj);
%  return
% OT_plot_sol(n,k, rho0_proj);
% hold on
% OT_plot_sol(n,k, rho1_proj);
% return
% %%%%%%%%%%%%%%%%%%%



%%%%%% initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Phi_xi_0 = zeros(gdim,1);

%%%%%%%%%%%%%%%
% alpha_xy = @(x,y) 0.5*(rho_0(x,y) + rho_1(x,y));
% f_source = @(x,y) rho_1(x,y) - rho_0(x,y);
% g_bc = @(x, y) 0;
%
% Phi_xi_0 = solve_poisson(f_source, g_bc, n, sigma, eps, alpha_xy);
%
%
% Phi_xi = Phi_xi_0;
%
% Alpha_xi = beta_xi + Phi_xi;
%
%
% %%%%%% convexification
% g_bc = @(xi1,xi2) -( beta(xi1,xi2) );  %with b.c homogeneous Direchlet
% acc = solve_concav(g_bc, -Alpha_xi, n, sigma, eps, c); %%convexification = -ac
%
% Alpha_xi = -acc;
% Phi_xi = Alpha_xi - beta_xi;

%%%%%%%%%%%%%%%

Phi_xi = Phi_xi_0;
Alpha_xi = beta_xi + Phi_xi;


res=1; tol=1e-4;
it=1; E=0; dE=0;
GRAD=zeros(gdim,2);
l2error = 1;
while (res > tol)
    det_xi = computeDet_ddl(n,k,Alpha_xi);     %%det(I + \nabla^2_xi \Phi(\xi))
    
    %plot_sol2(n,k,det_xi)
    fprintf('min_det=%5e\n', min(det_xi));
    
    if(~isempty(find(det_xi<0, 1)))
        fprintf('Warning: det<0\n'); %break;
        det_xi(det_xi<0) = 0;
    end
    
    
    rho0_tilde_xi = OT2_compute_rho0Tilde_xi(n,k, rho_1, Phi_xi).*det_xi;
    
    
    
    vol =  integrate(rho0_tilde_xi,n);
    rho0_tilde_xi = 1/vol * rho0_tilde_xi * vol0;
    %OT_plot_sol(n,k, rho0_tilde_xi);
    
    
    %%% solve poisson equation
    f_src = @(xi1,xi2) -compute_sol(xi1, xi2, rho0_tilde_xi, n, k) + rho_0(xi1, xi2);
    g_bc = @(xi1, xi2) 0;
    
    phi_xi = solve_poisson(f_src, g_bc, n, sigma, eps);
    
    %%% update Phi_xi
    dPhi_xi = -alpha*phi_xi;
    Phi_xi = Phi_xi + dPhi_xi;
    
    
    
    
    Alpha_xi = beta_xi + Phi_xi;
    %%% convexification (-concavization)
    g_bc = @(xi1,xi2) - compute_sol(xi1,xi2, Alpha_xi, n, k) ;
    acc = solve_concav(g_bc, -Alpha_xi, n, sigma, eps, c); %%convexification = -ac
    
    Alpha_xi = -acc;
    Phi_xi = Alpha_xi - beta_xi;
    
    
    
    %%% check energy minimization
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
    
    
    
    cost =  OT_compute_cost(rho_0, Phi_xi, n);
    GRADold = GRAD;
    GRAD    = computeGrad_ddl(n,k,Phi_xi);
    res     = sqrt(sum((GRAD(:,1)-GRADold(:,1)).^2+(GRAD(:,2)-GRADold(:,2)).^2))/sqrt(sum(GRAD(:,1).^2+GRAD(:,2).^2));
    
    l2error = integrate((rho0_proj-rho0_tilde_xi).^2,n)^0.5;
    fprintf('volume=%5e, cost =%5e, res_grad=%5e, l2-error=%5e, counter=%i\n', vol, cost, res, l2error, it); %%exact cost=0.02
    it = it +1;
    
    
end










