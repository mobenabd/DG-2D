clear; clc; close all;

%%%%%%%%%%%%% Solve optimal transport problem %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global k
n     =  8;
sigma =  10; eps   =  1;
k    =  4; gdim = n^2*(k+1)^2;
c    = 1e4;
alpha = 0.01;
itmax = 200;


%%%%%%% densities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TC0 %%%%%%%%%%%%
nu0 = 0.4; nu1 = 0.42;
rho_0 = @(x,y)  1/0.03141592 * exp(-100*((x-nu0)^2 + (y-nu0)^2));
rho_1 = @(x,y)  1/0.03141592 * exp(-100*((x-nu1)^2 + (y-nu1)^2)); 

%%%%%%%%%%%%%%%%%%%%

rho0_proj = computeDirectProjection(n,k,rho_0);
rho1_proj = computeDirectProjection(n,k,rho_1);


%%%%%% initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lambda_x = zeros(gdim,1);
lambda0 = @(x,y) 0.6 * ( x.^3/3 - x.^2/2 + y.^3/3 - y.^2/2) ;
lambda_x = computeDirectProjection(n,k,lambda0) + 0.6/6;
lambda_x0 = lambda_x;

beta     = @(x,y) 0.5*(x.^2 + y.^2) ;
beta_x = computeDirectProjection(n,k,beta);

Alpha_x = beta_x - lambda_x;

X0old = getPhysicalNodes(n);
it = 1; 
while (it < itmax)    
    det_x = computeDet_ddl(n,k, Alpha_x);     %%det(I - \nabla^2_x \lambda(\x))
    
    if(~isempty(find(det_x<0, 1)))
        fprintf('Warning: det<0\n'); break;
    end
    
    rho0_tilde_x = rho1_proj./det_x;   
    
    if (it==2)
        break;
        %%% continue to compute rho0_tilde_xi_1
    end
    %use \lambda and X(xi) exact for it=1
    [rho0_tilde_xi, X0new] = OT_compute_rho0Tilde_xi(n,k, rho0_tilde_x, lambda_x, X0old); 
    X0old = X0new;
    
    rho0_tilde_xi(rho0_tilde_xi<0) = 0; % positive density
    
    vol =  integrate(rho0_tilde_xi,n);
    rho0_tilde_xi = 1/vol * rho0_tilde_xi;
    %OT_plot_sol(n,k, rho0_tilde_xi);

    
    

    %%% solve poisson equation
    f_src = @(xi1,xi2) -compute_sol(xi1, xi2, rho0_tilde_xi, n, k) + rho_0(xi1, xi2);
    g_bc = @(xi1, xi2) 0;
    
    phi_xi = solve_poisson(f_src, g_bc, n, sigma, eps);
    
    %%% update lambda
    dlambda_x = -alpha*OT_compute_Dlambda(n,phi_xi, lambda_x);
    lambda_x = lambda_x + dlambda_x;
    %lambda_x = lambda_x - integrate(lambda_x,n);
    

    
    Alpha_x = beta_x - lambda_x;

            
    
    rho0_tilde_xi_0 = rho0_tilde_xi; 
    det_x_0 = det_x;
    
    
    it = it +1 ;
     
end


%%% compute rho0_tilde_xi_1
[rho0_tilde_xi_1, ~] = OT_compute_rho0Tilde_xi_TEST(n,k, rho0_tilde_x, lambda_x, X0old);
rho0_tilde_xi_1(rho0_tilde_xi_1<0) = 0; % positive density
vol =  integrate(rho0_tilde_xi_1,n);


rho0_tilde_xi_1 = 1/vol * rho0_tilde_xi_1;


%%% calcul de \delta\Tilde{\rho_0}
d_rho0Tilde_x = OT_compute_d_rho0Tilde_x(n,k, rho1_proj, dlambda_x)./det_x_0;
[d_rho0Tilde_xi, ~] = OT_compute_rho0Tilde_xi(n,k, d_rho0Tilde_x, lambda_x0, X0old); 

U1 = rho0_tilde_xi_1-rho0_tilde_xi_0;
U2 = d_rho0Tilde_xi;

normL2 = sqrt(sum((U1-U2).^2))/sqrt(sum(U1.^2));

normLinf = max(abs(U1-U2));


fprintf("L2-norm=%f | Linf-norm=%f\n", normL2, normLinf);


OT_plot_sol(n,k, U1);
hold on
OT_plot_sol(n,k, d_rho0Tilde_x);


