clear; clc; close all;
%%%%%%%%%%%%% Solve optimal transport problem %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global k
n     =  8;
sigma =  10; eps   =  1;
k    =  4; gdim = n^2*(k+1)^2;
c    = 1e3;
alpha = 0.01;

%%%%%%% densities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho0 = @(x,y) 1.98944 * (exp(-12.5*((x-0.5).^2 + (y+0.3).^2)) ...
    + exp(-12.5*((x+0.5).^2 + (y+0.3).^2)));
rho1 = @(x,y) 3.97887*exp(-12.5*(x.^2 + y.^2));

rho_0 = @(x,y) 16*rho0(4*x-2, 4*y-2);   %%defined on (0,1)^2
rho_1 = @(x,y) 16*rho1(4*x-2, 4*y-2);   %%defined on (0,1)^2

rho0_proj = computeDirectProjection(n,k,rho_0);
rho1_proj = computeDirectProjection(n,k,rho_1);

%%%%%% initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lambda0 = @(x,y) 0.3*(x.^2 + y.^2);
lambda0 = @(x,y) 0;
lambda_x = computeDirectProjection(n,k,lambda0);

beta =  @(x,y) 0.5*(x.^2 + y.^2);
beta_x = computeDirectProjection(n,k,beta); %% 0.5(x^2 + y^2) projected

Alpha_x = beta_x - lambda_x;

it = 0;
while (it < 20)
    %det_x = OT_computeDet(n,k,lambda_x, 1);     %%det(I - \nabla^2_x \lambda(\x))
    det_x = computeDet_ddl(n,k,Alpha_x);     %%det(I - \nabla^2_x \lambda(\x))
    
    rho0_tilde_x = rho1_proj./det_x;
    M = OT_Xxi(n,k,lambda_x);  %%M = (x,y, xi_1, xi_2)
    rho0_tilde_xi = OT_compute_rho0Tilde_xi(n,k, rho0_tilde_x, M);
    
    
    OT_plot_sol(n,k, det_x);
    if ( integrate(rho0_tilde_xi,n) < 0.5)
        return;
    end
       
    %%% solve poisson equation
    f_src = @(xi1,xi2) -compute_sol(xi1, xi2, rho0_tilde_xi, n, k) + rho_0(xi1, xi2);
    g_bc = @(xi1, xi2) 0;
    
    phi_xi = solve_poisson(f_src, g_bc, n, sigma, eps);
    
    %%% update lambda
    dlambda_x = -alpha*OT_compute_Dlambda(n,phi_xi, lambda_x);
    
    lambda_x = lambda_x + dlambda_x;
    

    %%% convexification 
    Alpha_x = beta_x - lambda_x;
    g_bc = @(x,y) -(0.5*(x.^2 + y.^2) - compute_sol(x, y, lambda_x, n, k));
    f_src = @(x, y) 0;
    %ac = -computeDirectProjection(n,k,alpha_x);
    %OT_plot_sol(n,k, ac);
    acc = solve_obstacle(f_src, g_bc, -Alpha_x, n, sigma, eps, c); %%convexification = -ac
    
    Alpha_x = -acc;
    lambda_x = beta_x - Alpha_x;
    
   
    
    it = it +1 ;
    fprintf('volume=%f, counter=%i\n', integrate(rho0_tilde_xi,n), it);
    %OT_plot_sol(n,k, -acc);
    %OT_plot_sol(n,k, rho0_tilde_xi);
end

%OT_plot_sol(n,k, rho0_tilde_xi);

