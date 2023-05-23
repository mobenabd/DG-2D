clear; clc; close all;
%%%%%%%%%%%%% Solve optimal transport problem %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global k
n     =  6;
sigma =  10; eps   =  1;
k    =  4; gdim = n^2*(k+1)^2;
c    = 1e4;
alpha = 0.03;
itmax = 200;
b1 = 0.9; b2 = 0.2;
a1 = 1.25; a2 = 0.8;

%%%%%%% densities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TC1 %%%%%%%%%%%%
%rho0 = @(x,y) 3.97887 * exp(-12.5*(x.^2 + (y+0.1).^2));
%rho1 = @(x,y) 4.59441 * exp(-16.66666*(x.^2 + x.*(y-0.1) + (y-0.1).^2));
%%%%%%%%%%%%%%%%%%%%

%%% TC2 %%%%%%%%%%%%
% rho0 = @(x,y) 1.98944 * (exp(-12.5*((x-0.5).^2 + (y+0.3).^2)) ...
%     + exp(-12.5*((x+0.5).^2 + (y+0.3).^2)));
% rho1 = @(x,y) 3.97887*exp(-12.5*(x.^2 + y.^2));
%%%%%%%%%%%%%%%%%%%

%rho_0 = @(x,y) 16*rho0(4*x-2, 4*y-2);   %%defined on (0,1)^2
%rho_1 = @(x,y) 16*rho1(4*x-2, 4*y-2);   %%defined on (0,1)^2


%%% TC0 %%%%%%%%%%%%
nu0 = 0.4; nu1 = 0.42;
rho_0 = @(x,y)  1/0.03141592 * exp(-100*((x-nu0)^2 + (y-0.4)^2));
rho_1 = @(x,y)  1/0.03141592 * exp(-100*((x-nu1)^2 + (y-0.42)^2)); 

% q_func = @(z) (-1/(8*pi)*z^2 + 1/(256*pi^3) + 1/(32*pi))*cos(8*pi*z) + 1/(32*pi^2)*z*sin(8*pi*z);
% q_prim = @(z) (228155022448185*sin(8*pi*z))/72057594037927936 - (5734161139222659*z*cos(8*pi*z))/72057594037927936 + 8*pi*sin(8*pi*z)*((5734161139222659*z^2)/144115188075855872 - 371634248860825663/36893488147419103232) + (228155022448185*z*pi*cos(8*pi*z))/9007199254740992;
% q_sec = @(z) (228155022448185*pi*cos(8*pi*z))/4503599627370496 - (5734161139222659*cos(8*pi*z))/72057594037927936 + 64*pi^2*cos(8*pi*z)*((5734161139222659*z^2)/144115188075855872 - 371634248860825663/36893488147419103232) - (228155022448185*z*pi^2*sin(8*pi*z))/1125899906842624 + (5734161139222659*z*pi*sin(8*pi*z))/4503599627370496;
% 
% 
% rho_0 = @(x,y) 1 + 4*(q_sec(x)*q_func(y) + q_func(x)*q_sec(y)) + 16*(q_func(x)*q_func(y)*q_sec(x)*q_sec(y) ...
%         - q_prim(x)^2*q_prim(y)^2  );
% rho_1 = @(x,y) 1;
%%%%%%%%%%%%%%%%%%%%

rho0_proj = computeDirectProjection(n,k,rho_0);
rho1_proj = computeDirectProjection(n,k,rho_1);


%%%%%% initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lambda0 = @(x,y) 0.3*(x.^2 + y.^2);
%lambda0 = @(x,y) 0;
%lambda_x = computeDirectProjection(n,k,lambda0);
lambda_x = zeros(gdim,1);


beta_old = @(x,y) 0.5*(x.^2 + y.^2);
beta     = @(x,y) 0.5*(x.^2 + y.^2) - 1/3;
beta_x = computeDirectProjection(n,k,beta_old) - 1/3;

%beta_old =  @(x,y)  0.5*((x-0.5).^2 + (y-0.5).^2);
%beta =  @(x,y) 0.5*((x-0.5).^2 + (y-0.5).^2) - 2/3 *(0.5^3);
%beta_x = computeDirectProjection(n,k,beta_old) - 2/3 *(0.5^3); 

Alpha_x = beta_x - lambda_x;

it = 1; E = 0;
while (it < itmax)    
    det_x = computeDet_ddl(n,k,Alpha_x);     %%det(I - \nabla^2_x \lambda(\x))
    
    if(~isempty(find(det_x<0, 1)))
        fprintf('Warning: det<0\n'); return;
    end
    
    rho0_tilde_x = rho1_proj./det_x;    
    rho0_tilde_xi = OT_compute_rho0Tilde_xi(n,k, rho0_tilde_x, lambda_x);   %M = OT_Xxi(n,k,lambda_x);  %%M = (x,y, xi_1, xi_2)
    
    
    vol =  integrate(rho0_tilde_xi,n);
     
    rho0_tilde_xi = 1/vol * rho0_tilde_xi;
    OT_plot_sol(n,k, rho0_tilde_xi);


    %%% solve poisson equation
    f_src = @(xi1,xi2) -compute_sol(xi1, xi2, rho0_tilde_xi, n, k) + rho_0(xi1, xi2);
    g_bc = @(xi1, xi2) 0;
    
    phi_xi = solve_poisson(f_src, g_bc, n, sigma, eps);
    
    %%% update lambda
    dlambda_x = -alpha*OT_compute_Dlambda(n,phi_xi, lambda_x);
    lambda_x = lambda_x + dlambda_x;
    lambda_x = lambda_x - integrate(lambda_x,n);
    

    
    Alpha_x = beta_x - lambda_x;
    %%% convexification 
      g_bc = @(x,y) -( beta(x,y) - compute_sol(x, y, lambda_x, n, k) );
      acc = solve_concav(g_bc, -Alpha_x, n, sigma, eps, c); %%convexification = -ac
      
      Alpha_x = -acc;
      lambda_x = beta_x - Alpha_x;
      
    
    
    %%% check energy minimization
    Eold = E;
    [E, dE] = OT_compute_energy(phi_xi, rho_1, lambda_x, n);
    if (it>1)
       fprintf('E=%f, dE=%f, alpha=%f | (Eold-E) - alpha*dE=%f\n', E, dE, alpha, Eold-E- alpha*dE);
       %%% update step size
       if ( E-Eold > -b1*alpha*dE && alpha > 0.01 )
            alpha = a2*alpha;
       elseif ( E-Eold < -b2*alpha*dE )
            alpha = a1*alpha;
       end
    end
    
    

    %cost = abs(OT_compute_cost(rho_0, lambda_x, Alpha_x, n) - 2*(nu0-nu1)^2);
    cost =  OT_compute_cost(rho_0, lambda_x, Alpha_x, n);
    fprintf('volume=%f, cost error=%f, counter=%i\n', vol, cost, it);
    it = it +1 ;
    
    

end



