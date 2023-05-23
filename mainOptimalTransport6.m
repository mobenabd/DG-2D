clear; clc; close all;
alpha_list  = [1e-6 1.5e-6 2e-6 4e-6 7e-6 1e-5 2e-5 4e-5 7e-5 1e-4, 1.3e-4, 1.5e-4, 1.75e-4, 2e-4, 2.5e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4, 8e-4, 9e-4, 1e-3,...
    1.2e-3, 1.5e-3, 2e-3, 3e-3, 4e-3, 5e-3, ...
  6e-3, 7e-3, 8e-3, 9e-3, 1e-2, 1.3e-2, 1.5e-2, 1.7e-2,...
  2e-2, 2.3e-2, 2.4e-2, 2.5e-2, 2.6e-2];



quotient = NaN(42,1);

%%%%%%%%%%%%% Solve optimal transport problem %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% test with initializing \delta\Phi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global k
n     =  16;
sigma =  10; eps   =  1;
k    =  3; gdim = n^2*(k+1)^2;
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


beta_old = @(x,y) 0.5*(x.^2 + y.^2);
beta_xi = computeDirectProjection(n,k,beta_old);

%%%%%% initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Phi_xi = zeros(gdim,1);
Phi_xi0 = @(x,y) -0.6 * ( x.^3/3 - x.^2/2 + y.^3/3 - y.^2/2) ;
Phi_xi_0 = computeDirectProjection(n,k,Phi_xi0) - 0.6/6;

for kk=1:42
alpha = alpha_list(kk);

Phi_xi = Phi_xi_0;

Alpha_xi = beta_xi + Phi_xi;

it = 1; E = 0; dE=0;
while (it < itmax)    
    det_xi = computeDet_ddl(n,k,Alpha_xi);     %%det(I + \nabla^2_xi \Phi(\xi))
    
    if(~isempty(find(det_xi<0, 1)))
        fprintf('Warning: det<0, kk=%i\n', kk); %break;
    end
    
    rho0_tilde_xi = OT2_compute_rho0Tilde_xi(n,k, rho_1, Phi_xi).*det_xi;  
    

    
    vol =  integrate(rho0_tilde_xi,n); rho0_tilde_xi = 1/vol * rho0_tilde_xi;
    %OT_plot_sol(n,k, rho0_tilde_xi);


    %%% solve poisson equation
    f_src = @(xi1,xi2) -compute_sol(xi1, xi2, rho0_tilde_xi, n, k) + rho_0(xi1, xi2);
    g_bc = @(xi1, xi2) 0;
    
    phi_xi = solve_poisson(f_src, g_bc, n, sigma, eps);
    
    %%% update Phi_xi
    dPhi_xi = -alpha*phi_xi;
    Phi_xi = Phi_xi + dPhi_xi;
    
    

    
    Alpha_xi = beta_xi + Phi_xi;
    %%% convexification 
      %g_bc = @(x,y) -( beta(x,y) - compute_sol(x, y, lambda_x, n, k) );
      %acc = solve_concav(g_bc, -Alpha_x, n, sigma, eps, c); %%convexification = -ac
      
      %Alpha_x = -acc;
      %lambda_x = beta_x - Alpha_x;
      
    
    
    %%% check energy minimization
    Eold = E;
    dEold = dE;
    [E, dE] = OT2_compute_energy(phi_xi, rho_1, Phi_xi, n);
    if (it>1)
       fprintf('E=%f | Eold=%f | (Eold-E)/(alpha*dE)=%f \n', E, Eold, (Eold-E)/(alpha*dEold));
    end
    
    
    it = it +1; 
    
    
    if (it==3)
        break;
    end
    
end

kk
quotient(kk) = (Eold-E)/(alpha*dEold);
end


