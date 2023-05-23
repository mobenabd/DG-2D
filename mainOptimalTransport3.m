clear; clc; close all;
alpha_list  = [1e-6 1.5e-6 2e-6 4e-6 7e-6 1e-5 2e-5 4e-5 7e-5 1e-4, 1.3e-4, 1.5e-4, 1.75e-4, 2e-4, 2.5e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4, 8e-4, 9e-4, 1e-3,...
    1.2e-3, 1.5e-3, 2e-3, 3e-3, 4e-3, 5e-3, ...
  6e-3, 7e-3, 8e-3, 9e-3, 1e-2, 1.3e-2, 1.5e-2, 1.7e-2,...
  2e-2, 2.3e-2, 2.4e-2, 2.5e-2, 2.6e-2];

%alpha_list  = [1e-7 3e-7 5e-7 6e-7 8e-7 1e-6];


quotient = NaN(33,1);

%%%%%%%%%%%%% Solve optimal transport problem %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global k
n     =  10;
sigma =  10; eps   =  1;
k    =  4; gdim = n^2*(k+1)^2;
c    = 1e4;
alpha = 0.01;
itmax = 200;
b1 = 0.9; b2 = 0.2;
a1 = 1.25; a2 = 0.8;

%%%%%%% densities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TC0 %%%%%%%%%%%%
nu0 = 0.4; nu1 = 0.42;
rho_0 = @(x,y)  1/0.03141592 * exp(-100*((x-nu0)^2 + (y-nu0)^2));
rho_1 = @(x,y)  1/0.03141592 * exp(-100*((x-nu1)^2 + (y-nu1)^2)); 

%%%%%%%%%%%%%%%%%%%%

rho0_proj = computeDirectProjection(n,k,rho_0);
rho1_proj = computeDirectProjection(n,k,rho_1);

for kk=1:42
alpha = alpha_list(kk);
%%%%%% initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lambda_x = zeros(gdim,1);
lambda0 = @(x,y) 0.6 * ( x.^3/3 - x.^2/2 + y.^3/3 - y.^2/2) ;
lambda_x = computeDirectProjection(n,k,lambda0) + 0.6/6;


beta_old = @(x,y) 0.5*(x.^2 + y.^2);
beta     = @(x,y) 0.5*(x.^2 + y.^2) - 1/3;
beta_x = computeDirectProjection(n,k,beta_old) - 1/3;

Alpha_x = beta_x - lambda_x;

X0old = getPhysicalNodes(n);
it = 1; E = 0; dE=0;
while (it < itmax)    
    det_x = computeDet_ddl(n,k,Alpha_x);     %%det(I - \nabla^2_x \lambda(\x))
    
    if(~isempty(find(det_x<0, 1)))
        fprintf('Warning: det<0\n'); break;
    end
    
    rho0_tilde_x = rho1_proj./det_x;   
    [rho0_tilde_xi, X0new] = OT_compute_rho0Tilde_xi(n,k, rho0_tilde_x, lambda_x, X0old); 
    X0old = X0new;
    
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
    %%% convexification 
      %g_bc = @(x,y) -( beta(x,y) - compute_sol(x, y, lambda_x, n, k) );
      %acc = solve_concav(g_bc, -Alpha_x, n, sigma, eps, c); %%convexification = -ac
      
      %Alpha_x = -acc;
      %lambda_x = beta_x - Alpha_x;
      
    
    
    %%% check energy minimization
    Eold = E;
    dEold = dE;
    [E, dE] = OT_compute_energy(phi_xi, rho_1, lambda_x, n);
    if (it>1)
       fprintf('E=%f | Eold=%f | (Eold-E)/(alpha*dE)=%f \n', E, Eold, (Eold-E)/(alpha*dE));
       %%% update step size
%        if ( E-Eold > -b1*alpha*dE && alpha > 0.01 )
%             alpha = a2*alpha;
%        elseif ( E-Eold < -b2*alpha*dE )
%             alpha = a1*alpha;
%        end
    end
    
    

    %cost = abs(OT_compute_cost(rho_0, lambda_x, Alpha_x, n) - 2*(nu0-nu1)^2);
    %cost =  OT_compute_cost(rho_0, lambda_x, Alpha_x, n);
    %fprintf('volume=%f, cost error=%f, counter=%i\n', vol, cost, it);
    it = it +1 ;
    
    
    if (it==3)
        break;
    end
    
end

kk
quotient(kk) = (Eold-E)/(alpha*dE);
end


