clear; clc; close all;
%%%%%%%%%%%%% Solve optimal transport problem %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% test with initializing \delta\Phi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%profile on
global k
n     =  8;
sigma =  10; eps   =  1;
k    =  4; gdim = n^2*(k+1)^2;
c    = 1e4;
alpha = 0.5;%1e-2;
itmax = 200;
b1 = 0.25; b2 = 0.75;
a1 = 0.8; a2 = 1.25;

%%%%%%% densities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TC0 %%%%%%%%%%%%
nu0 = 0.4; nu1 = 0.6;
%rho_0 = @(x,y)  1/0.03141592 * exp(-100*((x-nu0)^2 + (y-nu0)^2));
%rho_1 = @(x,y)  1/0.03141592 * exp(-100*((x-nu1)^2 + (y-nu1)^2)); 

%rho_0 = @(x,y)  1/0.03141592 * exp(-100*((x-nu0)^2 + (y-nu0)^2));
%rho_1 = @(x,y)  1/0.03141592 * exp(-100*((x-nu1)^2 + (y-nu0)^2));


rho_0 = @(x,y)  1/0.531416*(exp(-100*((x-nu0)^2 + (y-nu0)^2)) + 0.5);
rho_1 = @(x,y)  1/0.531416*(exp(-100*((x-nu1)^2 + (y-nu0)^2)) + 0.5);

%%%%%%%%%%%%%%%%%%%%

rho0_proj = computeDirectProjection(n,k,rho_0);
rho1_proj = computeDirectProjection(n,k,rho_1);

beta     = @(x,y) 0.5*((x-0.5).^2 + (y-0.5).^2);% -1/3;
beta_old = @(x,y) 0.5*((x-0.5).^2 + (y-0.5).^2);
beta_xi = computeDirectProjection(n,k,beta_old);% -1/3;

%%%%%% initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Phi_xi_0 = zeros(gdim,1);
%Phi_xi0 = @(x,y) -0.7 * ( x.^3/3 - x.^2/2 + y.^3/3 - y.^2/2 -1/6) ;
%Phi_xi0 = @(x,y) 0.5*(x.^2+y.^2);
%Phi_xi0 = @(x,y) (nu1-nu0)*(x+y)   ;
%Phi_xi_0 = computeDirectProjection(n,k,Phi_xi0) ;

alpha_xy = @(x,y) 0.5*(rho_0(x,y) + rho_1(x,y));
f_source = @(x,y) rho_1(x,y) - rho_0(x,y);
g_bc = @(x, y) 0;

Phi_xi_0 = solve_poisson(f_source, g_bc, n, sigma, eps, alpha_xy);

%Phi_xi = Phi_xi_0 -beta_xi;
Phi_xi = Phi_xi_0;

Alpha_xi = beta_xi + Phi_xi;
%Alpha_xi = Phi_xi_0;

%%%%%% convexification 
g_bc = @(xi1,xi2) -( beta(xi1,xi2) );  %with b.c homogeneous Direchlet 
acc = solve_concav(g_bc, -Alpha_xi, n, sigma, eps, c); %%convexification = -ac
      
Alpha_xi = -acc;
Phi_xi = Alpha_xi - beta_xi;

%%%%%%%%%%%%%%%


%%%%%% restore initial densities %%%
% rho_0 = @(x,y)  1/0.03141592 * exp(-100*((x-nu0)^2 + (y-nu0)^2));
% rho_1 = @(x,y)  1/0.03141592 * exp(-100*((x-nu1)^2 + (y-nu0)^2));
% rho0_proj = computeDirectProjection(n,k,rho_0);
% rho1_proj = computeDirectProjection(n,k,rho_1);
%%%%%%%%%%%%%%%%%%%%%%ùù


it = 1; E = 0; dE=0;
GRAD = zeros(gdim,2);
ll_k = 1;
while (it < itmax)    
    det_xi = computeDet_ddl(n,k,Alpha_xi);     %%det(I + \nabla^2_xi \Phi(\xi))
    
    if(~isempty(find(det_xi<0, 1)))
        fprintf('Warning: det<0\n'); %break;
    end
    
    %rho_1_k = @(x,y) rho_1(x,y) + ll_k;
    %rho0_tilde_xi = OT2_compute_rho0Tilde_xi(n,k, rho_1_k, Phi_xi).*det_xi;  
    rho0_tilde_xi = OT2_compute_rho0Tilde_xi(n,k, rho_1, Phi_xi).*det_xi;    %%modified 
    
    

    vol =  integrate(rho0_tilde_xi,n); 
    %rho0_tilde_xi = 1/vol * rho0_tilde_xi; %%modified 
    OT_plot_sol(n,k, rho0_tilde_xi);


    %%% solve poisson equation
    %rho_0_k = @(x,y) rho_0(x,y) + ll_k;
    %f_src = @(xi1,xi2) -compute_sol(xi1, xi2, rho0_tilde_xi, n, k) + rho_0_k(xi1, xi2); %%modified 
    f_src = @(xi1,xi2) -compute_sol(xi1, xi2, rho0_tilde_xi, n, k) + rho_0(xi1, xi2);
    g_bc = @(xi1, xi2) 0;
    
    phi_xi = solve_poisson(f_src, g_bc, n, sigma, eps);
    
    %%% update Phi_xi
    dPhi_xi = -alpha*phi_xi;
    Phi_xi = Phi_xi + dPhi_xi;
    %lambda_x = lambda_x - integrate(lambda_x,n);
    

    
    Alpha_xi = beta_xi + Phi_xi;
    %%% convexification (-concavization)
      %g_bc = @(xi1,xi2) -( beta(xi1,xi2) + compute_sol(xi1,xi2, Phi_xi, n, k) );
      g_bc = @(xi1,xi2) -( beta(xi1,xi2) );  %with b.c homogeneous Direchlet 
      acc = solve_concav(g_bc, -Alpha_xi, n, sigma, eps, c); %%convexification = -ac
      
      Alpha_xi = -acc;
      Phi_xi = Alpha_xi - beta_xi;
      
    
    
    %%% check energy minimization
    Eold = E;
    dEold = dE;
    %[E, dE] = OT2_compute_energy(phi_xi, rho_1_k, Phi_xi, n);
    [E, dE] = OT2_compute_energy(phi_xi, rho_1, Phi_xi, n); %%modified 
    if (it>1)
       fprintf('E=%f | dE=%f | alpha=%f| (Eold-E)/(alpha*dE)=%f \n', E, dE, alpha, (Eold-E)/(alpha*dEold));
       %%% update step size
         if ( E-Eold > -b1*alpha*dEold )
              alpha = a1*alpha;
         elseif ( E-Eold < -b2*alpha*dEold )
              alpha = a2*alpha;
         end
    end
    
    

    %cost =  OT_compute_cost(rho_0_k, Phi_xi, n);
    cost =  OT_compute_cost(rho_0, Phi_xi, n); %%modified 
    GRADold = GRAD; 
    GRAD    = computeGrad_ddl(n,k,Phi_xi);
    res     = sqrt(sum((GRAD(:,1)-GRADold(:,1)).^2+(GRAD(:,2)-GRADold(:,2)).^2))/sqrt(sum(GRAD(:,1).^2+GRAD(:,2).^2));
    
    fprintf('volume=%f, cost =%f, res_grad=%f, counter=%i\n', vol, cost, res, it); %%exact cost=0.02
    it = it +1; 
    
    
    
    %ll_k = 0.9*ll_k;
end




