clear; clc; close all;
%%%%%%%%%%%%% Solve optimal transport problem %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global k
n     =  10;
sigma =  10; eps   =  1;
k    =  4; gdim = n^2*(k+1)^2;
c    = 1e4;
alpha = 3e-3;
itmax = 200;
b1 = 1; b2 = 0.2;
a1 = 1.25; a2 = 0.8;

%%%%%%% densities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TC0 %%%%%%%%%%%%
ss = 1/200;
nu0 = 0.4; nu1 = 0.6;
rho_0 = @(x,y)  1/(2*pi*ss) * exp(-0.5*(ss^-1)*((x-nu0)^2 + (y-nu0)^2));
rho_1 = @(x,y)  1/(2*pi*ss) * exp(-0.5*(ss^-1)*((x-nu1)^2 + (y-nu1)^2)); 


rho0_proj = computeDirectProjection(n,k,rho_0);
rho1_proj = computeDirectProjection(n,k,rho_1);

%  OT_plot_sol(n,k, rho0_proj+rho1_proj);
%  integrate(rho0_proj,n)
%  integrate(rho1_proj,n)
% return

%%%%%% initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lambda0 = @(x,y) 0.1*(x + y - 1);
%lambda_x = computeDirectProjection(n,k,lambda0);
lambda_x = zeros(gdim,1);


beta_old = @(x,y) 0.5*(x.^2 + y.^2);
beta     = @(x,y) 0.5*(x.^2 + y.^2) - 1/3;
beta_x   = computeDirectProjection(n,k,beta_old) - 1/3;


Alpha_x = beta_x - lambda_x;

X0old = getPhysicalNodes(n);
it = 1; E = 0; GRAD = zeros(gdim,2);
while (it < itmax)    
    det_x = computeDet_ddl(n,k,Alpha_x);     %%det(I - \nabla^2_x \lambda(\x))
    
    if(~isempty(find(det_x<0, 1)))
        fprintf('Warning: det<0\n'); return;
    end
    
    rho0_tilde_x = rho1_proj./det_x;   
    [rho0_tilde_xi, X0new] = OT_compute_rho0Tilde_xi(n,k, rho0_tilde_x, lambda_x, X0old); 
    X0old = X0new;
    
    vol =  integrate(rho0_tilde_xi,n);
    
    rho0_tilde_xi = (1/vol) * rho0_tilde_xi;
    OT_plot_sol(n,k, rho0_tilde_xi);


    %if (it==3)
    %   break;
    %end
    
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
      g_bc = @(x,y) -( beta(x,y) - compute_sol(x, y, lambda_x, n, k) );
      acc = solve_concav(g_bc, -Alpha_x, n, sigma, eps, c); %%convexification = -ac
      
      Alpha_x = -acc;
      lambda_x = beta_x - Alpha_x;
      
    
    
    %%% check energy minimization
    Eold = E;
    [E, dE] = OT_compute_energy(phi_xi, rho_1, lambda_x, n);
    quotient = (Eold-E)/(alpha*dE);
    if (it>1)
       fprintf('E=%f, Eold=%f, alpha=%f | (Eold-E)/(alpha*dE)=%f\n', E, Eold, alpha, quotient);
       %%% update step size
%        if ( quotient > b1 && alpha > 0.01 )
%             alpha = a2*alpha;
%        elseif ( quotient < b2 )
%             alpha = a1*alpha;
%        end
    end
    
    

    %cost = abs(OT_compute_cost(rho_0, lambda_x, Alpha_x, n) - 2*(nu0-nu1)^2);
    cost    =  OT_compute_cost(rho_0, lambda_x, Alpha_x, n);
    GRADold = GRAD; 
    GRAD    = computeGrad_ddl(n,k,lambda_x);
    res     = sqrt(sum((GRAD(:,1)-GRADold(:,1)).^2+(GRAD(:,2)-GRADold(:,2)).^2))/sqrt(sum(GRAD(:,1).^2+GRAD(:,2).^2));
    fprintf('volume=%f, cost=%f, res_grad=%f, counter=%i\n', vol, cost, res, it);
    it = it +1 ;
    
    

end



