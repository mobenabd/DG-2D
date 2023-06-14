function U = solve_poisson(f, g, n, sigma, eps, alpha_xy)
%%%%%%%%%%%%%%%  Solve -Δu = f in Ω = [0,1].^2 %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%          u = g on ∂Ω               %%%%%%%%%%%%%%%%%%%%%
%  !! define k as global varibale before calling this function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k     =  getGlobal_k();

f_t = @(x,y,t) f(x,y); 
g_t = @(x,y,t) g(x,y); 


%%% Solve system
bc_t = 0;  %%=1 if Dirichlet b.c, Neumann otherwise
epsilon = 1e-6;


if ~exist('alpha_xy', 'var')
    [A, Mass] = MatricesSystem(n,sigma,eps,k, bc_t);
    b = SourceBCSystem(n,sigma,eps,k,g_t,f_t,1., bc_t);
else
    [A, Mass] = MatricesSystem(n,sigma,eps,k, bc_t, alpha_xy);
    b = SourceBCSystem(n,sigma,eps,k,g_t,f_t,1., bc_t, alpha_xy);
end

if (bc_t == 0)
    K = A+epsilon*Mass;
    U = K\b;
    
    
    vol = integrate(U,n);
    U = U-vol;

elseif (bc_t==1)
    U = A\b;
end

%fprintf("cond2 = %f\n", cond(K));
%plot_solHeat(n,k,0, U, g_t);


end

