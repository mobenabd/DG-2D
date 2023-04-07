function U = solve_poisson(f, g, n, sigma, eps)
%%%%%%%%%%%%%%%  Solve -Δu = f in Ω = [0,1].^2 %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%          u = g on ∂Ω               %%%%%%%%%%%%%%%%%%%%%
%  !! define k as global varibale before calling this function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k     =  getGlobal_k();

f_t = @(x,y,t) f(x,y); 
g_t = @(x,y,t) g(x,y); 


%%% Solve system
[A,~] = MatricesSystem(n,sigma,eps,k);
b = SourceBCSystem(n,sigma,eps,k,g_t,f_t,1.);

U = A\b;

%plot_solHeat(n,k,0, U, g_t);


end

