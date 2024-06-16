function U = solve_poisson2(f, g, n, sigma, eps, Ldecomp, Udecomp, Mass)
%%%%%%%%%%%%%%%  Solve -Δu = f in Ω = [0,1].^2 %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  the Laplacian matrix is A = Ldecomp * Udecomp
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
    b = SourceBCSystem(n,sigma,eps,k,g_t,f_t,1., bc_t);
else
    b = SourceBCSystem(n,sigma,eps,k,g_t,f_t,1., bc_t, alpha_xy);
end

if (bc_t == 0)    
    Ldecomp2 = Ldecomp+epsilon*Mass;
    Udecomp2 = Udecomp+epsilon*Mass;
    %Y = Ldecomp2\b;
    %U = Udecomp2\Y;
    opts = struct('LT', {true});
    Y = linsolve(Ldecomp2,b,opts);
    opts = struct('UT', {true});
    U = linsolve(Udecomp2,Y,opts);
    
    vol = integrate(U,n);
    U = U-vol;

elseif (bc_t==1)
    Y = Ldecomp\b;
    U = Udecomp\Y;
end


%fprintf("cond2 = %f\n", cond(K));
%plot_solHeat(n,k,0, U, g_t);


end

