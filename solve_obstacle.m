function U = solve_obstacle(fxy, gxy,g_projected, n, sigma, eps, c)
%%%%%%%%%%%%% Solve min(-Δu-f, u-g)=0 in Ω = [0,1].^2 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  with linearized obstacle problem %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% using DG methods with Lgarange basis on quadrature nodes %%%%
%%%%%%%%%%%%%         u = u_bc on ∂Ω              %%%%%%%%%%%%%%%%%%%%%%%%
%  !! define k as global varibale before calling this function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 1e-3;
k     =  getGlobal_k();
gdim = n^2*(k+1)^2;


f    =  @(x,y,t) fxy(x,y);
g = @(x,y,t) gxy(x,y);
%u_final = @(x,y,t) final_sol(x,y);

%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%
[A,~] = MatricesSystem(n,sigma,eps,k);


%g_projected = computeDirectProjection(n,k,g);

U = g_projected;


%b = SourceBCSystem(n,sigma,eps,k, f,f,0.0);
%f_projected = Mass\b;
%plot_sol(n,k, f_projected);

f_projected = zeros(gdim,1) ;
 
DET = NaN(gdim,1);
Eig_m = DET;

res = 1;
t=1;
while (res > tol )
    Lap = computeLap2_ddl(n,k,U);   %%compute determinant
    %DET = computeDet_ddl(n,k,U);    %%compute determinant
    %Eig_m = computeDet_ddl(n,k,U,1);  %%compute largest eigenvalue
    Mk =  ACM_GlobalMass(Eig_m, DET, n,k,c,U,Lap,g_projected, f_projected);
    RHS = ACM_RHS(Eig_m, DET, n,k,c,U,Lap,sigma,eps,g,g_projected, f_projected);
    
    
    Uold = U;
    U = (A+Mk)\RHS;

        

    
    res = sqrt(sum((U-Uold).^2))/sqrt(sum(U.^2));
    
    t = t + 1;
    
end
fprintf('Obstacle problem: res=%f, nbr itr=%i\n', res, t-1);


%plot_sol(n,k, U, g);
end

