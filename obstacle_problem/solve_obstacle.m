function U = solve_obstacle(gxy,g_projected, n, sigma, eps, c, f)
%%%%%%%%%%%%% Solve min(-det(u)-f, u-g)=0 in Ω = [x0,xN].^2 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%  with linearized obstacle problem %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% using DG methods with Lagrange basis on quadrature nodes %%%%
%%%%%%%%%%%%%         u = gxy (concave) on ∂Ω              %%%%%%%%%%%%%%%
% !! define k as global varibale before calling this function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Eig_m = 0; %% no concavization
tol   = 1e-6;
k     =  getGlobal_k();
gdim  = n^2*(k+1)^2;
gxyt =  @(x,y,t) gxy(x,y);  %%Dirichlet boundary condition

%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%
[A,~] = MatricesSystem(n,sigma,eps,k);

U = g_projected;


DET = NaN(gdim,1);              %%compute determinant
Lap = computeLap2_ddl(n,k,U);   %%compute laplacian
f_projected = f*ones(gdim,1) ; 


res = 1; t=1;
while (res > tol)
    Mk =  ACM_GlobalMass(Eig_m, DET, n,k,c,U,Lap,g_projected, f_projected);
    RHS = ACM_RHS(Eig_m, DET, n,k,c,U,Lap,sigma,eps,gxyt,g_projected, f_projected);
    
    Uold = U;
    U = (A+Mk)\RHS;

    plot_sol(n,k, U);
    
    Lap = computeLap2_ddl(n,k,U);   %%compute laplacian u^{k-1}  
    res = sqrt(sum((U-Uold).^2))/sqrt(sum(U.^2));
    fprintf('Obstacle problem: res=%.15f, itr=%i\n', res, t);
    
    
    t = t + 1;
end


%plot_sol3(n,k, U);

%fprintf('Obstacle problem: res=%f, nbr itr=%i\n', res, t-1);

end

