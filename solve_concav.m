function U = solve_concav(gxy,g_projected, n, sigma, eps, c)
%%%%%%%%%%%%% Solve min(-det(u)-f, u-g)=0 in Ω = [0,1].^2 %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  with linearized obstacle problem %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  + Monge-Ampere Equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% using DG methods with Lagrange basis on quadrature nodes %%%%
%%%%%%%%%%%%%         u = u_bc (concave) on ∂Ω              %%%%%%%%%%%%%%%
% !! define k as global varibale before calling this function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 1e-4;
k     =  getGlobal_k();
gdim = n^2*(k+1)^2;

U = g_projected;


DET = computeDet_ddl(n,k,U);    %%compute determinant

%if (isempty(find(DET<0.5, 1)))   %%function is already concave
%    return;
%end
%return

Lap = computeLap2_ddl(n,k,U);   %%compute laplacian



g = @(x,y,t) gxy(x,y);

%beta =  @(x,y,t) -0.5*((x-0.5).^2 + (y-0.5).^2);% + 0.1*x + 0.1*y -0.2;  %%Dirichlet boundary condition

beta =  @(x,y,t) gxy(x,y);  %%Dirichlet boundary condition

%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%
neumm = 1;
[A,Mass] = MatricesSystem(n,sigma,eps,k, neumm);

if (neumm == 0)
    A = A + 1e-5*Mass;
end

%f_projected = 0.1*ones(gdim,1);
%f_projected = zeros(gdim,1) ;


%f_projected = 0.5*ones(gdim,1) ; 
f_projected = 0.001*ones(gdim,1) ; 


res = 1;
t=1;
Eig_m = 0;
while (res > tol)
    Mk =  ACM_GlobalMass(Eig_m, DET, n,k,c,U,Lap,g_projected, f_projected);
    RHS = ACM_RHS(Eig_m, DET, n,k,c,U,Lap,sigma,eps,beta,g_projected, f_projected, neumm);
    
    Uold = U;
    U = (A+Mk)\RHS;


    
    %plot_sol(n,k, U, g);

    
    res = sqrt(sum((U-Uold).^2))/sqrt(sum(U.^2));
    fprintf('Obstacle problem: res=%f, itr=%i, neg=%f\n', res, t,length(find(DET<0)));
    
    
    t = t + 1;

    Lap = computeLap2_ddl(n,k,U);   %%compute laplacian u^{k-1}
    DET = computeDet_ddl(n,k,U);    %%compute determinant u^{k-1}
    
    
    if (res < tol && Eig_m == 0 && ~isempty(find(DET<0, 1)) && neumm==1)
        Eig_m = 1;
        c = 1;
        res = 1;
        fprintf('Obstacle problem: Regularization step\n');
    end
    
    if (Eig_m == 1 && isempty(find(DET<0, 1)))
        break;
    end
end


if (neumm == 0)
    vol = integrate(U,n);
    U = U-vol;
end

%plot_sol(n,k, U, g);

%fprintf('Obstacle problem: res=%f, nbr itr=%i\n', res, t-1);


%plot_sol(n,k, U, g);
end

