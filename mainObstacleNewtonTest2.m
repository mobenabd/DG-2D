clear; clc; close all;
%%%%%%%%%%%%% Solve min(-Δu, u-g)=0 in Ω = [0,1].^2 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  with semismooth newton algorithm      %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% using DG methods with Lgarange basis on quadrature nodes %%%%
%%%%%%%%%%%%%         u = u_final on ∂Ω              %%%%%%%%%%%%%%%%%%%%%%
% This main is used for tests without exact solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global k
n     =  10;
sigma =  100;
eps   =  -1;
k     =  4;
global x0 xN hxy
x0=0; xN=1; hxy=(xN-x0)/n; h = hxy;


f    =  @(x,y,t) x.*0  ;
%f    =  @(x,y,t) source(x,y);
g = @(x,y,t) obstacle(x,y);
u_final = @(x,y,t) final_sol(x,y);

%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%
[A,Mass] = MatricesSystem(n,sigma,eps,k);
invMass = diag(diag(Mass).^-1);  %%=M^-1 where M is diago,al

b = SourceBCSystem(n,sigma,eps,k, g,f,0.0);
g_projected = diag(invMass).*b;

% cond(A)
% cond(Mass)
% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gdim = n^2*(k+1)^2;
G = zeros(gdim,1);
dG = zeros(gdim,gdim);
Uold = zeros(gdim,1);
U = g_projected;

%plot_sol2(n,k, U);

res = 1;
t=1;
%U = ones(gdim,1);
b = SourceBCSystem(n,sigma,eps,k,u_final,f,1.);
while (res > 1e-8)
%%%%%compute G[u^k] = -Lap(u^k) or a(u^k)=-Lap(u^k)-f  %%%%%%%%%
    %use compute_Lap() to copmpute Laplacian
    G = A*U - b; 
    G = diag(invMass).*G; %%% sans changement de varibale
    
    %G = -computeLap2_ddl(n,k,U);  %% using stronf laplacian

    
    %Lap = compute_Lap(A,Mass,U,b,f,n,k);
    %G = computeDet_ddl(n,k,Lap,U);
    %plot_sol(n,k, LAP)
    
%%%%% compute G[u^k] = min(a(u^k), b(u^k)) && Da && Db where %a(u^k)=-Lap(u^k), b(u^k)=u^k-g
    [G,Da,Db] = updateNewton_obstacle(n,k, G, U,g_projected, g);
     
%%%%% Compute jacobian matrix dG
    dG = Da*invMass*A + Db;  %%% sans changement de varibale
    %dG = Da*A + Db;         %%% avec changement de varibale
    
    Uold = U;
    
    %U = -gmres(dG,G,[], 1e-12, 1000);
    U = -dG\G;
    U = Uold + U;
    
    
    plot_sol(n,k, U, g)%, u_final);
    
    
%%%%%%% plot determinant %%%%%%%%%
    %Lap = compute_Lap(A,Mass,U,b,f,n,k);
    %plot_det(n,k, U, Lap, @(x,y,t) determiant_exct(x,y));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %plot_sol2(n,k, Lap);
    %return
    
    res = h*sqrt(sum((U-Uold).^2));
    t = t + 1;
    fprintf('res=%f, counter=%i\n', res, t);
end
fprintf('\n')



%Lap = compute_Lap(A,Mass,U,b,f,n,k);

%plot_sol(n,k, Lap)%, u_final);

%plot_det(n,k, U, Lap); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = obstacle(tx,ty)

%out = sin(pi*tx).*cos(2*pi*tx).^2  + sin(pi*ty).*cos(2*pi*ty).^2 ;

% if (tx==0 || tx==1 || ty==0 || ty==1)
%     out = final_sol(tx,ty);
% end

out = sin(2*pi*tx)*sin(2*pi*ty);

end



function out = final_sol(tx,ty)

%out = sin(pi*tx)  + sin(pi*ty);

out = 0.;
end


function out = source(tx,ty)

out = pi^2*(sin(pi*tx) + sin(pi*ty));

end



