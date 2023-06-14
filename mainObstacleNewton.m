clear; clc; close all;
%%%%%%%%%%%%% Solve min(-Δu, u-g)=0 in Ω = [0,1].^2 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  with semismooth newton algorithm      %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% using DG methods with Lgarange basis on quadrature nodes %%%%
%%%%%%%%%%%%%         u = u_final on ∂Ω              %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global k
n     =  8;
global x0 xN hxy
x0=0; xN=1; hxy=(xN-x0)/n; h = hxy;
sigma =  100;
eps   =  -1;
k     =  4;


f    =  @(x,y,t) x.*0;
g = @(x,y,t) obstacle(x,y);
u_final = @(x,y,t) final_sol(x,y);

%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%
[A,Mass] = MatricesSystem(n,sigma,eps,k);
invMass = diag(diag(Mass).^-1);  %%=M^-1 where M is diago,al

b = SourceBCSystem(n,sigma,eps,k,g,f,0.0);
g_projected = diag(invMass).*b;

plot_sol(n,k, g_projected, g, u_final);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gdim = n^2*(k+1)^2;
G = zeros(gdim,1);
dG = zeros(gdim,gdim);
Uold = zeros(gdim,1);
U = g_projected;

res = 1;
t=1;
%U = ones(gdim,1);
b = SourceBCSystem(n,sigma,eps,k,u_final,f,1.);
while (res > 1e-8)
    %compute G[u^k] = -Lap(u^k): use compute_Lap() if a(u^k)=-Lap(u^k)-f and f is not null
    G = A*U - b; 
    %G = diag(invMass).*G; %%% sans changement de varibale
    
    %compute G[u^k] = min(a(u^k), b(u^k)) && Da && Db where %a(u^k)=-Lap(u^k), b(u^k)=u^k-g
    [G,Da,Db] = updateNewton_obstacle(n,k, G, U,g_projected, g);
     
    
    %  Compute jacobian matrix dG
    %dG = Da*invMass*A + Db;  %%% sans changement de varibale
    dG = Da*A + Db;
    
  
    Uold = U;
    
    %U = -gmres(dG,G,[], 1e-12, 1000);
    U = -dG\G;
    U = Uold + U;
    
    
    plot_sol(n,k, U, g, u_final);
    
    
    res = h*sqrt(sum((U-Uold).^2));
    t = t + 1;
    fprintf('res=%f, counter=%i\n', res, t);
  
end
fprintf('\n')

compute_error(U,@(x,y) final_sol(x,y),n,k)

%Lap = compute_Lap(A,Mass,U,b,f,n,k);
%plot_det(n,k, U, Lap);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = obstacle(tx,ty)
x = 4*tx-2;
y = 4*ty-2;


rxy = x.^2 + y.^2;
if (rxy <= 1 )
    out = sqrt(1 - rxy);
else
    out = -1;
end


if ((x==-2 ||  x==2 || y==-2 || y==2))
    out = final_sol(tx,ty);
end
end


function out = final_sol(tx,ty)
x = 4*tx-2;
y = 4*ty-2;

% syms q;
% eqn = q.^2 .*(1-log(q/2)) == 1;
% r = double(solve(eqn,q));

r = 0.69796514822337356642145262626753;

rxy = sqrt(x.^2 + y.^2);

if (rxy <= r )
    out = sqrt(1 - rxy.^2);
else
    out = -r^2*log(rxy/2)/sqrt(1-r^2);
end


end
