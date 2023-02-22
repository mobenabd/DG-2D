clear; clc; close all;
%%%%%%%%%%%%% Solve min(-Δu, u-g)=0 in Ω = [0,1].^2 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  with semismooth newton algorithm      %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% using DG methods with Lgarange basis on quadrature nodes %%%%
%%%%%%%%%%%%%         u = u_final on ∂Ω              %%%%%%%%%%%%%%%%%%%%%%
% This main is used for tests of convergence for u, lap(u) and det(u) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global k
n     =  10;
h     = 1/n;
sigma =  100;
eps   =  -1;
k     =  3;


%f    =  @(x,y,t) x.*0;
f    =  @(x,y,t) source(x,y);
g = @(x,y,t) obstacle(x,y);
u_final = @(x,y,t) final_sol(x,y);

%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%
[A,Mass] = MatricesSystem(n,sigma,eps,k);
invMass = diag(diag(Mass).^-1);  %%=M^-1 where M is diago,al

b = SourceBCSystem(n,sigma,eps,k, g,f,0.0);
g_projected = diag(invMass).*b;

%return
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
    %use compute_Lap() to copmpute Laplacian
    %compute G[u^k] = -Lap(u^k) or a(u^k)=-Lap(u^k)-f 
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
    
    
    %plot_sol(n,k,0., U, g, u_final);
    
    
    res = h*sqrt(sum((U-Uold).^2));
    t = t + 1;
    fprintf('res=%f, counter=%i\n', res, t);
  
end
fprintf('\n')

compute_error(U,@(x,y) final_sol(x,y),n,k)

Lap = compute_Lap(A,Mass,U,b,f,0.,n,k);
%plot_sol(n,k,0., Lap, g, @(x,y,t) determiant_exct(x,y));

plot_det(n,k, U, Lap, @(x,y,t) determiant_exct(x,y));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = obstacle(tx,ty)

out = sin(pi*tx).*cos(2*pi*tx).^2  + sin(pi*ty).*cos(2*pi*ty).^2 ;

if (tx==0 || tx==1 || ty==0 || ty==1)
    out = final_sol(tx,ty);
end

end



function out = final_sol(tx,ty)

out = sin(pi*tx)  + sin(pi*ty);

end


function out = source(tx,ty)

out = pi^2*(sin(pi*tx) + sin(pi*ty));

end


function out = determiant_exct(tx,ty)
out =  pi^4*sin(pi*tx).*sin(pi*ty) ;

%out = -0.5*(-source(tx,ty) + pi^2*abs( -sin(pi*tx) + sin(pi*ty)  ) ) + source(tx,ty) ;
end