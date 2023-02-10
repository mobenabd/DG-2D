clear; clc; close all;
%%%%%%%%%%%%%%%  Solve u_t-Δu = f in Ω = [0,1].^2 %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%         u = u_0  in Ω at t=0      %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%         u = g on ∂Ω               %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n     =  20;
h     = 1/n;
sigma =  100;
eps   =  1;
k     =  1;
dt    = h^2/4;

f    =  @(x,y,t) x.*0;
g = @(x,y,t) obstacle(x,y);
u_final = @(x,y,t) final_sol(x,y);

%%% Solve with SemiSmooth Newton %%%%%%%%%%%%%%%%%%%%%
[A,Mass] = MatricesSystem(n,sigma,eps,k);
%initial contition
b = SourceBCSystem(n,sigma,eps,k,g,f,0.0);
U = Mass\b;
plot_sol(n,k,0, U, g, u_final);

return


gdim = length(U);
G = zeros(gdim,1);
dG = zeros(gdim,gdim);

Uold = zeros(gdim,1);
g_projected = zeros(gdim,1);
g_projected = U;
res = 1;
t=1;
%U = ones(gdim,1);
while (res > 1e-8)
    b = SourceBCSystem(n,sigma,eps,k,g,f,t*dt);
    G = A*U - b;
    %compute H && Da && Db: %a=-lap(u)-f   %b=u-g
    [G,Da,Db] = updateNewton_obstacle(n,k, G, U, g);
    
    
    %  Compute dH(u^k)
    dG = Da*A + Db;
    
    Uold = U;
    U = -dG\G;
    U = Uold + U;
    
    
    
    %plot_sol(U, u_final,n,k,t*dt,1);
    plot_sol(n,k,t*dt, U, g, u_final);
    res = h*sqrt(sum((U-Uold).^2));
    
    
    t = t + 1;
    fprintf('res=%f, counter=%i\n', res, t);
  
end
fprintf('\n')



function out = obstacle(tx,ty)
x = 4*tx-2;
y = 4*ty-2;


rxy = x.^2 + y.^2;
if (rxy <= 1 )
    out = sqrt(1 - rxy);
else
    out = -1;
end

%return
if ((x==-2 ||  x==2 || y==-2 || y==2))
    out =  final_sol(tx,ty);
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
