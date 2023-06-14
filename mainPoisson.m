clear; clc; close all;
%%%%%%%%%%%%%%%  Solve u_t-Δu = f in Ω = [0,1].^2 %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%         u = u_0  in Ω at t=0      %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%         u = g on ∂Ω               %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n     =  16;
sigma =  100;
eps   =  -1;
global k
k     =  2;
global x0 xN hxy
x0 = -2; xN = 2; hxy = (xN-x0)/n;


%exact solution  + BC condition + initial condition
uexct =  @(x,y,t)  x.*y.*(1-x).*(1-y);
%uexct =  @(x,y,t) sin(2*pi*x).*sin(2*pi*y);

%source term function
f    =  @(x,y,t) 2*(x-x.*x + y-y.*y);
%f    =  @(x,y,t) 8*pi^2*(sin(2*pi*x).*sin(2*pi*y));


%%% Solve linear system %%%%%%%%%%%%%%%%%%%%%
[A,Mass] = MatricesSystem(n,sigma,eps,k);
b = SourceBCSystem(n,sigma,eps,k,uexct,f,1);
U = A\b;
plot_solHeat(n,k,0, U, uexct);



%%% compute weak laplacian
Lap = compute_Lap(A,Mass,U,b,f,n,k);

return
compute_error(U,@(x,y) uexct(x,y,0),n,k)
compute_error(Lap,@(x,y) -f(x,y,0) ,n,k)

