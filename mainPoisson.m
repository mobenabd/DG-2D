clear; clc; close all;
addpath './src_DG/'
%%%%%%%%%%%%%%%  Solve u_t-Δu = f in Ω = [x0,xN].^2 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%         u = u_0  in Ω at t=0      %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%         u = g on ∂Ω               %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global k x0 xN hxy
n     =  16;     %% grid size in 1D direction
x0 = 0; xN = 1; %% domain [x0,xN]^2
hxy = (xN-x0)/n; %% grid size step
sigma =  100;    %% DG penalization parameter
eps   =  -1;     %% DG symmetrization parameter
k     =  2;      %% polynomial degree 


%exact solution  + BC condition + initial condition
uexct =  @(x,y,t)  x.*y.*(1-x).*(1-y);
%uexct =  @(x,y,t) sin(2*pi*x).*sin(2*pi*y);

%source term function
f    =  @(x,y,t) 2*(x-x.*x + y-y.*y);
%f    =  @(x,y,t) 8*pi^2*(sin(2*pi*x).*sin(2*pi*y));


%%% Solve linear system and plot %%%%%%%%%%%%%%%%%%%%%
[A,Mass] = MatricesSystem(n,sigma,eps,k);
b = SourceBCSystem(n,sigma,eps,k,uexct,f,1);
U = A\b;
plot_solHeat(n,k,0, U, uexct);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% compute weak laplacian
Lap = compute_Lap(A,Mass,U,b,f,n,k);


compute_error(U,@(x,y) uexct(x,y,0),n,k)
compute_error(Lap,@(x,y) -f(x,y,0) ,n,k)

