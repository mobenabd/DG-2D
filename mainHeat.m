clear; clc; close all;
addpath './src_DG/'
%%%%%%%%%%%%%%%  Solve u_t-Δu = f in Ω = [x0,xN].^2 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%         u = u_0  in Ω at t=0      %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%         u = g on ∂Ω               %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global k x0 xN hxy
n     =  16;        %% grid size in 1D direction
x0 = 0; xN = 1;     %% domain [x0,xN]^2
hxy = (xN-x0)/n;    %% grid size step
sigma =  100;       %% DG penalization parameter
eps   =  -1;        %% DG symmetrization parameter
k     =  1;         %% polynomial degree 
dt = 0.1;           %% time step


%exact solution  + BC condition + initial condition
%uexct =  @(x,y,t)  x.*y.*(1-x).*(1-y).*cos(pi/2. *t);
uexct =  @(x,y,t) sin(2*pi*x).*sin(2*pi*y).*cos(pi/2. *t);

%source term function
%f    =  @(x,y,t) 2*(x-x.*x + y-y.*y)*cos(pi/2. *t) - pi/2.*sin(pi/2. *t).*x.*y.*(1-x).*(1-y);
f    =  @(x,y,t) 8*pi^2*(sin(2*pi*x).*sin(2*pi*y))*cos(pi/2. *t) - pi/2.*sin(pi/2. *t).*sin(2*pi*x).*sin(2*pi*y);


%%% Solve with Backward Euler %%%%%%%%%%%%%%%%%%%%%
[A,Mass] = MatricesSystem(n,sigma,eps,k);
%initial contition
b = SourceBCSystem(n,sigma,eps,k,uexct,f,0);
U = Mass\b;
plot_solHeat(n,k,0, U, uexct);



count=1;
while (count<10)
    b = SourceBCSystem(n,sigma,eps,k,uexct,f,count*dt);
    U = (Mass + dt*A)\(Mass*U+dt*b);
    plot_solHeat(n,k,count*dt, U, uexct);
    
    
    %%% compute and plot weak laplacian of u %%%%%
    %Lap = compute_Lap(A,Mass,U,b,f,count*dt,n,k);
    %plot_solHeat(n,k,count*dt, Lap, @(x,y,t) -f(x,y,t) - pi/2.*sin(pi/2. *t).*sin(2*pi*x).*sin(2*pi*y));
    
    
    count = count + 1;
    fprintf('Count = %i \n', count);
end
t=(count-1)*dt;


compute_error(U,@(x,y) uexct(x,y,t),n,k)


