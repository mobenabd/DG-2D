clear; clc; close all;
%%%%%%%%%%%%%%%  Solve u_t-Δu = f in Ω = [0,1].^2 %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%         u = u_0  in Ω at t=0      %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%         u = g on ∂Ω               %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n     =  20;
sigma =  100;
eps   =  -1;
k     =  1;
dt = 0.01;

%exact solution  + BC condition + initial condition
uexct =  @(x,y,t)  x.*y.*(1-x).*(1-y).*cos(pi/2. *t);
%source term function
f    =  @(x,y,t) 2*(x-x.*x + y-y.*y)*cos(pi/2. *t) - pi/2.*sin(pi/2. *t).*x.*y.*(1-x).*(1-y);


%%% Solve with Backward Euler %%%%%%%%%%%%%%%%%%%%%
[A,Mass] = MatricesSystem(n,sigma,eps,k);
%initial contition
b = SourceBCSystem(n,sigma,eps,k,uexct,f,0);
U = Mass\b;
plot_sol(U, uexct,n,k,0);

t=1;
while (t<200)
    b = SourceBCSystem(n,sigma,eps,k,uexct,f,t*dt);
    U = (Mass + dt*A)\(Mass*U+dt*b);
    
    plot_sol(U, uexct,n,k,t*dt);
    t = t + 1;
    %pause(0.1);
end

