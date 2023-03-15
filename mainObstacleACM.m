clear; clc; close all;
%%%%%%%%%%%%% Solve min(-Δu, u-g)=0 in Ω = [0,1].^2 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  with ACM method      %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% using DG methods with Lgarange basis on quadrature nodes %%%%
%%%%%%%%%%%%%         u = u_final on ∂Ω              %%%%%%%%%%%%%%%%%%%%%%
% This main is used for tests without exact solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global k
n     =  8;
h     = 1/n;
sigma =  1;
eps   =  1;
k    =  4;
c    = 1000;
gdim = n^2*(k+1)^2;


f    =  @(x,y,t) x.*0  ;
%f    =  @(x,y,t) source(x,y);
g = @(x,y,t) obstacle(x,y);
u_final = @(x,y,t) final_sol(x,y);

%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%
[A,Mass] = MatricesSystem(n,sigma,eps,k);

b = SourceBCSystem(n,sigma,eps,k, g,f,0.0);
g_projected = Mass\b;

Uold = zeros(gdim,1);
U = g_projected;
plot_sol(n,k, U, g);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lap = computeLap2_ddl(n,k,U);
Eig_m = computeDet_ddl(n,k,U, 1);  %%compute largest eigenvalue
f_projected = -Lap + Eig_m ;

f_projected = zeros(gdim,1);
 
% plot_sol(n,k, f_projected);
% return

res = 1;
t=1;
while (res > 1e-4)
    Lap = computeLap2_ddl(n,k,U);
    Mk =  ACM_GlobalMass(n,k,c,U,Lap,g_projected, f_projected);
    RHS = ACM_RHS(n,k,c,U,Lap,sigma,eps,g,g_projected, f_projected);
    
    Uold = U;
    U = (A+Mk)\RHS;
    plot_sol(n,k, U, g);
    %return
    
    
    res = sqrt(sum((U-Uold).^2))/sqrt(sum(U.^2));
    t = t + 1;
    fprintf('res=%f, counter=%i\n', res, t);
  
end
fprintf('\n')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = obstacle(tx,ty)


out = sin(2*pi*tx)*sin(2*pi*ty);

%out = (tx.^2 - ty.^2).*sin(pi*tx).*sin(2*pi*ty);

% if (out<0)
%      out = 0;
%  end

% x = 2*tx-1;
% y = 2*ty-1;
% 
% out = -(0.98 - x.^2 - y.^2)^2;


end




function out = final_sol(tx,ty)

out = 0.;

end




