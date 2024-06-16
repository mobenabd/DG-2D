clear; clc; close all;
addpath './src_DG/' './obstacle_problem'
%%%%%%%%%%%%% Concavization examples  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% using DG methods with Lgarange basis on quadrature nodes %%%%
%%%%%%%%%%%%%         u = g on ∂Ω              %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This main is used for tests without exact solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global k x0 xN hxy
n     =  8;         %% grid size in 1D direction
x0=-1; xN=1;        %% domain [x0,xN]^2
hxy=(xN-x0)/n;      %% grid size step
sigma =  10;        %% DG penalization parameter
eps   =  1;         %% DG symmetrization parameter
k    =  4;          %% polynomial degree
gdim = n^2*(k+1)^2; %% global dimention (degrees of freedom vector length)
c    = 1e4;         %% c parameter of the fixed point algorithm


f    = 0.1 ;                %% regulation term (source term)
g = @(x,y,t) obstacle(x,y); %% obstacle function
gxy = @(x,y) obstacle(x,y);
% L2 projection in the discrete space
g_projected = computeDirectProjection(n,k,g); 


U = solve_concav(gxy,g_projected, n, sigma, eps, c, f, true);
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = obstacle(tx,ty)
x = tx;
y = ty;
%%%%%%%%%% Concavization example 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu0 = 0;
out = -0.5*(x^2 + y^2) -0.5* exp(-50*((x-nu0)^2 + (y-nu0)^2));
%%%%%%%%%% Concavization example 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out = - (sin(pi*(x-y))./exp(2*(x.^2+y.^2)) + x.^2 + y.^2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return

%%%%%%%%%%% For Monge-Ampere solution %%%%%%%%%%%%%%
%out = abs(tx-0.5);
% x = 2*tx-1;
% y = 2*ty-1;
% out = exp((x.^2 + y.^2)/2);

%out = 2*sqrt(2)/3 * (tx.^2 + ty.^2).^(3/4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end








