function [dphxx, dphyy, dphxy] =  second_derv(index,n,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute global second deriviatives of basis functions on the refrence %%
%%% element [-1,1]*[-1,1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% using Lagrange shape functions on quadrature nodes %%%%%%%%%%%%%%%%%%%%
% index          : shape function index
% n              : Global dscretisation 
% k              : polynomial degree Q_k
% (dphxx, dphyy, dphxy) : = globla second deriviatives of phi(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mod(index, k+1) ==0 
    i = k;
else
    i =  mod(index, k+1) -1;
end
j = index/(k+1) - (i+1)/(k+1);


[Jac,~, ~] = jacobian_elem(1,n); %elem=1 (uniforme grid)
M = transpose(inv(Jac));

li_x = @(x) Lagrange(x,i+1,k,0);
lj_y = @(y) Lagrange(y,j+1,k,0);


dphxx = @(x,y) M(1,1)^2* Lagrange(x,i+1,k,2).* lj_y(y);
dphyy = @(x,y) M(2,2)^2* li_x(x).* Lagrange(y,j+1,k,2);
dphxy = @(x,y) M(1,1)*M(2,2)* Lagrange(x,i+1,k,1).* Lagrange(y,j+1,k,1);



end