function [phi, dphix, dphiy] = basis(index,n,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute local basis functions on the refrence element [-1,1]*[-1,1] %%%
%%% using Lagrange shape functions on quadrature nodes %%%%%%%%%%%%%%%%%%%%
% index          : shape function index
% n              : Global dscretisation (h=(xN-x0)/n)
% phi            : fonction de base sur l'element de reference
% (dphix, dphiy) : =Jac^(-T)*nabla.(phi) la derivée globale de phi
% det_Jac        : determinant de la jacobienne
% k              : polynomial degree Q_k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mod(index, k+1) ==0 
    i = k;
else
    i =  mod(index, k+1) -1;
end
j = index/(k+1) - (i+1)/(k+1);


li_x = @(x) Lagrange(x,i+1,k,0);
lj_y = @(y) Lagrange(y,j+1,k,0);

phi = @(x,y) li_x(x).* lj_y(y);
dphix_loc = @(x,y) Lagrange(x,i+1,k,1).* lj_y(y);
dphiy_loc = @(x,y) li_x(x).* Lagrange(y,j+1,k,1);




[Jac,~, ~] = jacobian_elem(1,n); %elem=1 (uniforme grid -> same jacobian)
M = transpose(inv(Jac));

%dphix  = @(x,y) M(1,1)*dphix_loc(x,y) + M(1,2)*dphiy_loc(x,y);
%dphiy = @(x,y) M(2,1)*dphix_loc(x,y)  + M(2,2)*dphiy_loc(x,y);


dphix  = @(x,y) M(1,1)*dphix_loc(x,y);
dphiy = @(x,y)  M(2,2)*dphiy_loc(x,y);

end