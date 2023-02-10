function [phi, dphix, dphiy] = basis(index,n,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute local basis functions on the refrence element [-1,1]*[-1,1] %%%
% index          : shape function index
% n              : Global dscretisation (h=1/n)
% phi            : fonction de base sur l'element de reference
% (dphix, dphiy) : =Jac^(-T)*nabla.(phi) la derivée globale de phi
% det_Jac        : determinant de la jacobienne
% k              : polynomial degree Q_k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% (with local deriviatives)
if (k==1)
    if (index==1)
        phi       = @(x,y) 1/4*(1-x).*(1-y);
        dphix_loc = @(x,y) -1/4*(1-y);
        dphiy_loc = @(x,y) -1/4*(1-x);
    elseif (index==2)
        phi       = @(x,y) 1/4*(1+x).*(1-y);
        dphix_loc = @(x,y) 1/4*(1-y);
        dphiy_loc = @(x,y) -1/4*(1+x);
    elseif (index==3)
        phi       = @(x,y) 1/4*(1+x).*(1+y);
        dphix_loc = @(x,y) 1/4*(1+y);
        dphiy_loc = @(x,y) 1/4*(1+x);
    elseif (index==4)
        phi   = @(x,y) 1/4*(1-x).*(1+y);
        dphix_loc = @(x,y) -1/4*(1+y);
        dphiy_loc = @(x,y) 1/4*(1-x);
    end
elseif (k==2)
    Xi = [-1 1 1 -1 0 0 -1 1 0]; Nu = [-1 -1 1 1 -1 1 0 0 0];
    xi = Xi(index);            nu = Nu(index);
    if (index<=4)
        phi       = @(x,y) 1/4 * x*xi.*(1+x*xi).*y*nu.*(1+nu*y);
        dphix_loc = @(x,y) 1/4*nu*xi*(1+2*xi*x).*(y+nu*y.^2);
        dphiy_loc = @(x,y) 1/4*nu*xi*(1+2*nu*y).*(x+xi*x.^2);
    elseif (index==5 || index==6)
        phi       = @(x,y) 1/2 * y*nu.*(1+y*nu).*(1-x.^2);
        dphix_loc = @(x,y) -y*nu.*(1+y*nu).*x;
        dphiy_loc = @(x,y) 1/2 * nu*(1-x.^2).*(1+2*nu*y);
    elseif (index==7 || index==8)
        phi       = @(x,y) 1/2 * x*xi.*(1+x*xi).*(1-y.^2);
        dphix_loc = @(x,y) 1/2 * xi*(1+2*xi*x).*(1-y.^2);
        dphiy_loc = @(x,y) -y*xi.*x.*(1+x*xi);
    elseif (index==9)
        phi       = @(x,y) (1-x.^2).*(1-y.^2);
        dphix_loc = @(x,y) -2*x.*(1-y.^2);
        dphiy_loc = @(x,y) -2*y.*(1-x.^2);
    end
elseif (k==3)
    Xi = [-1 1 1 -1 -1/3 1/3 1/3 -1/3 1 1 -1 -1 -1/3 1/3 1/3 -1/3]; 
    Nu = [-1 -1 1 1 -1 -1 1 1 -1/3 1/3 1/3 -1/3 -1/3 -1/3 1/3 1/3];
    xi = Xi(index);            nu = Nu(index);
    if (index <=4)
        phi       = @(x,y) 81/256 * (1+xi*x).*(1+nu*y).*(1/9-x.^2).*(1/9-y.^2);
        dphix_loc = @(x,y) 81/256 *(1+nu*y).*(1/9-y.^2).*(-2*x + 1/9*xi -3*xi*x.^2);
        dphiy_loc = @(x,y) 81/256 *(1+xi*x).*(1/9-x.^2).*(-2*y + 1/9*nu -3*nu*y.^2);
    elseif (index>=5 && index<=8)
        phi       = @(x,y) 243/256 * (1-x.^2).*(y.^2 - 1/9).*(1/3 + 3*xi*x).*(1+nu*y);
        dphix_loc = @(x,y) 243/256 * (3*xi -2/3*x -9*xi*x.^2).*(y.^2 - 1/9).*(1+nu*y);
        dphiy_loc = @(x,y) 243/256 * (2*y + 3*nu*y.^2 - 1/9*nu).*(1-x.^2).*(1/3 + 3*xi*x);
    elseif (index>=9 && index<=12)
        phi       = @(x,y) 243/256 * (1-y.^2).*(x.^2 - 1/9).*(1/3 + 3*nu*y).*(1+xi*x);
        dphix_loc = @(x,y) 243/256 * (2*x + 3*xi*x.^2 - 1/9*xi).*(1-y.^2).*(1/3 + 3*nu*y);
        dphiy_loc = @(x,y) 243/256 * (3*nu -2/3*y -9*nu*y.^2).*(x.^2 - 1/9).*(1+xi*x);
    elseif (index>=13 && index<=16)
        phi       = @(x,y) 729/256 * (1-x.^2).*(1-y.^2).*(1/3 + 3*xi*x).*(1/3 + 3*nu*y);
        dphix_loc = @(x,y) 729/256 * (3*xi -2/3*x -9*xi*x.^2).*(1-y.^2).*(1/3 + 3*nu*y);
        dphiy_loc = @(x,y) 729/256 * (3*nu -2/3*y -9*nu*y.^2).*(1-x.^2).*(1/3 + 3*xi*x);
    end
end


[Jac,~, ~] = jacobian_elem(1,n); %elem=1 (uniforme grid)
M = transpose(inv(Jac));

dphix  = @(x,y) M(1,1)*dphix_loc(x,y) + M(1,2)*dphiy_loc(x,y);
dphiy = @(x,y) M(2,1)*dphix_loc(x,y)  + M(2,2)*dphiy_loc(x,y);


end