function [A] = localMatrix_vol(n,k, num, alpha_xy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% compute local matrix A  (volume integrals) %%%%%%%%%%%%%%%%%%%%%%%%
% n        : Global discretisation 
% k        : polynomial degree
% num     :
% alpha_xy : diffusion function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloc = (k+1)^2;         %local dimension
A    = zeros(Nloc,Nloc);

[~,~, detJac] = jacobian_elem(1,n);    % (uniforme grid)

if exist('alpha_xy', 'var')
    % mapping diffusion function to refrence element (xi,nu) -> foF(xi, nu)) = f(x,y)
    alpha_mapped  =  func_mapping(alpha_xy,n,num);  
else
    alpha_mapped = @(x,y) 1;
end



for i=1:Nloc
    [~, dphix, dphiy] = basis(i,n,k);
    for j=1:Nloc
        [~, dphjx, dphjy] = basis(j,n,k);
        A(i,j) =  detJac * quadrature(@(x,y) alpha_mapped(x,y).*(dphix(x,y).*dphjx(x,y) ...
                                           + dphiy(x,y).*dphjy(x,y)) ,2);
    end
end

end