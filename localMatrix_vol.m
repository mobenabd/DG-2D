function [A] = localMatrix_vol(n,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% compute local matrix A  (volume integrals) %%%%%%%%%%%%%%%%%%%%%%%%
% n  : Global discretisation 
% k  : polynomial degree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloc = (k+1)^2;         %local dimension
A    = zeros(Nloc,Nloc);

[~,~, detJac] = jacobian_elem(1,n);    % (uniforme grid)

for i=1:Nloc
    [~, dphix, dphiy] = basis(i,n,k);
    for j=1:Nloc
        [~, dphjx, dphjy] = basis(j,n,k);
        A(i,j) =  detJac * quadrature(@(x,y) dphix(x,y).*dphjx(x,y) ...
                                           + dphiy(x,y).*dphjy(x,y) ,2);
    end
end

end