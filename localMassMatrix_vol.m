function [A] = localMassMatrix_vol(n,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% compute local Mass matrix A  (volume integrals) %%%%%%%%%%%%%%%%%%%
% n  : Global discretisation 
% k  : polynomial degree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloc = (k+1)^2;         %local dimension
A    = zeros(Nloc,Nloc);

[~,~, detJac] = jacobian_elem(1,n);    % (uniforme grid)

for i=1:Nloc
    [phi, ~, ~] = basis(i,n,k);
    for j=1:Nloc
        [phj, ~, ~] = basis(j,n,k);
        A(i,j) =  detJac * quadrature(@(x,y) phi(x,y).*phj(x,y),2);
    end
end

end