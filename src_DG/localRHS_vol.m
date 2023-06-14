function b = localRHS_vol(f,n,k,num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% compute local rhs b (volume integrals) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f   : source term function 
% n   : Global discretisation 
% k   : polynomial degree Q_k
% num : element number 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloc = (k+1)^2;         %local dimension
b    = zeros(Nloc,1);

% mapping source function to refrence element (xi,nu) -> foF(xi, nu)) = f(x,y)
source  =  func_mapping(f,n,num);
[~,~, detJac] = jacobian_elem(1,n);    % (uniforme grid)
for i=1:Nloc
    [phi,~, ~] = basis(i,n,k);
    b(i) = detJac * quadrature(@(x,y) source(x,y).*phi(x,y),2,3);
end
end