function b = localRHS_t0_vol(u0,n,k,num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% compute local rhs b (volume integrals) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% u0  : initial condition
% n   : Global discretisation 
% k   : polynomial degree Q_k
% num : element number 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloc = (k+1)^2;         %local dimension
b    = zeros(Nloc,1);

% mapping u0=f function to refrence element (xi,nu) -> foF(xi, nu)) = f(x,y)
u0Mapped  =  func_mapping(u0,n,num);
[~,~, detJac] = jacobian_elem(1,n);    % (uniforme grid)
for i=1:Nloc
    [phi,~, ~] = basis(i,n,k);
    b(i) = detJac * quadrature(@(x,y) u0Mapped(x,y).*phi(x,y),2);
end
end