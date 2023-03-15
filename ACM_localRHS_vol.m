function b = ACM_localRHS_vol(n,k,num,c,U,LAP,f_projected,g_projected)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% compute local mass matrix for ACM method %%%%%%%%%%%%%%%%%%%
% n  : Global discretisation 
% k  : polynomial degree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloc = (k+1)^2;         %local dimension
b    = zeros(Nloc,1);


Bk = ACM_bk(k,c,U,LAP, f_projected,g_projected,num);
%ak  =  func_mapping(ak_function,n,num);
[~,~, detJac] = jacobian_elem(1,n);    % (uniforme grid)

for i=1:Nloc
    [phi, ~, ~] = basis(i,n,k);
    b(i) =  detJac * ACM_quadrature(Bk, phi);
end

end