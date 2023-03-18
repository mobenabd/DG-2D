function [A] = ACM_localMassMatrix_vol(n,k,num,c,U,LAP,f_projected,g_projected, DET, Eig_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% compute local mass matrix for ACM method %%%%%%%%%%%%%%%%%%%
% n  : Global discretisation 
% k  : polynomial degree
% num: elements number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloc = (k+1)^2;         %local dimension
A    = zeros(Nloc,Nloc);


Ak = ACM_ak(k,c,U,LAP, f_projected,g_projected,num, DET, Eig_m);
%ak  =  func_mapping(ak_function,n,num);
[~,~, detJac] = jacobian_elem(1,n);    % (uniforme grid)

for i=1:Nloc
    [phi, ~, ~] = basis(i,n,k);
    for j=1:Nloc
        [phj, ~, ~] = basis(j,n,k);
        f = @(x,y) phi(x,y).*phj(x,y);
        A(i,j) =  detJac * ACM_quadrature(Ak,f);
    end
end

end