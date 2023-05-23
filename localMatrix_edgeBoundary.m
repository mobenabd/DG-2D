function [M11] = localMatrix_edgeBoundary(verHor,ort,n,k,eps,sigma, num, alpha_xy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% compute local matrix M_11 (integrals over boundary edge ) %%%%%%%%%%%%%%%%%%
% n        : Global discretisation 
% k        : polynomial degree Q_k
% eps      : symmetrization parameter
% sigma    : penelisation parameter
% verHor   : =1 if vertical interior edge. =0 if horizontale interior edge.
% ort      : =1 if left/bottom boundary, 0 if right/top boundary
% num      : 
% alpha_xy : diffusion function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('alpha_xy', 'var')
     % mapping diffusion function to refrence element (xi,nu) -> foF(xi, nu)) = f(x,y)
    alpha_mapped = func_mapping(alpha_xy,n,num);
else
    alpha_mapped = @(x,y) 1;
end



h = 1/n;
Nloc = (k+1)^2;
beta0 = 1;     %supper-penalization =1 (for now)
M11 = zeros(Nloc,Nloc);

if (ort==1)
    normal = -1;
elseif (ort==0)
    normal = 1;
end
ss = normal;


%[Jac,~, ~] = jacobian_elem(1,n) ; %(uniforme grid)
%detJac11 = Jac(1,1);
detJac11 = h/2;

for i=1:Nloc
    for j=1:Nloc
        [phi, dphix, dphiy] = basis(i,n,k);
        [phj, dphjx, dphjy] = basis(j,n,k);
        
        if (verHor ==1 )   %vertical edge on x=+/-1, normal = [+/-1 0];
           M11(i,j) = detJac11* quadrature(@(y) alpha_mapped(ss,y)*normal*(-dphjx(ss,y) .* phi(ss,y) ...
                       + eps*dphix(ss,y) .* phj(ss,y)) ...
                       + (2*sigma/h^beta0)*phi(ss,y).*phj(ss,y),1);
        elseif (verHor == 0)  %horizontal edge on y=+/-1, normal = [0 +/-1];
           M11(i,j) = detJac11*quadrature(@(x) alpha_mapped(x,ss)*normal*(-dphjy(x,ss) .* phi(x,ss) ...
                       + eps*dphiy(x,ss) .* phj(x,ss)) ...
                       + (2*sigma/h^beta0)*phi(x,ss).*phj(x,ss),1);
        end
    end
end

end