function b = localRHS_edgeBoundary(verHor,ort,num,n,k,eps,sigma,uexct, alpha_xy)
%%% routine: Non homegenous Dirichlet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n        : Global discretisation 
% k        : polynomial degree Q_k
% eps      : symmetrization parameter
% sigma    : penelisation parameter
% verHor   : =1 if vertical interior edge. =0 if horizontale interior edge.
% ort      : =1 if left/bottom boundary, 0 if right/top boundary
% alpha_xy : diffusion function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('alpha_xy', 'var')
     % mapping diffusion function to refrence element (xi,nu) -> foF(xi, nu)) = f(x,y)
    alpha_mapped = func_mapping(alpha_xy,n,num);
else
    alpha_mapped = @(x,y) 1;
end

[~, ~, h] = getGlobal_x0N();
Nloc = (k+1)^2;
beta0 = 1;       %supper-penalization =1 (for now)
b = zeros(Nloc,1);

if (ort==1)
    normal = -1;
elseif (ort==0)
    normal = 1;
end
ss = normal;


%[Jac,~, ~] = jacobian_elem(1,n) ; %(uniforme grid)
%detJac11 = Jac(1,1);
detJac11 = h/2;

uexct_loc  =  func_mapping(uexct,n,num);

for i=1:Nloc
    [phi, dphix, dphiy] = basis(i,n,k);
    if (verHor ==1 )   %vertical edge on x=+/-1, normal = [+/-1 0];
        b(i) = detJac11* quadrature(@(y) (normal*alpha_mapped(ss,y)*eps*dphix(ss,y)  ...
                                    + 2*sigma/(h^beta0)*phi(ss,y)).*uexct_loc(ss,y),1);
        
    elseif (verHor == 0)  %horizontal edge on y=+/-1, normal = [0 +/-1];
        b(i) = detJac11* quadrature(@(x) (normal*alpha_mapped(x,ss)*eps*dphiy(x,ss) ....
                                    + 2*sigma/(h^beta0)*phi(x,ss)).*uexct_loc(x,ss),1);
    end
end


end