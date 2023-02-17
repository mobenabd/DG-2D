function b = localRHS_edgeBoundary(verHor,ort,num,n,k,eps,sigma,uexct)
%%% routine: Non homegenous Dirichlet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n       : Global discretisation 
% k       : polynomial degree Q_k
% eps     : symmetrization parameter
% sigma   : penelisation parameter
% verHor  : =1 if vertical interior edge. =0 if horizontale interior edge.
% ort     : =1 if left/bottom boundary, 0 if right/top boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


h = 1/n;
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
        b(i) = detJac11* quadrature(@(y) (eps*normal*dphix(ss,y)  + sigma/(h^beta0)*phi(ss,y)).*uexct_loc(ss,y),1,5);
        
    elseif (verHor == 0)  %horizontal edge on y=+/-1, normal = [0 +/-1];
        b(i) = detJac11* quadrature(@(x) (eps*normal*dphiy(x,ss) + sigma/(h^beta0)*phi(x,ss)).*uexct_loc(x,ss),1,5);
    end
end


end