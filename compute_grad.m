function [ux, uy] = compute_grad(n,k,U,x,y)
%%%% evaluate grad(u) on point (x,y) %%%%%%%%%%%%%%%%
% using basis deriviatives %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloc = (k+1)^2;
h = 1/n;
j = floor(x/h)+1; i = floor(y/h)+1;
if (x==1)
    j = j-1;
end
if (y==1)
    i = i-1;
end
E = (i-1)*n+j;

[xi,nu] = inverse_mapp(x,y,E,n);


ux = 0; %%=\partial_{x}u
uy = 0; %%=\partial_{y}u

for i=1:Nloc
    ui =  U((E-1)*Nloc+i);
    [~, dphx, dphy] = basis(i,n,k);
    ux = ux + ui*dphx(xi,nu);
    uy = uy + ui*dphy(xi,nu);
end


end