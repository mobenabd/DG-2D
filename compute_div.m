function val = compute_div(n,k,U,x,y)
%%%% evaluate div(u) on point (x,y) %%%%%%%%%%%%%%%%
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


ux = 0; %%=\partial_{x}u_1
uy = 0; %%=\partial_{y}u_2

for i=1:Nloc
    uix =  U((E-1)*Nloc+i,1);
    uiy =  U((E-1)*Nloc+i,2);
    
    [~, dphx, dphy] = basis(i,n,k);
    
    ux = ux + uix*dphx(xi,nu);
    uy = uy + uiy*dphy(xi,nu);
end


val = ux + uy;

end