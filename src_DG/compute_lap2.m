function [val_lap] = compute_lap2(n,k,U,x,y)
%%%% evaluate u and lap(u) on point (x,y) %%%%%%%%%%%%%%%%
% using basis deriviatives %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloc = (k+1)^2;
E = get_elemPos(x,y,n);

[xi,nu] = inverse_mapp(x,y,E,n);

%val_sol = 0;
uxx = 0; %%=\partial_{xx}u
uyy = 0; %%=\partial_{yy}u

for i=1:Nloc
    ui =  U((E-1)*Nloc+i);
    [dphxx, dphyy, ~] =  second_derv(i,n,k);
    uxx = uxx + ui*dphxx(xi,nu);
    uyy = uyy + ui*dphyy(xi,nu);
    %[phi, ~, ~] = basis(i,n,k);
    %val_sol = val_sol +ui*phi(xi,nu);
end


val_lap = uxx + uyy;

end