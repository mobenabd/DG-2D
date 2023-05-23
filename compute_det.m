function val = compute_det(x,y,n,k,U, eig_m)
%%%% evaluate det(u) on point (x,y) %%%%%%%%%%%%%%%%
% Solution U ddls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lap: weak laplacian ddls %%%%%%%%%%%%%%%%%%%%%%%%%
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

uxx = 0; %%=\partial_{xx}u
uyy = 0; %%=\partial_{yy}u
uxy = 0; %%=\partial_{xy}u
for i=1:Nloc
    [dphxx, dphyy, dphxy] =  second_derv(i,n,k);
    
    uxx = uxx + U((E-1)*Nloc+i)*dphxx(xi,nu);
    uyy = uyy + U((E-1)*Nloc+i)*dphyy(xi,nu);
    uxy = uxy + U((E-1)*Nloc+i)*dphxy(xi,nu);
end



if (eig_m==1)
    %val = 0.5*(compute_sol(x,y,Lap,n,k) + sqrt( (uxx-uyy)^2 + 4*uxy^2 )); %largest eigenvalue
    val = 0.5*(uxx + uyy + sqrt( (uxx-uyy)^2 + 4*uxy^2 )); %largest eigenvalue
    return;
elseif (eig_m==2)
    val = 0.5*(uxx + uyy - sqrt( (uxx-uyy)^2 + 4*uxy^2 )); %smallest eigenvalue
    return;
end

%val1 = 0.5*(uxx + uyy + sqrt( (uxx-uyy)^2 + 4*uxy^2 )); %largest eigenvalue
%val2 = 0.5*(uxx + uyy - sqrt( (uxx-uyy)^2 + 4*uxy^2 )); %smallest eigenvalue

%val = val1*val2;

val = uxx*uyy - uxy^2; %%for testing without weak laplacian


end