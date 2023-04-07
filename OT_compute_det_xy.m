function val = OT_compute_det_xy(x,y,n,k,U, eig_m)
%%%% evaluate det(I +/- \nabla^2 U) on point (x,y) %%%%%%
%%%%%% if eig==1, compute det(I-\nabla^2 U)(x)
%%%%%% if eig==2, compute det(I+\nabla^2 U)(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    val = (1-uxx)*(1-uyy) - uxy^2;
else
    val = (1+uxx)*(1+uyy) - uxy^2;
end





end