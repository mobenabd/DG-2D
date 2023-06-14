function val = compute_sol(x,y,sol,n,k)
%%%% evaluate solution on point (x,y) %%%%%%%%%%

Nloc = (k+1)^2;

E = get_elemPos(x,y,n);

[xi,nu] = inverse_mapp(x,y,E,n);

val = 0;
 for i=1:Nloc
    [phi, ~, ~] = basis(i,n,k);
     val = val + sol((E-1)*Nloc+i)*phi(xi,nu);
 end

end