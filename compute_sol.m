function val = compute_sol(x,y,sol,n,k)
%%%% evaluate solution on point (x,y) %%%%%%%%%%

Nloc = (k+1)^2;
h = 1/n;

j = floor(x/h)+1;
i = floor(y/h)+1;
if (x==1)
    j = j-1;
end
if (y==1)
    i = i-1;
end


E = (i-1)*n+j;

[xi,nu] = inverse_mapp(x,y,E,n);

val = 0;
 for i=1:Nloc
    [phi, ~, ~] = basis(i,n,k);
     val = val + sol((E-1)*Nloc+i)*phi(xi,nu);
 end

end