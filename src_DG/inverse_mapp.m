function [xi,nu] = inverse_mapp(x,y,E,n)


[Jac,bE, ~] = jacobian_elem(E,n);

M = inv(Jac);
bE = M*bE;

xi = M(1,1)*x -bE(1);
nu = M(2,2)*y -bE(2);

end