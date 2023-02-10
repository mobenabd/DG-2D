function [x, y] = mapp_xy(xi, nu, elem, n)
%%%%%% Mapp (xi,yi) from refrence element to physical element %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Jac,b, ~] = jacobian_elem(elem,n);   


x = Jac(1,1)*xi + Jac(1,2)*nu + b(1);
y = Jac(2,1)*xi + Jac(2,2)*nu + b(2);


end