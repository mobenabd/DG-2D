function y = mapp_y(nu, elem, n)
%%%%%% Mapp yi from refrence element to physical element %%%%%%%%%
%%%%%% 1D mapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Jac,b, ~] = jacobian_elem(elem,n);   

y =  Jac(2,2)*nu + b(2);


end