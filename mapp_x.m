function x = mapp_x(xi, elem, n)
%%%%%% Mapp (xi) from refrence element to physical element %%%%%%%%%
%%%%%% 1D mapping 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Jac,b, ~] = jacobian_elem(elem,n);   


x = Jac(1,1)*xi  + b(1);


end