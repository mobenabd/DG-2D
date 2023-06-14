function func = func_mapping(f,n,num)   %func(xi,nu) = f(F(xi,nu))
%%% subroutine : define compose function (xi, nu) -> foF(xi,nu) %%%%%%%%%%%
%%% F: mapping function from refrence to physical element
% (general case)
% n            : global descritisation
% num           : element number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[Jac,b, ~] = jacobian_elem(num,n);    

%func = @(xi,nu) f( Jac(1,1)*xi + Jac(1,2)*nu + b(1), Jac(2,1)*xi + Jac(2,2)*nu + b(2) );

func = @(xi,nu) f( Jac(1,1)*xi + b(1), Jac(2,2)*nu + b(2) );

end

