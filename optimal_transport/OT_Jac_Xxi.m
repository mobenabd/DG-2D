function dx = OT_Jac_Xxi(xi, ss)
%%%%% compute the diagonal component of the jacobian %%%%%%%%%%%%%%%
%%%%% \partial X(xi) for X(xi)  known (the symetric test case) %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù

delta = (1+ss)^2 - 4*ss*xi;

dx = 1/sqrt(delta);

end