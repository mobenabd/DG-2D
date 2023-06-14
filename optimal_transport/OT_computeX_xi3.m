function [x,y] = OT_computeX_xi3(X0, Y0, xi_1,xi_2, LAMBDA, n,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute X(xi) with a Guass-Newton algo %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%x0 = [xi_1, xi_2];
x0 = [X0, Y0];

fun = @(X) func1(X, LAMBDA, xi_1,xi_2, n,k);

%options = optimoptions('fsolve','Display','none');
options=optimset('Display','none');% 'Algorithm', 'levenberg-marquardt');
%options.TolX = 1e-15;

[X, ~, EXITFLAG] = fsolve(fun,x0, options);

x = X(1);
y = X(2);

 if (EXITFLAG == - 2 )
     fprintf("WARNING: Not a root. x=%f, y=%f\n", x, y);
 end


end


function F = func1(X, LAMBDA, xi_1,xi_2,n,k)

% if (X(1)<0 || X(1)>1 || X(2)<0 || X(2)>1)
%     fprintf("WARNING: ux, uy out Omega \n");
%     %ux = 0; uy=0;
%    %ux = 0.1; uy=0.1;


test = 0;
if (X(1)<0 && X(2)>=0 && X(2)<=1 )
    X(1)
    [~, uy] = compute_grad(n,k,LAMBDA,0,X(2));
    ux = 0;
    test = 1;
elseif (X(1)>1 && X(2)>=0 && X(2)<=1 )
   	X(1)
    [~, uy] = compute_grad(n,k,LAMBDA,1,X(2));
    ux = 0;
    test = 1;
elseif (X(2)<0 && X(1)>=0 && X(1)<=1 )
    X(2)
    [ux, ~] = compute_grad(n,k,LAMBDA,X(1),0);
    uy = 0;
    test = 1;
elseif (X(2)>1 && X(1)>=0 && X(1)<=1 )
    X(2)
    [ux, ~] = compute_grad(n,k,LAMBDA,X(1),1);
    uy = 0;
    test = 1;
elseif ((X(1)<0 || X(1)>1 || X(2)<0 || X(2)>1) && test==0)
    ux=0;
    uy=0;
    test = 1;
end


if (test)
    fprintf("WARNING: ux, uy out Omega \n");
end

if (~test)
    [ux, uy] = compute_grad(n,k,LAMBDA,X(1),X(2));
end




F(1) = xi_1 + ux - X(1);
F(2) = xi_2 + uy - X(2);

end