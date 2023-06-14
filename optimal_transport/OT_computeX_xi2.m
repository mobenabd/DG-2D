function [x,y] = OT_computeX_xi2(xi_1,xi_2, LAMBDA, n,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute X(xi) with a fixe point algorithme %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iter=0;
x = xi_1; y=xi_2;
res = 1;
% while(res>1e-3 && iter<50)
%     xold = x; yold=y;
%     if (x>=0 && x<=1 && y>=0 && y<=1)
%         xold = x; yold=y;
%         [ux, uy] = compute_grad(n,k,LAMBDA,x,y); %%\nabla_x \lambda(x)
%         x = xi_1 + ux;
%         y = xi_2 + uy;
%         
%         res = sqrt((x-xold)^2 + (y-yold)^2)/sqrt(x^2 + y^2);
%     else
%         x = xold;
%         y = yold;
%         return;
%     end
%     iter = iter+1;
% end

s= 0.7;

while(res>1e-4 && iter<100)
    if (x>=0 && x<=1 && y>=0 && y<=1)
        xold = x; yold=y;
        [ux, uy] = compute_grad(n,k,LAMBDA,x,y);
        [uxx, uyy, ~] = compute_secondDerv(n,k,LAMBDA,x,y); %%changed without validation
        x = x - (xi_1 + ux -x)/(uxx-1);
        y = y - (xi_2 + uy -y)/(uyy-1);
        
        x = s*xold + (1-s)*x;
        y = s*yold + (1-s)*y;
        
        res = sqrt((x-xold)^2 + (y-yold)^2)/sqrt(x^2 + y^2);
    else
        fprintf("WARNING: Newton diverged, x=%f, y=%f, xold=%f, yold=%f\n", x,y, xold, yold);
        x = xold;
        y = yold;
        return;
    end
    %fprintf("xi1=%f, xi2=%f, x=%f, y=%f, res=%f, iter=%f\n",xi_1, xi_2, x, y, res, iter);
    iter = iter+1;
end


end