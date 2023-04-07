function PROJ = OT_compute_Dlambda(n,PHI, LAMBDA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% compute Dlambda(x) = phi(x - \nabla_x\lambda(x))
%%%%%% on quadrature nodes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k    =  getGlobal_k();
Nloc = (k+1)^2;
gdim = n^2*Nloc;
PROJ = zeros(gdim,1);

[nodes, ~] = getWeightsNodes(k+1);


kk = 0;
for elem=1:n^2
    for i=1:Nloc
        ie = i + kk;
        
        if mod(i, k+1) ==0
            ii = k;
        else
            ii =  mod(i, k+1) -1;
        end
        jj = i/(k+1) - (ii+1)/(k+1);
        
        [x, y] = mapp_xy(nodes(ii+1), nodes(jj+1), elem, n); %%physical qudrature nodes
        
        
        [ux, uy] = compute_grad(n,k,LAMBDA,x,y); %%\nabla_x \lambda(x)
        
        xi1 = x - ux;
        xi2 = y - uy;
        
        test = 0;
        xiold = inf;
        if (xi1 >1)
            xiold = xi1;
            test = 1;
            xi1 = 1;
        elseif (xi1 <0)
            xiold = xi1;
            test = 1;
            xi1 = 0;
        end
        if(test)
            fprintf('Warning: xi1=%f\n', xiold);
        end
        
        test = 0;
        xiold = inf;
        if (xi2 >1)
            xiold = xi2;
            test = 1;
            xi2 = 1;
        elseif (xi2 <0)
            xiold = xi2;
            test = 1;
            xi2 = 0;
        end
        if(test)
            fprintf('Warning: xi2=%f\n', xiold);
        end
        
        PROJ(ie) = compute_sol(xi1, xi2, PHI, n, k);
        
    end
    kk = kk + Nloc;
end


end