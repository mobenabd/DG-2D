function [PROJ, X0new] = OT_compute_rho0Tilde_xi(n,k, rho0_x, M, X0old)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Compute rho_tilde_(xi) knowing rho_tilde(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloc = (k+1)^2;
gdim = n^2*Nloc;
PROJ = zeros(gdim,1);
X0new = zeros(gdim,2);

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
        
        [xi1, xi2] = mapp_xy(nodes(ii+1), nodes(jj+1), elem, n); %%physical quadrature node
        
        %[x, y] = OT_computeX_xi(xi1,xi2, M);       % with interpolation
        
        %[x,y] = OT_computeX_xi2(xi1,xi2, M, n,k);  % with relaxed Newton, M=lambda_x
        
        %[x, y] = OT_Xxi2(xi1, xi2);  %%for lambda exact known
        
        
        %%%% Default option
        x0 = X0old(ie,1); y0 = X0old(ie,2); 
        [x, y] = OT_computeX_xi3(x0, y0, xi1, xi2, M, n, k);  % with Guass-Newton, M=lambda_x
       
         
        test = 0;
        xold = inf;
        if (x >1)
            test = 1;
            xold = x;
            x = 1;
        elseif (x <0)
            test = 1;
            xold = x;
            x = 0;
        end
        if(test)
            fprintf('WARNING X(xi): xold=%f\n', xold);
        end
        
        test = 0;
        yold = inf;
        if (y >1)
            test = 1;
            yold = y;
            y = 1;
        elseif (y <0)
            test = 1;
            yold = y;
            y = 0;
        end
        if(test)
            fprintf('WARNING X(xi): yold=%f\n', yold);
        end
        
        X0new(ie,1) = x;   X0new(ie,2) = y;
        
        PROJ(ie) = compute_sol(x,y,rho0_x,n,k);

    end
    kk = kk + Nloc;
end

%%% positive density
PROJ(PROJ<0) = 0;
end
