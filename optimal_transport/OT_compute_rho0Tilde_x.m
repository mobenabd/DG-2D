function [status, PROJ] = OT_compute_rho0Tilde_x(n,k, rho_1, Alpha_x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Compute rho_tilde_(x) = \rho_1 / det(\nabla^2 \alpha(x))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloc = (k+1)^2;
gdim = n^2*Nloc;
PROJ = zeros(gdim,1);

[nodes, ~] = getWeightsNodes(k+1);

status = true;

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
        
        [x, y] = mapp_xy(nodes(ii+1), nodes(jj+1), elem, n); %%physical quadrature node
        
        
        det = compute_det(x,y,n,k,Alpha_x, 0);
        
        if (det <=0)
            fprintf('Warning: det=%f\n', det);
            status = false;
            return
        end
        PROJ(ie) = rho_1(x,y)/det;
        

    end
    kk = kk + Nloc;
end


end
