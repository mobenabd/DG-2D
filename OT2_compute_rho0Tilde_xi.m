function PROJ = OT2_compute_rho0Tilde_xi(n,k, rho_1, Phi_xi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Compute rho_tilde_(xi) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        
        [xi1, xi2] = mapp_xy(nodes(ii+1), nodes(jj+1), elem, n); %%physical quadrature node
        
        
        [ux, uy] = compute_grad(n,k, Phi_xi,xi1,xi2);
        
        x = xi1 + ux;
        y = xi2 + uy;
         
        
        PROJ(ie) = rho_1(x,y); %*det_xi vectorized

    end
    kk = kk + Nloc;
end

%%% positive density
PROJ(PROJ<0) = 0;
end
