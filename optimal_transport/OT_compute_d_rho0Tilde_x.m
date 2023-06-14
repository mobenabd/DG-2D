function PROJ = OT_compute_d_rho0Tilde_x(n,k, rho1_proj, dlambda_x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Compute  Div \Tilde{\rho_0} (x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloc = (k+1)^2;
gdim = n^2*Nloc;
PROJ = zeros(gdim,1);
[nodes, ~] = getWeightsNodes(k+1);



GRAD_dlambda = computeGrad_ddl(n,k,dlambda_x);

U = rho1_proj.*GRAD_dlambda;


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
        
        [x, y] = mapp_xy(nodes(ii+1), nodes(jj+1), elem, n);
        
        val = compute_div(n,k,U,x,y);
        
        
        PROJ(ie) = val;
        
    end
    kk = kk + Nloc;
end

end
