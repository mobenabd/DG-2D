function PhyNodes = getPhysicalNodes(n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Compute quadrature physical nodes coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = getGlobal_k();
Nloc = (k+1)^2;
gdim = n^2*Nloc;
PhyNodes = NaN(gdim,2);

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
        
        [x, y] = mapp_xy(nodes(ii+1), nodes(jj+1), elem, n); %%physical quadrature node
           
        PhyNodes(ie,1) = x;   PhyNodes(ie,2) = y;
        
    end
    kk = kk + Nloc;
end


end
