function GRAD = computeGrad_ddl(n,k,U)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Compute strong gradient on quadrature points %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Nloc = (k+1)^2;
gdim = n^2*Nloc;
GRAD = zeros(gdim,2);

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
        
        [x, y] = mapp_xy(nodes(ii+1), nodes(jj+1), elem, n);
        
        [ux, uy] = compute_grad(n,k,U,x,y);
        
        %[ux, uy, ~] = compute_secondDerv(n,k,U,x,y);
        
        GRAD(ie,1) = ux;
        GRAD(ie,2) = uy;
        
    end
    kk = kk + Nloc;
end


end