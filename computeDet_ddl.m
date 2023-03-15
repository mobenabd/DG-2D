function DET = computeDet_ddl(n,k,U, eig_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Compute strong determinant on quadrature points %%%%%%%%
%%%%%% if eig==1, compute the largest eigenvalue in each quadrature point
%%%%%% instead of the determinant %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloc = (k+1)^2;
gdim = n^2*Nloc;
DET = zeros(gdim,1);

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
        
        DET(ie) = compute_det(x,y,n,k,U, eig_m);
        
    end
    kk = kk + Nloc;
end


end