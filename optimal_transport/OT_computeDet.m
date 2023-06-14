function DET = OT_computeDet(n,k,U, eig_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Compute strong determinant on quadrature points %%%%%%%%
%%%%%% if eig==1, compute det(I-\nabla^2 U)(x)
%%%%%% if eig==2, compute det(I+\nabla^2 U)(x)
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
        
        DET(ie) = OT_compute_det_xy(x,y,n,k,U, eig_m);
        
    end
    kk = kk + Nloc;
end


end