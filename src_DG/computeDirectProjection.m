function PROJ = computeDirectProjection(n,k,f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Cpmpute L2 projection directly with values of quadrature nodes %%%%
%%%%%% this works only for nodale basis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        
        [x, y] = mapp_xy(nodes(ii+1), nodes(jj+1), elem, n);
        
        PROJ(ie) = f(x,y);
    end
    kk = kk + Nloc;
end


end