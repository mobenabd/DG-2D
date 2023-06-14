function M = OT_Xxi(n,k,LAMBDA)  %%M = (x,y, xi_1, xi_2)

Nloc = (k+1)^2;
gdim = n^2*Nloc;
M = zeros(gdim,4);

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
        
        
        [ux, uy] = compute_grad(n,k,LAMBDA,x,y); %%\nabla_x \lambda(x)
        
        M(ie,1) = x; M(ie,2) = y;
        M(ie,3) = x -  ux; M(ie,4) = y -  uy; 
    end
    kk = kk + Nloc;
end


end