function U = regularize(n,k,U, DET)
%%% Regularize DDLs after convexification %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Nloc = (k+1)^2;

[nodes, ~] = getWeightsNodes(k+1);

beta =  @(x,y) 0.5*((x-0.5).^2 + (y-0.5).^2);
BETA =  computeDirectProjection(n,k,beta);

kk = 0;
for elem=1:n^2
    for i=1:Nloc
        ie = i + kk;
        
        if (DET(ie) <=0)
            if mod(i, k+1) ==0
                ii = k;
            else
                ii =  mod(i, k+1) -1;
            end
            jj = i/(k+1) - (ii+1)/(k+1);
            
            [x, y] = mapp_xy(nodes(ii+1), nodes(jj+1), elem, n);
            
            %lambda = 0.5*compute_det(x,y,n,k,U, 2); %largest eigenvalue
            
            %lap = compute_lap2(n,k,U,x,y);
            lambda =1;
            %U(ie) = U(ie) + 0.5*lambda*(x^2 + y^2);
            lambda = 2*U(ie)/BETA(ie);
            U(ie) = lambda/2 * (x^2 + y^2); 
        end
        
    end
    kk = kk + Nloc;
end


end