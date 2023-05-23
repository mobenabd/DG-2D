function MESH = OT_test(n,k, M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Compute x knowing xi on quad nodes %%
%%%% compare with exact formula %%%%%%%%%%

Nloc = (k+1)^2;
gdim = n^2*Nloc;
MESH = zeros(gdim,4);

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
        
        [xi_1, xi_2] = mapp_xy(nodes(ii+1), nodes(jj+1), elem, n); %%physical quadrature node
        
        %[x_app,y_app] = OT_computeX_xi(xi_1,xi_2, M);      %%with interpolation
        %[x_app,y_app] = OT_computeX_xi2(xi_1,xi_2, M, n,k);  %%with relaxed Newton, M=lambda_x
        [x_app,y_app] = OT_computeX_xi3(xi_1,xi_2, M, n,k); %%with Gauss-Newton, M=lambda_x
        
        [x, y] = OT_Xxi2(xi_1, xi_2);
        
        MESH(ie,1) = x; MESH(ie,2) = x_app;
        MESH(ie,3) = y; MESH(ie,4) = y_app;
        

    end
    kk = kk + Nloc;
end


end
