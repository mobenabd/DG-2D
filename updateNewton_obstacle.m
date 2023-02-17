function [G, Da, Db]= updateNewton_obstacle(n,k, G, U, g_projected, obstacle_g)
%%%%%%%%% Update Newton vector A*U-rhs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% using the obstacle g  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% this uses 'obstacle_g' function if it's avaible %%%%%%%%%%%%%%%%%
%%%%%%%%% otherwise: use the L2 projection 'g_projected' %%%%%%%%%%%%%%%%%%
% elem: element number
% n   : global discretisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloc  = (k+1)^2;                %local dimension matrix
gdim  = Nloc*n^2;               %global dimension matrix

Da = zeros(gdim, gdim);
Db = zeros(gdim, gdim);

[nodes, ~] = getWeightsNodes(k+1);


kk = 0;
for elem=1:n^2
    %for now its works only for Lagrange basis
    for i=1:Nloc
        ie = i + kk;
        diffr = U(ie) - g_projected(ie);
        
        if exist('obstacle_g', 'var')
            if mod(i, k+1) ==0
                ii = k;
            else
                ii =  mod(i, k+1) -1;
            end
            jj = i/(k+1) - (ii+1)/(k+1);
            [xi, yi] =  mapp_xy(nodes(ii+1), nodes(jj+1), elem, n);
            diffr = U(ie) - obstacle_g(xi, yi);
        end
        
        if (diffr <= G(ie))
            Db(ie,ie) = 1;
            %Da(ie,ie) = 0;
            G(ie) = diffr;
        else
            Da(ie,ie) = 1;
            %Db(ie,ie) = 0;
        end
    end
    
    kk = kk + Nloc;
end


end