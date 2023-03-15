function plot_sol2(n,k,U)
%%% Plot DDLs solution of the obstacle problem %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Nloc  = k+1;                %local dimension vector
%gdim  = Nloc*n^2;              %global dimension matrix


[nodes, ~] = getWeightsNodes(k+1);


x = NaN(Nloc*n,1);
y = NaN(Nloc*n,1);
M = NaN(Nloc*n,Nloc*n);


kk = 0;
nn = 0;
for elem=1:n
    for i=1:k+1
        ie = i + kk;
        x(ie) = mapp_x(nodes(i), elem, n);
        y(ie) = mapp_y(nodes(i), elem+nn, n);
        %M(elem+ii, elem+jj) = U(ie);
    end
    kk = kk + Nloc;
    nn = nn + n;
end

[X,Y] = meshgrid(x,y);
for i = 1:length(x)
    for j=1:length(y)
        M(i,j) = compute_sol(x(i),y(j),U,n,k);
    end
end


        
%subplot(121)
f = figure(1);
f.Position(3:4) = [1100 600];
surf(X,Y,M)%,'EdgeColor','none')
colorbar
title('Sol Approx'); xlabel('x'); ylabel('y'); zlabel('u')


end