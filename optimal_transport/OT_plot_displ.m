function OT_plot_displ(n,k,U,Phi_xi)

%%%%%%%%%% Plot contours %%%%%%%%%%
[x0, xN, ~] = getGlobal_x0N();

dx= (xN-x0)*0.01;

x = x0:dx:xN; y = x0:dx:xN;
[XX,YY] = meshgrid(x,y);
M = zeros(length(x),length(y));

for i=1:length(y)
    for j = 1:length(x)
        M(i,j) =compute_sol(x(j),y(i),U,n,k);
    end
end

%%%%%%%%%% Plot vector gradient %%%%%%%%%%

len = 40;
x = x0:(xN-x0)/len:xN-(xN-x0)/len; y = x0:(xN-x0)/len:xN-(xN-x0)/len;


Ux = NaN(len^2,1); Uy = NaN(len^2,1);
X  = NaN(len^2,1); Y  = NaN(len^2,1);

for i=1:len
    for j = 1:len        
        idx = (i-1)*len + j;
        [ux, uy] = compute_grad(n,k, Phi_xi, x(i), y(j));
        X(idx)  = x(i);
        Y(idx)  = y(j);
        Ux(idx) = ux;
        Uy(idx) = uy;  
    end
end


f = figure(3);
f.Position(3:4) = [700 500];


contourf(XX,YY,M,20,'EdgeColor','none')
%colormap jet
colorbar

hold on

%quiver(X,Y,Ux,Uy,0, 'Color','k')
quiver(X,Y,Ux,Uy, 'Color','k')
axis([x0 xN x0 xN])
xlabel('x'); ylabel('y'); zlabel('u')
hold off













end