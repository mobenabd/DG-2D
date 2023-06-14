function plot_sol3(n,k, U, U0)
%%% Plot solution of the obstacle problem %%%%%%%%%%%%%
%%% plt=1 for 1D plot DDLs
%%% plt=2 for 2D plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



dx= 0.01;

[x0, xN, ~] = getGlobal_x0N();

% if  ~exist('uexct','var')
%     uexct = g;
% end


x = x0:dx:xN;
M = zeros(length(x),1);
Mexct = U0;
for i = 1:length(x)
    M(i)     = compute_sol(x(i),x(i),U,n,k);
end
    


f = figure (1);
%f.Position
f.Position(3:4) = [1100 600];
subplot(122)
plot(x, M ,'DisplayName', 'Sol approx', 'Marker','.','MarkerSize', 10,'LineWidth',1)
hold on

plot(x, Mexct,'DisplayName', 'Obstacle (projected)', 'Marker','.','MarkerSize', 10,'LineWidth',1)

legend('-DynamicLegend')
%axis([0 1 -1.5 1.5])
title('y=x')
hold off
    


%dx= dx/2;
dx = 4*dx;
x = x0:dx:xN;
y = x0:dx:xN;
[X,Y] = meshgrid(x,y);
M = zeros(length(x),length(y));
for i=1:length(y)
    for j = 1:length(x)
        M(i,j) =compute_sol(x(j),y(i),U,n,k);
    end
end

        
subplot(121)
mesh(X,Y,M)%,'EdgeColor','none')
colorbar
title('Sol Approx'); xlabel('x'); ylabel('y'); zlabel('u')

 %pause(5);
end