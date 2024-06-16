function plot_solHeat(n,k,t, U, uexct)
%%% Plot solution of Heat equation %%%%%%%%%%%%%%%%%%%%%%
%%% plt=1 for 1D plot at y=0.6 and exact solution exists
%%% plt=2 for 2D plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




[x0, xN, ~] = getGlobal_x0N();
dx= (xN-x0)*0.01/2;
%%% Exact solution known %%%%%%%
if  exist('uexct','var')
    x = x0:dx:xN;
    M = zeros(length(x),1);
    Mexct = zeros(length(x),1);
    for i = 1:length(x)
        M(i)     = compute_sol(x(i),0.6,U,n,k);
        Mexct(i) = uexct(x(i),0.6,t);
    end
end

f = figure (1);
%f.Position
f.Position(3:4) = [1100 600];
if  exist('uexct','var')
    subplot(122)
    plot(x, M ,'DisplayName', 'Sol approx', 'Marker','.','MarkerSize', 10,'LineWidth',1)
    hold on
    if exist('uexct','var')
        plot(x, Mexct, 'LineWidth',1 , 'DisplayName', 'Sol exact')
    else
        plot(x, Mexct,'LineWidth',1 ,'DisplayName', 'Obstacle')
    end
    legend('-DynamicLegend')
    %axis([0 1 -1.5 1.5])
    title('at y=0.6')
    hold off
end

dx = 2*dx;
x = x0:dx:xN;
y = x0:dx:xN;
[X,Y] = meshgrid(x,y);
M = zeros(length(x),length(y));
for i=1:length(y)
    for j = 1:length(x)
        M(i,j) =compute_sol(x(j),y(i),U,n,k);
    end
end

if  exist('uexct','var')
    subplot(121)
end
surf(X,Y,M);%,'EdgeColor','none')
colorbar
title('Sol Approx'); xlabel('x'); ylabel('y'); zlabel('u')
%pause(5);
end
