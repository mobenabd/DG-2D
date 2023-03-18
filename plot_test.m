function plot_test(n,k, U, g, uexct)
%%% Plot solution of the obstacle problem %%%%%%%%%%%%%
%%% plt=1 for 1D plot on y=0.5 and exact solution exists
%%% plt=2 for 2D plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%dx = 1/n;

dx= 0.01;


% if  ~exist('uexct','var')
%     uexct = g;
% end


x = 0:dx:1;
M = zeros(length(x),1);
Mexct = zeros(length(x),1);
for i = 1:length(x)
    M(i)     = compute_sol(x(i),x(i),U,n,k);
end
    
if  (~exist('uexct','var') && exist('g','var'))
    for i = 1:length(x)
        Mexct(i) = g(x(i),x(i),0);
    end
elseif ( exist('uexct','var'))
    for i = 1:length(x)
        Mexct(i) = uexct(x(i),x(i),0);
    end
end



f = figure (1);
%f.Position
%f.Position(3:4) = [1100 600];
%subplot(122)
if exist('g','var')
    plot(x, M ,'DisplayName', 'Sol approx', 'Marker','.','MarkerSize', 10,'LineWidth',1)
else
    plot(x, M ,'DisplayName', 'Obstacle (projected)', 'Marker','.','MarkerSize', 10,'LineWidth',1, 'Color', "#D95319")
end
hold on
if exist('uexct','var')
    plot(x, Mexct, 'LineWidth',1 , 'DisplayName', 'Sol exacte')
elseif (~exist('uexct','var') &&  exist('g','var'))
   plot(x, Mexct,'LineWidth',1 ,'DisplayName', 'Sol exacte du P-O')
end
legend('-DynamicLegend')
%axis([0 1 -1.5 1.5])
title('y=x')
hold off
    
return
%dx= 1/n;
%dx= dx/2;
dx = 4*dx;
x = 0:dx:1;
y = 0:dx:1;
[X,Y] = meshgrid(x,y);
M = zeros(length(x),length(y));
for i = 1:length(x)
    for j=1:length(y)
        M(i,j) = compute_sol(x(i),y(j),U,n,k);
    end
end

        
% subplot(121)
% mesh(X,Y,M)%,'EdgeColor','none')
% colorbar
% title('Sol Approx'); xlabel('x'); ylabel('y'); zlabel('u')

 %pause(5);
end