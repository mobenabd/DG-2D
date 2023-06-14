function U0 = plot_sol(n,k, U, g, uexct)
%%% Plot solution of the obstacle problem %%%%%%%%%%%%%
%%% plt=1 for 1D plot on y=0.5 and exact solution exists
%%% plt=2 for 2D plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[x0, xN, ~] = getGlobal_x0N();
dx= (xN-x0)*0.01;


% if  ~exist('uexct','var')
%     uexct = g;
% end


x = x0:dx:xN;
M = zeros(length(x),1);
Mexct = zeros(length(x),1);
for i = 1:length(x)
    M(i)     = compute_sol(x(i),x(i),U,n,k);
end
    
if  (~exist('uexct','var') && exist('g','var'))
    for i = 1:length(x)
        Mexct(i) = g(x(i),x(i),0);
        %compute_sol(x(i),x(i),g,n,k);
    end
elseif ( exist('uexct','var'))
    for i = 1:length(x)
        Mexct(i) = uexct(x(i),x(i),0);
    end
end

U0 = M;

f = figure (1);
%f.Position
f.Position(3:4) = [1100 600];
subplot(122)
if exist('g','var')
    plot(x, M ,'DisplayName', 'Sol approx', 'Marker','.','MarkerSize', 10,'LineWidth',1)
else
    plot(x, M ,'DisplayName', 'Obstacle (projected)', 'Marker','.','MarkerSize', 10,'LineWidth',1, 'Color', "#D95319")
end
hold on
if exist('uexct','var')
    plot(x, Mexct, 'LineWidth',1 , 'DisplayName', 'Sol exacte')
elseif (~exist('uexct','var') &&  exist('g','var'))
   plot(x, Mexct,'LineWidth',1 ,'DisplayName', 'Obstacle')
end
legend('-DynamicLegend')
%axis([0 1 -1.5 1.5])
title('y=x')
hold off
    


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