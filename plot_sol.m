function plot_sol(n,k,t, U, g, uexct)
%%% Plot solution of the obstacle problem %%%%%%%%%%%%%
%%% plt=1 for 1D plot on y=0.5 and exact solution exists
%%% plt=2 for 2D plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%dx = 1/n;

dx= 0.01;


if  ~exist('uexct','var')
    uexct = g;
end


x = 0:dx:1;
M = zeros(length(x),1);
Mexct = zeros(length(x),1);
for i = 1:length(x)
    M(i)     = compute_sol(x(i),0.5,U,n,k);
    Mexct(i) = uexct(x(i),0.5,t);
end
    

f = figure (1);
%f.Position
f.Position(3:4) = [1100 600];
subplot(122)
plot(x, M ,'DisplayName', 'Sol approx', 'Marker','.','MarkerSize', 10,'LineWidth',1)
hold on
if exist('uexct','var')
    plot(x, Mexct, 'LineWidth',1 , 'DisplayName', 'Sol exacte')
else
   plot(x, Mexct,'LineWidth',1 ,'DisplayName', 'Obstacle')
end
legend('-DynamicLegend')
%axis([0 1 -1.5 1.5])
title('at y=0.5')
hold off
    

%dx= 1/n;
%dx= dx/2;
dx = dx*2;
x = 0:dx:1;
y = 0:dx:1;
[X,Y] = meshgrid(x,y);
M = zeros(length(x),length(y));
for i = 1:length(x)
    for j=1:length(y)
        M(i,j) = compute_sol(x(i),y(j),U,n,k);
    end
end

        
subplot(121)
surf(X,Y,M,'EdgeColor','none')
colorbar
title('Sol Approx'); xlabel('x'); ylabel('y'); zlabel('u')

 %pause(5);
end