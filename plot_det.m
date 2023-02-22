function plot_det(n,k,U, Lap, uexct)
%%% Plot solution determiant %%%%%%%%%%%%%%%%%%%%
%%% plt=1 for 1D plot on y=0.5 and exact solution exists
%%% plt=2 for 2D plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%dx = 1/n;

dx= 0.01;



x = 0:dx:1;
M = zeros(length(x),1);
Mexct = zeros(length(x),1);
for i = 1:length(x)
    %M(i)     = compute_sol(x(i),0.5,U,n,k);
    M(i) =  compute_det(n,k,Lap,U,x(i),0.5);
    Mexct(i) = uexct(x(i),0.5,0);
end
    

figure
subplot(131)
plot(x, M ,'DisplayName', 'det approx', 'Marker','.','MarkerSize', 10,'LineWidth',1)
hold on
plot(x, Mexct, 'LineWidth',1 , 'DisplayName', 'det exacte')

legend('-DynamicLegend')
%axis([0 1 -1.5 1.5])
title('at y=0.5')
hold off
    

%dx= 1/n;
dx= dx*2;
x = 0:dx:1;
y = 0:dx:1;
[X,Y] = meshgrid(x,y);
M = zeros(length(x),length(y));
Mexct = zeros(length(x),length(y));
for i = 1:length(x)
    for j=1:length(y)
        %M(i,j) = compute_sol(x(i),y(j),U,n,k);
        M(i,j) =  compute_det(n,k,Lap,U,x(i),y(j));
        Mexct(i,j) = uexct(x(i),y(j),0);
    end
end

        
subplot(132)
surf(X,Y,Mexct, 'EdgeColor','none')%, 'FaceAlpha',0.5)
title('det exacte')
subplot(133)
surf(X,Y,M, 'EdgeColor','none')%, 'FaceAlpha',0.5);
title('det approx')
colorbar


end