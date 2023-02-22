function plot_det(n,k,U, Lap, det_exct)
%%% Plot solution determiant %%%%%%%%%%%%%%%%%%%%
%%% plt=1 for 1D plot on y=0.5 and exact solution exists
%%% plt=2 for 2D plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%dx = 1/n;

dx= 0.01;



x = 0:dx:1;
M = zeros(length(x),1);
for i = 1:length(x)
    M(i) =  compute_det(n,k,Lap,U,x(i),0.5);
end

cc = 2;
if  exist('det_exct','var')
    Mexct = zeros(length(x),1);
    for i = 1:length(x)
        Mexct(i) = det_exct(x(i),0.5,0);
    end
    cc = 3;
end



f2 = figure(2);
f2.Position(3:4) = [1100 600];
subplot(1,cc,1)
plot(x, M ,'DisplayName', 'det approx', 'Marker','.','MarkerSize', 10,'LineWidth',1)
if  exist('det_exct','var')
    hold on
    plot(x, Mexct, 'LineWidth',1 , 'DisplayName', 'det exacte')
end
legend('-DynamicLegend')
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
        M(i,j) =  compute_det(n,k,Lap,U,x(i),y(j));
    end
end


if  exist('det_exct','var')
    Mexct = zeros(length(x),length(y));
    for i = 1:length(x)
        for j=1:length(y)
            Mexct(i,j) = det_exct(x(i),y(j),0);
        end
    end
end



subplot(1,cc,2)
surf(X,Y,M, 'EdgeColor','none')%, 'FaceAlpha',0.5);
title('det approx')

if  exist('det_exct','var')
    subplot(1,cc,3)
    surf(X,Y,Mexct, 'EdgeColor','none')%, 'FaceAlpha',0.5)
    title('det exacte')
    colorbar
    
    
    err2 = dx*sqrt(sum(sum((Mexct-M).^2)));
    fprintf('|e|_2 (det) = %.10f\n',err2);
    
end


end