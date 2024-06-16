function plot_sol3(n,k, U)
%%% Plot solution of the obstacle problem %%%%%%%%%%%%%
%% for 2D plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x0, xN, ~] = getGlobal_x0N();
%dx= 0.01;
%dx= dx/2;
%dx = 4*dx;

dx = (xN-x0)*0.02;
x = x0:dx:xN; y = x0:dx:xN;
[X,Y] = meshgrid(x,y);
M = zeros(length(x),length(y));
for i=1:length(y)
    for j = 1:length(x)
        M(i,j) =compute_sol(x(j),y(i),U,n,k);
    end
end

f = figure(3);
f.Position(3:4) = [700 500]; 
%subplot(122)
%surf(X,Y,M)%,'FaceAlpha','0.')%'FaceLighting','gouraud','EdgeColor','black')
surf(X,Y,M, 'EdgeColor','none', 'FaceColor','interp', 'FaceLighting','gouraud')
view(2)
colormap hsv
 colorbar
h=colorbar;
t=get(h,'Limits');
T=linspace(t(1),t(2),5);
T = [T 0];
T = sort(T);
set(h,'Ticks',T);
TL=arrayfun(@(x) sprintf('%.2f',x),T,'un',0);
set(h,'TickLabels',TL);

xlabel('x'); ylabel('y'); 
end