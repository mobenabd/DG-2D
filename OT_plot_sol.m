function OT_plot_sol(n,k,U,numFig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('numFig', 'var')
    numFig = 2;
end 


dx= 0.01;




f = figure(numFig);
f.Position(3:4) = [700 500];

%dx = 4*dx;
x = 0:dx:1;
y = 0:dx:1;
[X,Y] = meshgrid(x,y);
M = zeros(length(x),length(y));

for i=1:length(y)
    for j = 1:length(x)
        M(i,j) =compute_sol(x(j),y(i),U,n,k);
    end
end



%surf(X,Y,M)%,'EdgeColor','none')
mesh(X,Y,M)
%contourf(X,Y,M,15,'EdgeColor','none')
colorbar
xlabel('x'); ylabel('y'); zlabel('u')

%hold on
%pause(5);
end