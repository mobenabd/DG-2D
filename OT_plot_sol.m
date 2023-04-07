function OT_plot_sol(n,k,U,g)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



dx= 0.01;




f = figure (1);
f.Position(3:4) = [1100 600];

%dx = 4*dx;
x = 0:dx:1;
y = 0:dx:1;
[X,Y] = meshgrid(x,y);
M = zeros(length(x),length(y));
if exist('g','var')
    for i = 1:length(x)
        for j=1:length(y)
            M(i,j) = g(x(i),y(j));
        end
    end
else
    for i = 1:length(x)
        for j=1:length(y)
            M(i,j) =compute_sol(x(i),y(j),U,n,k);
        end
    end
end

   
mesh(X,Y,M)%,'EdgeColor','none')
colorbar
xlabel('x'); ylabel('y'); zlabel('u')

%hold on
%pause(5);
end