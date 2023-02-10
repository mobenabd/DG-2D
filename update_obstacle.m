function U= update_obstacle(n,k, U, obstacle_g)
%%%%%%%%% Update Lagrange ddls  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% using the obstacle g  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% elem: element number
% n   : global discretisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloc  = (k+1)^2;                %local dimension matrix
%gdim  = Nloc*n^2;               %global dimension matrix
h = 1/n;


kk = 0;
for elem=1:n^2
    %indice (i,j) de l'element elm
    if (mod(elem,n) ~= 0)
        i_elem = fix(elem/n) + 1;   j_elem = mod(elem,n);
    else % (j=n)
        i_elem = elem/n;    j_elem = n;
    end
    %coordonnée cartesienne des sommets
    x1 = (j_elem-1)*h;  y1 = (i_elem-1)*h;
    x2 = x1 + h;    y2 = y1;
    x4 = x1;        y4 = y1 + h;
    x3 = x2;        y3 = y4;
    %fprintf('x1=%f, y1=%f, x2=%f, y2=%f, x3=%f, y3=%f, x4=%f, y4=%f\n',x1,y1,x2,y2,x3,y3,x4,y4);
    x_coord = [x1 x2 x3 x4];    y_coord = [y1 y2 y3 y4];
    
    %for now its works only for Nloc=4 
    for i=1:Nloc 
        ie = i + kk;
        if (U(ie)  <= obstacle_g(x_coord(i), y_coord(i)))
            U(ie) =  obstacle_g(x_coord(i), y_coord(i));
        end
    end
 
    kk = kk + Nloc;
end


end