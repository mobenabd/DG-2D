function [G, Da, Db]= updateNewton_obstacle(n,k, G, U, obstacle_g)
%%%%%%%%% Update Newton vextor A*U-rhs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% using the obstacle g  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% elem: element number
% n   : global discretisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloc  = (k+1)^2;                %local dimension matrix
gdim  = Nloc*n^2;               %global dimension matrix
h = 1/n;

Da = zeros(gdim, gdim);
Db = zeros(gdim, gdim);

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
    
    
    %x_coord = [x1 x2 x3 x4];    y_coord = [y1 y2 y3 y4];
    if (k==1 || k==2 )
        x_coord =  [-1 1 1 -1 0 0 -1 1 0]; 
        y_coord = [-1 -1 1 1 -1 1 0 0 0];
    elseif (k==3)
        x_coord = [-1 1 1 -1 -1/3 1/3 1/3 -1/3 1 1 -1 -1 -1/3 1/3 1/3 -1/3]; 
        y_coord = [-1 -1 1 1 -1 -1 1 1 -1/3 1/3 1/3 -1/3 -1/3 -1/3 1/3 1/3];
    end
    
    %for now its works only for Nloc=4 
    for i=1:Nloc 
        ie = i + kk;
        [xi, yi] =  mapp_xy(x_coord(i), y_coord(i), elem, n);
        diffr = U(ie) - obstacle_g(xi, yi);
        %fprintf('elem=%i, i=%i x=%f, y=%f\n',elem,i,xi, yi);
        
        %diffr = U(ie) - obstacle_g(x_coord(i), y_coord(i));
        %fprintf('elem=%i, i=%i x=%f, y=%f\n',elem,i,x_coord(i), y_coord(i));
        
        if (diffr <= G(ie))
            Db(ie,ie) = 1;
            Da(ie,ie) = 0;
            G(ie) = diffr;
        else
            Da(ie,ie) = 1;
            Db(ie,ie) = 0;
        end
    end
 
    kk = kk + Nloc;
end


end