function [Jac,bE, det_Jac] = jacobian_elem(elem,n)
%%%%%%%%% compute jacobian matrix and determinant  %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Affine map function define as F(x) = Jac*x+bE  %%%%%%%%%%%%%%%%%%
% elem: element number 
% n   : global discretisation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x0, ~, h] = getGlobal_x0N();


%indice (i,j) de l'element elm
if (mod(elem,n) ~= 0)
    i_elem = fix(elem/n) + 1;
    j_elem = mod(elem,n);
else % (j=n)
    i_elem = elem/n;
    j_elem = n;
end

%coordonnée cartesienne des sommets
x1 =  x0 + (j_elem-1)*h; y1 =  x0 + (i_elem-1)*h;
x2 = x1 + h;    y2 = y1;
x4 = x1;        y4 = y1 + h;     


%fprintf('x1=%f, y1=%f, x2=%f, y2=%f, x4=%f, y4=%f\n',x1,y1,x2,y2,x4,y4);

%la matrice jacobienne
Jac      = zeros(2,2);    
Jac(1,1) = x2 - x1;
Jac(1,2) = x4 - x1;
Jac(2,1) = y2 - y1;
Jac(2,2) = y4 - y1;
Jac = 0.5 * Jac;

%determinant
det_Jac = det(Jac);

% vecteur à l'origine
bE = zeros(2,1);
bE(1) = x2 + x4;
bE(2) = y2 + y4;
bE = 0.5 * bE;

end