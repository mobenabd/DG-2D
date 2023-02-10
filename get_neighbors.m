function [E1, E2] = get_neighbors(N1,verHor,n)
% E1,E2   : face neighbors 
% n       : Global discretisation 
% verHor  : =1 if vertical interior edge. =0 if horizontale interior edge.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% get (i,j) index
[i, j] = getIndex(N1,n);
elem = (i-1)*n + j;   %element number having N1 as left bottom vertix

% define which edge (vertical or horizontale)
E1 = elem;
if (verHor==1)      %vertical, normal = [-1 0];
    E2 = elem -1;
elseif (verHor==0)  %horizontal, normal = [0 -1];
    E2 = E1  -n;
end

end