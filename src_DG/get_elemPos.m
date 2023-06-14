function E = get_elemPos(x,y,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get element number at position (x,y) in the domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x0, xN, h] = getGlobal_x0N();
j = floor((x-x0)/h)+1;
i = floor((y-x0)/h)+1;
if (x==xN)
    j = j-1;
end
if (y==xN)
    i = i-1;
end

E = (i-1)*n+j;
end