function b = localRHS_edgeBoundaryNeumann(verHor,ort,num,n,k,fN)
%%% routine:  Neumann buondary condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n       : Global discretisation 
% k       : polynomial degree Q_k
% verHor  : =1 if vertical interior edge. =0 if horizontale interior edge.
% ort     : =1 if left/bottom boundary, 0 if right/top boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[~, ~, h] = getGlobal_x0N();
Nloc = (k+1)^2;
b = zeros(Nloc,1);

if (ort==1)
    ss = -1;
elseif (ort==0)
    ss = 1;
end


if (verHor ==1 )   %vertical edge on x=+/-1
    func = @(x,y) ss*fN(x);
elseif (verHor == 0)  %horizontal edge on y=+/-1
    func = @(x,y) ss*fN(y);
end



%[Jac,~, ~] = jacobian_elem(1,n) ; %(uniforme grid)
%detJac11 = Jac(1,1);
detJac11 = h/2;

uexct_loc  =  func_mapping(func,n,num);

for i=1:Nloc
    [phi, ~, ~] = basis(i,n,k);
    if (verHor ==1 )   %vertical edge on x=+/-1, normal = [+/-1 0];
        b(i) = detJac11* quadrature(@(y) phi(ss,y).*uexct_loc(ss,y),1,5);
        
    elseif (verHor == 0)  %horizontal edge on y=+/-1, normal = [0 +/-1];
        b(i) = detJac11* quadrature(@(x) phi(x,ss).*uexct_loc(x,ss),1,5);
    end
end


end