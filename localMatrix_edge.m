function [M11, M22, M12, M21] = localMatrix_edge(verHor,n,k,eps,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% compute local matrix M_ij (integrals over egdes) %%%%%%%%%%%%%%%%%
% n       : Global discretisation 
% k       : polynomial degree  Q_k  
% eps     : symmetrization parameter
% sigma   : penelisation parameterc
% verHor  : =1 if vertical interior edge. =0 if horizontale interior edge.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h     = 1/n;
Nloc  = (k+1)^2;
beta0 = 1;     %supper-penalization =1 (for now)


M11 = zeros(Nloc,Nloc);
M22 = zeros(Nloc,Nloc);
M12 = zeros(Nloc,Nloc);
M21 = zeros(Nloc,Nloc);


%[Jac,~, ~] = jacobian_elem(1,n) ; %(interior cells && uniforme grid)
%detJac11 = Jac(1,1);
detJac11 = h/2;

for i=1:Nloc
    for j=1:Nloc
        [phi, dphix, dphiy] = basis(i,n,k);
        [phj, dphjx, dphjy] = basis(j,n,k);
        
        if (verHor ==1)   %vertical edge on x=-1, normal = [-1 0];
           M11(i,j) = -detJac11*0.5 * quadrature(@(y) -dphjx(-1,y) .* phi(-1,y),1) ...
                       +detJac11*(eps/2)* quadrature(@(y) -dphix(-1,y) .* phj(-1,y),1) ...
                       +detJac11*(sigma/h^beta0) *quadrature(@(y) phi(-1,y).*phj(-1,y),1);
                   
           M22(i,j) = detJac11*0.5 * quadrature(@(y) -dphjx(1,y) .* phi(1,y),1) ...
                       -detJac11*(eps/2)* quadrature(@(y) -dphix(1,y) .* phj(1,y),1) ...
                       + detJac11*(sigma/h^beta0) *quadrature(@(y) phi(1,y).*phj(1,y),1);
                   
           M12(i,j)  = -detJac11*0.5 * quadrature(@(y) -dphjx(1,y) .* phi(-1,y),1) ...
                       -detJac11*(eps/2)* quadrature(@(y) -dphix(-1,y) .* phj(1,y),1) ...
                       - detJac11*(sigma/h^beta0) *quadrature(@(y) phj(1,y).*phi(-1,y),1);
                   
           M21(i,j)  = detJac11*0.5 * quadrature(@(y) -dphjx(-1,y) .* phi(1,y),1) ...
                       +detJac11*(eps/2)* quadrature(@(y) -dphix(1,y) .* phj(-1,y),1) ...
                       - detJac11*(sigma/h^beta0) *quadrature(@(y) phj(-1,y).*phi(1,y),1);
                   
        elseif (verHor == 0)  %horizontal edge on y=-1, normal = [0 -1];
           M11(i,j) = -detJac11*0.5 * quadrature(@(x) -dphjy(x,-1) .* phi(x,-1),1) ...
                       +detJac11*(eps/2)* quadrature(@(x) -dphiy(x,-1) .* phj(x,-1),1) ...
                       + detJac11*(sigma/h^beta0) *quadrature(@(x) phi(x,-1).*phj(x,-1),1);
                   
           M22(i,j) = detJac11*0.5 * quadrature(@(x) -dphjy(x,1) .* phi(x,1),1) ...
                       -detJac11*(eps/2)* quadrature(@(x) -dphiy(x,1) .* phj(x,1),1) ...
                       + detJac11*(sigma/h^beta0) *quadrature(@(x) phi(x,1).*phj(x,1),1);
                   
           M12(i,j) = -detJac11*0.5 * quadrature(@(x) -dphjy(x,1) .* phi(x,-1),1) ...
                       -detJac11*(eps/2)* quadrature(@(x) -dphiy(x,-1) .* phj(x,1),1) ...
                       -detJac11*(sigma/h^beta0) *quadrature(@(x) phj(x,1).*phi(x,-1),1);
                  
           M21(i,j) = detJac11*0.5 * quadrature(@(x) -dphjy(x,-1) .* phi(x,1),1) ...
                       +detJac11*(eps/2)*quadrature(@(x) -dphiy(x,1) .* phj(x,-1),1) ...
                       -detJac11*(sigma/h^beta0) *quadrature(@(x) phj(x,-1).*phi(x,1),1);
        end 
    end
end


end