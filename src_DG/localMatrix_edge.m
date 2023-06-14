function [M11, M22, M12, M21] = localMatrix_edge(verHor,n,k,eps,sigma, num, alpha_xy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% compute local matrix M_ij (integrals over egdes) %%%%%%%%%%%%%%%%%
% n        : Global discretisation 
% k        : polynomial degree  Q_k  
% eps      : symmetrization parameter
% sigma    : penelisation parameterc
% verHor   : =1 if vertical interior edge. =0 if horizontale interior edge.
% num      : 
% alpha_xy : diffusion function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('alpha_xy', 'var')
     % mapping diffusion function to refrence element (xi,nu) -> foF(xi, nu)) = f(x,y)
    alpha_mapped = func_mapping(alpha_xy,n,num);
else
    alpha_mapped = @(x,y) 1;
end

[~, ~, h] = getGlobal_x0N();
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
           M11(i,j) =  detJac11*quadrature(@(y) alpha_mapped(-1,y).*(-0.5*-dphjx(-1,y) .* phi(-1,y) ...
                       +(eps/2)* -dphix(-1,y) .* phj(-1,y)) ...
                       +(sigma/h^beta0) *phi(-1,y).*phj(-1,y),1);
                   
           M22(i,j) = detJac11* quadrature(@(y) alpha_mapped(-1,y).*(0.5 *-dphjx(1,y) .* phi(1,y) ...
                       -(eps/2)* -dphix(1,y) .* phj(1,y)) ...
                       +(sigma/h^beta0) *phi(1,y).*phj(1,y),1);
                   
           M12(i,j)  = detJac11* quadrature(@(y) alpha_mapped(-1,y).*(-0.5 *-dphjx(1,y) .* phi(-1,y) ...
                       -(eps/2)* -dphix(-1,y) .* phj(1,y)) ...
                       -(sigma/h^beta0) *phj(1,y).*phi(-1,y),1);
                   
           M21(i,j)  = detJac11* quadrature(@(y) alpha_mapped(-1,y).*(0.5*-dphjx(-1,y) .* phi(1,y) ...
                       +(eps/2)* -dphix(1,y) .* phj(-1,y)) ...
                       - (sigma/h^beta0) * phj(-1,y).*phi(1,y),1);
                   
        elseif (verHor == 0)  %horizontal edge on y=-1, normal = [0 -1];
           M11(i,j) = detJac11*quadrature(@(x) alpha_mapped(x,-1).*(-0.5 *-dphjy(x,-1) .* phi(x,-1) ...
                       +  (eps/2)*-dphiy(x,-1) .* phj(x,-1)) ...
                       +  (sigma/h^beta0) *phi(x,-1).*phj(x,-1),1);
                   
           M22(i,j) = detJac11* quadrature(@(x) alpha_mapped(x,-1).*(0.5*-dphjy(x,1) .* phi(x,1) ...
                       +  -(eps/2)*-dphiy(x,1) .* phj(x,1)) ...
                       +  (sigma/h^beta0)*phi(x,1).*phj(x,1),1);
                   
           M12(i,j) = detJac11*quadrature(@(x) alpha_mapped(x,-1).*(-0.5*-dphjy(x,1) .* phi(x,-1) ...
                        -(eps/2)*-dphiy(x,-1) .* phj(x,1)) ...
                        -(sigma/h^beta0)*phj(x,1).*phi(x,-1),1);
                                
           M21(i,j) = detJac11* quadrature(@(x) alpha_mapped(x,-1).*(0.5 *-dphjy(x,-1) .* phi(x,1) ...
                       + (eps/2)*-dphiy(x,1) .* phj(x,-1)) ...
                       -(sigma/h^beta0)*phj(x,-1).*phi(x,1),1);
        end 
    end
end


end