function [E, dE] = OT2_compute_energy(phi_xi, rho_1, Phi_xi, n)

k = getGlobal_k();

f1 = @(xi1,xi2) func1(xi1, xi2, phi_xi, n,k);
E =  0.5*integrate2(f1, n);


f2 = @(xi1,xi2) func2(xi1,xi2, phi_xi, rho_1, Phi_xi, n, k);
dE = integrate2(f2, n);

end


function val = func1(xi1, xi2, phi_xi, n,k)

[ux, uy] = compute_grad(n,k,phi_xi,xi1,xi2);

val = ux^2 + uy^2; 

end

function val = func2(xi1,xi2, phi_xi, rho_1, Phi_xi, n, k)

[ux, uy] = compute_grad(n,k,Phi_xi,xi1,xi2);
x = xi1 + ux;
y = xi2 + uy;

[uxi1, uxi2] = compute_grad(n,k,phi_xi,xi1,xi2);


[uxx, uyy, uxy] = compute_secondDerv(n,k,Phi_xi,xi1,xi2);

a11 = 1 + uxx;
a22 = 1 + uyy;
a12 = uxy;

% det_xi = (a11*a22 - a12^2);
% 
% coeff = 1/det_xi;

b11 = a22;
b22 = a11;
b12 = -a12;

%x^T * A^-1 * x
val = ( uxi1*(b11*uxi1 + b12*uxi2) + uxi2*(b12*uxi1 + b22*uxi2));%  *coeff;  


val = rho_1(x,y)*val;%*det_xi;


end