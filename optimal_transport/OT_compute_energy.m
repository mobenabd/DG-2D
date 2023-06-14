function [E, dE] = OT_compute_energy(phi_xi, rho_1, lambda_x, n)

k = getGlobal_k();
f1 = @(x,y) func1(x,y, phi_xi, n,k);
f2 = @(x,y) func2(x,y, phi_xi, rho_1, lambda_x, n,k);

E =  0.5*integrate2(f1, n);

dE = integrate2(f2, n);

end


function val = func1(x,y, phi_xi, n,k)

[ux, uy] = compute_grad(n,k,phi_xi,x,y);

val = ux^2 + uy^2; 

end

function val = func2(x,y, phi_xi, rho_1, lambda_x, n,k)

[ux, uy] = compute_grad(n,k,lambda_x,x,y);
xi1 = x - ux;
xi2 = y - uy;

[ux, uy] = compute_grad(n,k,phi_xi,xi1,xi2);

%val = rho_1(x,y) * (ux^2 + uy^2); 



%%%%%%%% the right way %%%%%%

[uxx, uyy, uxy] = compute_secondDerv(n,k,lambda_x,x,y);


val = rho_1(x,y) * ( (ux*(1-uxx) + uy*(-uxy))^2 + (ux*(-uxy) + uy*(1-uyy) )^2 ); 


end