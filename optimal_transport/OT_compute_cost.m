function cost = OT_compute_cost(rho_0, Phi_xi, n)
%%%%% compute transport cost (L2 wasserstein distance) %%%%%%%%%%%
k = getGlobal_k();
f1 = @(xi1,xi2) func1(xi1,xi2, rho_0, Phi_xi, n,k);


cost =  integrate2(f1, n);


end


function val = func1(xi1,xi2, rho_0, Phi_xi, n,k)

[ux, uy] = compute_grad(n,k,Phi_xi,xi1,xi2);


val = (ux^2 + uy^2)*rho_0(xi1,xi2);

end