function val = integrate(RHO,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Integrate a function defined on quadrature points over the domain [0,1]^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = getGlobal_k();
Nloc = (k+1)^2;                %local dimension matrix
[~,~, detJac] = jacobian_elem(1,n);    % (uniforme grid)

kk = 0;
val = 0;
for ielm=1:n^2
    for i=1:Nloc
        ie = i + kk;
        [phi, ~, ~] = basis(i,n,k);
        val =val + RHO(ie) * detJac * quadrature(phi, 2);
    end

    
    kk = kk + Nloc;
end

end