function Mass = ACM_GlobalMass(n,k,c,U,LAP,g_projected, f_projected)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get global mass matrix M^k for ACM method 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nloc = (k+1)^2;                %local dimension matrix
gdim = Nloc*n^2;               %global dimension matrix
    
Mass = zeros(gdim,gdim);

if (~exist('f_projected', 'var'))
    f_projected = zeros(gdim,1);
end

%% Assemble volume contrubution
kk = 0;
for ielm=1:n^2
    Masselem = ACM_localMassMatrix_vol(n,k,ielm,c,U,LAP,f_projected,g_projected);
    for i=1:Nloc
        ie = i + kk;
        for j=1:Nloc
            je = j + kk;
            Mass(ie, je) = Masselem(i,j);
        end
    end
    kk = kk + Nloc;
end

end
