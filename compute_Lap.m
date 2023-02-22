function Lap = compute_Lap(A,Mass,U,b,f,n,k)
%%%%% Compute weak Laplacian %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nloc = (k+1)^2;                %local dimension matrix
gdim = Nloc*n^2;               %global dimension matrix 
f_source = @(x,y) f(x,y,0.);    %source term function

%%% volume contribution of source term
F = zeros(gdim,1);

kk = 0;
for ielm=1:n^2
    belem = localRHS_t0_vol(f_source,n,k,ielm);
    for i=1:Nloc
        ie = i + kk;
        F(ie) = belem(i);
    end
    kk = kk + Nloc;
end

RHS  = -F - A*U + b;
Lap = Mass\RHS;

%Lap = pcg(Mass,RHS, 1e-12, 1000);


end