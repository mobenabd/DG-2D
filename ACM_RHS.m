function b = ACM_RHS(n,k,c,U,LAP,sigma,eps,uexct,g_projected,f_projected)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nloc = (k+1)^2;                %local dimension matrix
gdim = Nloc*n^2;               %global dimension matrix      
b    = zeros(gdim,1);


f_bc = @(x,y) uexct(x,y,0.);  %boundary condition function = exact solution at t=0 (stationary)

if (~exist('f_projected', 'var'))
    f_projected = zeros(gdim,1);
end

%% Assemble volume contrubution
kk = 0;
for ielm=1:n^2
    belem = ACM_localRHS_vol(n,k,ielm,c,U,LAP,f_projected,g_projected);
    for i=1:Nloc
        ie = i + kk;
        b(ie) = belem(i);
    end
    kk = kk + Nloc;
end


%% horizontal boundary edge
verHor = 0;
% Bottom:
ort = 1;
for N1=1:n
    [E1, ~] = get_neighbors(N1,verHor,n);
    belem = localRHS_edgeBoundary(verHor,ort,E1,n,k,eps,sigma,f_bc);
    for i=1:Nloc
        ie = i + (E1-1)*Nloc;
        b(ie) = b(ie) + belem(i);
    end
end

% Top:
ort = 0;
for N1=(n+1)^2-n:(n+1)^2-1
    [~, E1] = get_neighbors(N1,verHor,n);
    belem = localRHS_edgeBoundary(verHor,ort,E1,n,k,eps,sigma,f_bc);
    for i=1:Nloc
        ie = i + (E1-1)*Nloc;
        b(ie) = b(ie) + belem(i);
    end
    
end


%% vertical boundary edge
verHor = 1;
% left:
ort = 1;
for N1=1:n+1:n^2
    [E1, ~] = get_neighbors(N1,verHor,n);
    belem = localRHS_edgeBoundary(verHor,ort,E1,n,k,eps,sigma,f_bc);
    for i=1:Nloc
        ie = i + (E1-1)*Nloc;
        b(ie) = b(ie) + belem(i);
    end
end

% right:
ort = 0;
for N1=n+1:n+1:(n+1)^2 - (n+1)
    [~, E1] = get_neighbors(N1,verHor,n+1);
    belem = localRHS_edgeBoundary(verHor,ort,E1,n,k,eps,sigma,f_bc);
    for i=1:Nloc
        ie = i + (E1-1)*Nloc;
        b(ie) = b(ie) + belem(i);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
%disp(b);
end