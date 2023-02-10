function b = SourceBCSystem(n,sigma,eps,k,uexct,f,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nloc = (k+1)^2;                %local dimension matrix
gdim = Nloc*n^2;               %global dimension matrix      
b    = zeros(gdim,1);

f_source = @(x,y) f(x,y,t); %source term function
f_bc = @(x,y) uexct(x,y,t);  %boundary condition function = exact solution





%% Assemble volume contrubution
if (t==0.0)
    kk = 0;
    for ielm=1:n^2
        belem = localRHS_t0_vol(f_bc,n,k,ielm);
        for i=1:Nloc
            ie = i + kk;
            b(ie) = belem(i);
        end
        kk = kk + Nloc;
    end
    return;
end

kk = 0;
for ielm=1:n^2
    belem = localRHS_vol(f_source,n,k,ielm);
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