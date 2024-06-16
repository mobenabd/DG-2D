function b = SourceBCSystem(n,sigma,eps,k,uexct,f,t, bc_t, alpha_xy)
% alpha_xy : diffusion function
% bc_t     : type of bc: Dirichlet: bc_t=1 | Neumann: bc_t=0 (=1 by default)
if ~exist('bc_t', 'var')
    bc_t = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nloc = (k+1)^2;                %local dimension matrix
gdim = Nloc*n^2;               %global dimension matrix
b    = zeros(gdim,1);

f_source = @(x,y) f(x,y,t); %source term function
f_bc = @(x,y) uexct(x,y,t);  %boundary condition function 



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

if (bc_t==1)
    %% horizontal boundary edge
    verHor = 0;
    % Bottom:
    ort = 1;
    for N1=1:n
        [E1, ~] = get_neighbors(N1,verHor,n);
        if exist('alpha_xy', 'var')
            belem = localRHS_edgeBoundary(verHor,ort,E1,n,k,eps,sigma,f_bc, alpha_xy);
        else
            belem = localRHS_edgeBoundary(verHor,ort,E1,n,k,eps,sigma,f_bc);
        end
        for i=1:Nloc
            ie = i + (E1-1)*Nloc;
            b(ie) = b(ie) + belem(i);
        end
    end

       
    % Top:
    ort = 0; verHor = 0;
    for N1=(n+1)^2-n:(n+1)^2-1
        [~, E1] = get_neighbors(N1,verHor,n);
        if exist('alpha_xy', 'var')
            belem = localRHS_edgeBoundary(verHor,ort,E1,n,k,eps,sigma,f_bc, alpha_xy);
        else
            belem = localRHS_edgeBoundary(verHor,ort,E1,n,k,eps,sigma,f_bc);
        end
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
        if exist('alpha_xy', 'var')
            belem = localRHS_edgeBoundary(verHor,ort,E1,n,k,eps,sigma,f_bc, alpha_xy);
        else
            belem = localRHS_edgeBoundary(verHor,ort,E1,n,k,eps,sigma,f_bc);
        end
        for i=1:Nloc
            ie = i + (E1-1)*Nloc;
            b(ie) = b(ie) + belem(i);
        end
    end
    
    % right:
    ort = 0; verHor = 1;
    for N1=n+1:n+1:(n+1)^2 - (n+1)
        [~, E1] = get_neighbors(N1,verHor,n+1);
        if exist('alpha_xy', 'var')
            belem = localRHS_edgeBoundary(verHor,ort,E1,n,k,eps,sigma,f_bc, alpha_xy);
        else
            belem = localRHS_edgeBoundary(verHor,ort,E1,n,k,eps,sigma,f_bc);
        end
        for i=1:Nloc
            ie = i + (E1-1)*Nloc;
            b(ie) = b(ie) + belem(i);
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (bc_t==0)
    %% horizontal boundary edge
    verHor = 0;
    % Bottom:
    ort = 1;
    for N1=1:n
        [E1, ~] = get_neighbors(N1,verHor,n);
        belem = localRHS_edgeBoundaryNeumann(verHor,ort,E1,n,k,f_bc);
        for i=1:Nloc
            ie = i + (E1-1)*Nloc;
            b(ie) = b(ie) + belem(i);
        end
    end
    
    % Top:
    ort = 0;
    for N1=(n+1)^2-n:(n+1)^2-1
        [~, E1] = get_neighbors(N1,verHor,n);
        belem = localRHS_edgeBoundaryNeumann(verHor,ort,E1,n,k,f_bc);
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
        belem = localRHS_edgeBoundaryNeumann(verHor,ort,E1,n,k,f_bc);
        for i=1:Nloc
            ie = i + (E1-1)*Nloc;
            b(ie) = b(ie) + belem(i);
        end
    end
    
    % right:
    ort = 0;
    for N1=n+1:n+1:(n+1)^2 - (n+1)
        [~, E1] = get_neighbors(N1,verHor,n+1);
        belem = localRHS_edgeBoundaryNeumann(verHor,ort,E1,n,k,f_bc);
        for i=1:Nloc
            ie = i + (E1-1)*Nloc;
            b(ie) = b(ie) + belem(i);
        end
    end
end



end