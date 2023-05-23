function [A, Mass] = MatricesSystem(n,sigma,eps,k, bc_t, alpha_xy)
% alpha_xy : diffusion function
% bc_t     : type of bc: Dirichlet: bc_t=1 | Neumann: bc_t=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('bc_t', 'var')
    bc_t = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nloc = (k+1)^2;                %local dimension matrix
gdim = Nloc*n^2;               %global dimension matrix
A    = zeros(gdim,gdim);
Mass = zeros(gdim,gdim);

%% Assemble volume contrubution
Masselem = localMassMatrix_vol(n,k);
kk = 0;
if ~exist('alpha_xy', 'var')
    Aelem = localMatrix_vol(n,k);
end
for ielm=1:n^2
    if exist('alpha_xy', 'var')
        Aelem = localMatrix_vol(n,k,ielm, alpha_xy);
    end
    for i=1:Nloc
        ie = i + kk;
        for j=1:Nloc
            je = j + kk;
            A(ie,je)     = Aelem(i,j);
            Mass(ie, je) = Masselem(i,j);
        end
    end
    kk = kk + Nloc;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble edges contribution
% (loop over intriot edges. using vertex N1 as intermediate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Horizontale interior edge :
verHor = 0;
if ~exist('alpha_xy', 'var')
    [M11, M22, M12, M21] = localMatrix_edge(verHor,n,k,eps,sigma);
end
for N1=n+2:(n+1)^2 - (n+1)
    if (mod(N1,n+1)~=0)
        [E1, E2] = get_neighbors(N1,verHor,n);
        
        if exist('alpha_xy', 'var')
            [M11, M22, M12, M21] = localMatrix_edge(verHor,n,k,eps,sigma, E1, alpha_xy);
        end
        
        % M11
        for i=1:Nloc
            ie = i + (E1-1)*Nloc;
            for j=1:Nloc
                je = j + (E1-1)*Nloc;
                A(ie,je) = A(ie,je) + M11(i,j);
            end
        end
        % M22
        for i=1:Nloc
            ie = i + (E2-1)*Nloc;
            for j=1:Nloc
                je = j + (E2-1)*Nloc;
                A(ie,je) = A(ie,je) + M22(i,j);
            end
        end
        
        %M12
        for i=1:Nloc
            ie = i + (E1-1)*Nloc;
            for j=1:Nloc
                je = j + (E2-1)*Nloc;
                A(ie,je) = A(ie,je) + M12(i,j);
            end
        end
        
        %M21
        for i=1:Nloc
            ie = i + (E2-1)*Nloc;
            for j=1:Nloc
                je = j + (E1-1)*Nloc;
                A(ie,je) = A(ie,je) + M21(i,j);
            end
        end
    end
end

if (bc_t==1)
    % horizontal boundary edge
    % Bottom:
    ort = 1; verHor = 0;
    if ~exist('alpha_xy', 'var')
        M11 = localMatrix_edgeBoundary(verHor,ort,n,k,eps,sigma);
    end
    for N1=1:n
        [E1, ~] = get_neighbors(N1,verHor,n);
        
        if exist('alpha_xy', 'var')
            M11 = localMatrix_edgeBoundary(verHor,ort,n,k,eps,sigma, E1, alpha_xy);
        end
        
        % M11
        for i=1:Nloc
            ie = i + (E1-1)*Nloc;
            for j=1:Nloc
                je = j + (E1-1)*Nloc;
                A(ie,je) = A(ie,je) + M11(i,j);
            end
        end
    end
    
    % Top:
    ort = 0; verHor = 0;
    if ~exist('alpha_xy', 'var')
        M11 = localMatrix_edgeBoundary(verHor,ort,n,k,eps,sigma);
    end
    for N1=(n+1)^2-n:(n+1)^2-1
        [~, E1] = get_neighbors(N1,verHor,n);
        
        if exist('alpha_xy', 'var')
            M11 = localMatrix_edgeBoundary(verHor,ort,n,k,eps,sigma, E1, alpha_xy);
        end
        
        % M11
        for i=1:Nloc
            ie = i + (E1-1)*Nloc;
            for j=1:Nloc
                je = j + (E1-1)*Nloc;
                A(ie,je) = A(ie,je) + M11(i,j);  % normal_top = -normal_interior
            end
        end
        
    end
end

%% vertical interior edge :
verHor = 1;
if ~exist('alpha_xy', 'var')
    [M11, M22, M12, M21] = localMatrix_edge(verHor,n,k,eps,sigma);
end
for N1=2:(n+1)^2 - (n+1)
    if (mod(N1,n+1)~=0 && mod(N1-1,n+1))
        [E1, E2] = get_neighbors(N1,verHor,n);
        
        if exist('alpha_xy', 'var')
            [M11, M22, M12, M21] = localMatrix_edge(verHor,n,k,eps,sigma, E1, alpha_xy);
        end
        
        % M11
        for i=1:Nloc
            ie = i + (E1-1)*Nloc;
            for j=1:Nloc
                je = j + (E1-1)*Nloc;
                A(ie,je) = A(ie,je) + M11(i,j);
            end
        end
        % M22
        for i=1:Nloc
            ie = i + (E2-1)*Nloc;
            for j=1:Nloc
                je = j + (E2-1)*Nloc;
                A(ie,je) = A(ie,je) + M22(i,j);
            end
        end
        
        %M12
        for i=1:Nloc
            ie = i + (E1-1)*Nloc;
            for j=1:Nloc
                je = j + (E2-1)*Nloc;
                A(ie,je) = A(ie,je) + M12(i,j);
            end
        end
        
        %M21
        for i=1:Nloc
            ie = i + (E2-1)*Nloc;
            for j=1:Nloc
                je = j + (E1-1)*Nloc;
                A(ie,je) = A(ie,je) + M21(i,j);
            end
        end
    end
end

if (bc_t==1)
    % vertical boundary edge
    % left:
    ort = 1; verHor = 1;
    if ~exist('alpha_xy', 'var')
        M11 = localMatrix_edgeBoundary(verHor,ort,n,k,eps,sigma);
    end
    for N1=1:n+1:n^2
        [E1, ~] = get_neighbors(N1,verHor,n);
        
        if exist('alpha_xy', 'var')
            M11 = localMatrix_edgeBoundary(verHor,ort,n,k,eps,sigma, E1, alpha_xy);
        end
        
        % M11
        for i=1:Nloc
            ie = i + (E1-1)*Nloc;
            for j=1:Nloc
                je = j + (E1-1)*Nloc;
                A(ie,je) = A(ie,je) + M11(i,j);
            end
        end
    end
    
    % right:
    ort = 0; verHor = 1;
    if ~exist('alpha_xy', 'var')
        M11 = localMatrix_edgeBoundary(verHor,ort,n,k,eps,sigma);
    end
    for N1=n+1:n+1:(n+1)^2 - (n+1)
        [~, E1] = get_neighbors(N1,verHor,n+1);
        
        if exist('alpha_xy', 'var')
            M11 = localMatrix_edgeBoundary(verHor,ort,n,k,eps,sigma, E1, alpha_xy);
        end
        
        % M11
        for i=1:Nloc
            ie = i + (E1-1)*Nloc;
            for j=1:Nloc
                je = j + (E1-1)*Nloc;
                A(ie,je) = A(ie,je) + M11(i,j);
            end
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% spy(A);
% grid on
%disp(A);

end