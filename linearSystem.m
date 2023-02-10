function [A, Mass ,b] = linearSystem(n,sigma,eps,k,uexct,f)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nloc = (k+1)^2;                %local dimension matrix
gdim = Nloc*n^2;               %global dimension matrix
A    = zeros(gdim,gdim);      
b    = zeros(gdim,1);

Mass = zeros(gdim,gdim);

%% Assemble volume contrubution
Masselem = localMassMatrix_vol(n,k);
Aelem = localMatrix_vol(n,k);
kk = 0;
for ielm=1:n^2
    belem = localRHS_vol(f,n,k,ielm);
    for i=1:Nloc
        ie = i + kk;
        for j=1:Nloc
            je = j + kk;
            A(ie,je) = Aelem(i,j);
            Mass(ie, je) = Masselem(i,j);
        end
        b(ie) = belem(i);
    end
    kk = kk + Nloc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble edges contribution
% (loop over intriot edges. using vertex N1 as intermediate) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Horizontale interior edge :
verHor = 0;
[M11, M22, M12, M21] = localMatrix_edge(verHor,n,k,eps,sigma);
for N1=n+2:(n+1)^2 - (n+1)
    if (mod(N1,n+1)~=0)
        [E1, E2] = get_neighbors(N1,verHor,n);
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

% horizontal boundary edge
% Bottom:
ort = 1;
M11 = localMatrix_edgeBoundary(verHor,ort,n,k,eps,sigma);
for N1=1:n
    [E1, ~] = get_neighbors(N1,verHor,n);
    belem = localRHS_edgeBoundary(verHor,ort,E1,n,k,eps,sigma,uexct);
    % M11
    for i=1:Nloc
        ie = i + (E1-1)*Nloc;
        b(ie) = b(ie) + belem(i);
        for j=1:Nloc
            je = j + (E1-1)*Nloc;
            A(ie,je) = A(ie,je) + M11(i,j);
        end
    end
end


% Top:
ort = 0;
M11 = localMatrix_edgeBoundary(verHor,ort,n,k,eps,sigma);
for N1=(n+1)^2-n:(n+1)^2-1
    [~, E1] = get_neighbors(N1,verHor,n);
    belem = localRHS_edgeBoundary(verHor,ort,E1,n,k,eps,sigma,uexct);
    % M11
    for i=1:Nloc
        ie = i + (E1-1)*Nloc;
        b(ie) = b(ie) + belem(i);
        for j=1:Nloc
            je = j + (E1-1)*Nloc;
            A(ie,je) = A(ie,je) + M11(i,j);  % normal_top = -normal_interior 
        end
    end
    
end


%% vertical interior edge :
verHor = 1;
[M11, M22, M12, M21] = localMatrix_edge(verHor,n,k,eps,sigma);
for N1=2:(n+1)^2 - (n+1)
    if (mod(N1,n+1)~=0 && mod(N1-1,n+1))
        [E1, E2] = get_neighbors(N1,verHor,n);
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

% vertical boundary edge
% left:
ort = 1;
M11 = localMatrix_edgeBoundary(verHor,ort,n,k,eps,sigma);
for N1=1:n+1:n^2
    [E1, ~] = get_neighbors(N1,verHor,n);
    belem = localRHS_edgeBoundary(verHor,ort,E1,n,k,eps,sigma,uexct);
    % M11
    for i=1:Nloc
        ie = i + (E1-1)*Nloc;
        b(ie) = b(ie) + belem(i);
        for j=1:Nloc
            je = j + (E1-1)*Nloc;
            A(ie,je) = A(ie,je) + M11(i,j);
        end
    end
end

% right:
ort = 0;
M11 = localMatrix_edgeBoundary(verHor,ort,n,k,eps,sigma);
for N1=n+1:n+1:(n+1)^2 - (n+1)
    [~, E1] = get_neighbors(N1,verHor,n+1);
    belem = localRHS_edgeBoundary(verHor,ort,E1,n,k,eps,sigma,uexct);
    % M11
    for i=1:Nloc
        ie = i + (E1-1)*Nloc;
        b(ie) = b(ie) + belem(i);
        for j=1:Nloc
            je = j + (E1-1)*Nloc;
            A(ie,je) = A(ie,je) + M11(i,j);
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% spy(A);
% grid on
%disp(A);
%disp(b);

%%Solve
%U = A\b;

end