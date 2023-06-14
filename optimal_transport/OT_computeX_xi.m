function [x,y] = OT_computeX_xi(xi_1,xi_2, M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute X(xi) with bilinear interpolation %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

npt = 2; %%nombre de points d'interpolation

norm = sqrt((xi_1 - M(:,3)).^2 + (xi_2 - M(:,4)).^2);

bx = zeros(4,1);
by = zeros(4,1);

XI = zeros(2,2);
for k=1:npt
    [~, i] = min(norm);
    
    XI(k,1) = M(i,3); %xi_1
    XI(k,2) = M(i,4); %xi_2
    
    bx(k) = M(i,1); 
    by(k) = M(i,2); 
    norm(i) = inf;
end


% x = A(1,2);
% y = A(1,3);

XI1 = XI(:,1);
if (all(diff(sort(XI1(XI1 ~= 0)))) && xi_1> min(XI1) &&  xi_1 < max(XI1))
    x = 0;
    for i=1:npt
         val = 1;
         for k=1:npt
             if (k ~= i)
                 val = val*(XI1(k) - xi_1)/(XI1(k) - XI1(i));
             end
         end
         val = val* bx(i);
         x = x + val;
    end
    %x = bx(1)*(XI(2,1) - xi_1)/(XI(2,1) - XI(1,1))  + bx(2)*(XI(1,1) - xi_1)/(XI(1,1) - XI(2,1));
else
   % x = 0.5*(bx(1) + bx(2));
    x = bx(1);
end

XI2 = XI(:,2);
if (all(diff(sort(XI2(XI2 ~= 0)))) && xi_2> min(XI2) &&  xi_2 < max(XI2))
    y = 0;
    for i=1:npt
         val = 1;
         for k=1:npt
             if (k ~= i)
                 val = val*(XI2(k) - xi_2)/(XI2(k) - XI2(i));
             end
         end
         val = val* by(i);
         y = y + val;
    end
    %y = by(1)*(XI(2,2) - xi_2)/(XI(2,2) - XI(1,2))  + by(2)*(XI(1,2) - xi_2)/(XI(1,2) - XI(2,2));
else
    %y =   0.5*(by(1) + by(2));
    y =   by(1);
end