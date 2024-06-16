function compute_error(U,u_exct,n,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute L2 error with respect to the exact solution %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x0, xN, ~] = getGlobal_x0N();
dx= (xN-x0)*0.01;

x = x0:dx:xN;
y = x;


M = zeros(length(x),length(y));
Mexct = zeros(length(x),length(y));
for i = 1:length(x)
    for j=1:length(y)
        M(i,j) = compute_sol(x(i),y(j),U,n,k);
        Mexct(i,j) = u_exct(x(i),y(j));
    end
end

err2 = dx*sqrt(sum(sum((Mexct-M).^2)));
fprintf('|e|_2 = %.10f\n',err2);


errInf = max(abs(Mexct-M), [], 'all');

fprintf('|e|_inf = %.10f\n', errInf);

end