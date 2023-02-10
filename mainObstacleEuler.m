clear; clc; close all;
%%%%%%%%%%%%%%%  Solve u_t-Δu = f in Ω = [0,1].^2 %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%         u = u_0  in Ω at t=0      %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%         u = g on ∂Ω               %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n     =  15;
h     = 1/n;
sigma =  100;
eps   =  -1;
k     =  1;
dt    = h^2/4;

f    =  @(x,y,t) x.*0;
g = @(x,y,t) obstacle(x,y);
u_final = @(x,y,t) final_sol(x,y);

%%% Solve with Backward Euler %%%%%%%%%%%%%%%%%%%%%
[A,Mass] = MatricesSystem(n,sigma,eps,k);
%initial contition
b = SourceBCSystem(n,sigma,eps,k,g,f,0);
U = Mass\b;
plot_sol(U, g,n,k,0,2);
%plot_sol(U, u_final,n,k,0,1);
%return

Uold = zeros(length(U),1);
res = 1;
t=1;
%fprintf('Counter: ')
while (res > 1e-8)
    Uold = U;
    b = SourceBCSystem(n,sigma,eps,k,g,f,t*dt);
    U = (Mass + dt*A)\(Mass*U+dt*b);
    
    U= update_obstacle(n,k, U, g);
    
    plot_sol(U, g,n,k,t*dt,2);
    %plot_sol(U, u_final,n,k,t*dt,1);
    res = h*sqrt(sum((U-Uold).^2));
    t = t + 1;
    fprintf('res=%f\n', res);
    
    %     for j=0:log10(t-1)
    %         fprintf('\b'); % delete previous counter display
    %     end
    %     fprintf('%i',t);
    %pause(0.1);
end
fprintf('\n')



function out = obstacle(tx,ty)
x = 4*tx-2;
y = 4*ty-2;

rxy = x.^2 + y.^2;
if (rxy < 1 )
    out = sqrt(1 - rxy);
else
    out = 0;
end


if ((tx==1 ||  tx==0 || ty==1 || ty==0))
    r = 0.69796514822337356642145262626753;
    rxy = sqrt(x.^2 + y.^2);
    if (rxy < r )
        out = sqrt(1 - rxy.^2);
    else
        out = -r^2*log(rxy/2)/sqrt(1-r^2);
    end
end
end


function out = final_sol(tx,ty)
x = 4*tx-2;
y = 4*ty-2;

% syms q;
% eqn = q.^2 .*(1-log(q/2)) == 1;
% r = double(solve(eqn,q));

r = 0.69796514822337356642145262626753;

rxy = sqrt(x.^2 + y.^2);

if (rxy < r )
    out = sqrt(1 - rxy.^2);
else
    out = -r^2*log(rxy/2)/sqrt(1-r^2);
end


end
