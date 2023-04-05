clear; clc; close all;
%%%%%%%%%%%%% Solve min(-Δu, u-g)=0 in Ω = [0,1].^2 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  with ACM method      %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% using DG methods with Lgarange basis on quadrature nodes %%%%
%%%%%%%%%%%%%         u = u_final on ∂Ω              %%%%%%%%%%%%%%%%%%%%%%
% This main is used for tests without exact solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global k
n     =  8;
h     = 1/n;
sigma =  50;
eps   =  -1;
k    =  3;
c    = 1e3;
gdim = n^2*(k+1)^2;


f    =  @(x,y,t) x*0 ;
%f    =  @(x,y,t) source(x,y);
g = @(x,y,t) obstacle(x,y);
u_final = @(x,y,t) final_sol(x,y);

%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%
[A,Mass] = MatricesSystem(n,sigma,eps,k);

Uold = zeros(gdim,1);

    %%%% adding noise %%%%%%%
    %g_projected = computeDirectProjection(n,k,u_final);
    %noiseAmplitude = 0.05;
    %g_projected = g_projected + noiseAmplitude * 2*(rand(gdim,1)-0.5);
    %U0 = plot_sol(n,k, g_projected);
    %%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% without noise %%%%%%%
    %b = SourceBCSystem(n,sigma,eps,k, g,f,0.0);
    %g_projected = Mass\b;
    g_projected = computeDirectProjection(n,k,g);
    %%%%%%%%%%%%%%%%%%

U = g_projected;
plot_sol(n,k, U, g);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lap = computeLap2_ddl(n,k,U);
% Eig_m = computeDet_ddl(n,k,U, 1);  %%compute largest eigenvalue
% f_projected = -Lap + Eig_m ;

%b = SourceBCSystem(n,sigma,eps,k, f,f,0.0);
%f_projected = Mass\b;
%plot_sol(n,k, f_projected);

f_projected = zeros(gdim,1) ;
 


res = 1;
t=1;
while (res > 1e-6 )
    Lap = computeLap2_ddl(n,k,U);   %%compute determinant
    DET = computeDet_ddl(n,k,U);    %%compute determinant
    Eig_m = computeDet_ddl(n,k,U,1);  %%compute largest eigenvalue
    Mk =  ACM_GlobalMass(Eig_m, DET, n,k,c,U,Lap,g_projected, f_projected);
    RHS = ACM_RHS(Eig_m, DET, n,k,c,U,Lap,sigma,eps,g,g_projected, f_projected);
    
    
    Uold = U;
    U = (A+Mk)\RHS;
    %plot_sol(n,k, Lap);
    plot_sol(n,k, U, g);
    %plot_sol3(n,k, U, U0)
    


    %%%%%% For Monge-Ampere %%%%%%%%%
    %RHS_ma = MA_RHS(Eig_m, DET, n,k,c,U,Lap,sigma,eps,g,g_projected, f_projected);
    %U = A\RHS_ma;
    %plot_sol(n,k, U, g);
    %return
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    res = sqrt(sum((U-Uold).^2))/sqrt(sum(U.^2));
    fprintf('res=%f, counter=%i\n', res, t);
    t = t + 1;
    
end
fprintf('\n')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = obstacle(tx,ty)
out = sin(2*pi*tx)*sin(2*pi*ty);


% x = 2*tx-1;
% y = 2*ty-1;
% out = -max([0.5*x, y-0.5, 2*x+y-1, -5*x+y-4, 2-10*x.^2 - 10*y.^2]);

%out = - (sin(pi*(x-y))./exp(2*(x.^2+y.^2)) + x.^2 + y.^2);




%out = -max([0.5*x, y-0.5, 2*x+y-1, -5*x+y-4]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x = 2*tx-1;
% y = 2*ty-1;
% 
% out = -(0.98 - x.^2 - y.^2)^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x = 4*tx-2;
% y = 4*ty-2;
% 
% rxy = x.^2 + y.^2;
% if (rxy <= 1 )
%     out = sqrt(1 - rxy);
% else
%     out = -1;
% end
% 
% if ((x==-2 ||  x==2 || y==-2 || y==2))
%     out = final_sol(tx,ty);
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%% For Monge-Ampere solution %%%%%%%%%%%%%%
%out = abs(tx-0.5);
% x = 2*tx-1;
% y = 2*ty-1;
% out = exp((x.^2 + y.^2)/2);

%out = 2*sqrt(2)/3 * (tx.^2 + ty.^2).^(3/4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end




function out = final_sol(tx,ty)
%out = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = 4*tx-2;
y = 4*ty-2;

r = 0.69796514822337356642145262626753;

rxy = sqrt(x.^2 + y.^2);

if (rxy <= r )
    out = sqrt(1 - rxy.^2);
else
    out = -r^2*log(rxy/2)/sqrt(1-r^2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end



function out = source(tx,ty)

%%%%%%%%%%% For Monge-Ampere solution %%%%%%%%%%%%%%
% x = 2*tx-1;
% y = 2*ty-1;
% out = 16*(1 + x.^2 + y.^2).*exp(x.^2 + y.^2);

%out = 1/sqrt(tx.^2 + ty.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out = 0;
end



