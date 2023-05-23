clear; clc; close all;
%%%%%%%%%%%%% check \Delta \Tilde{\rho_0} expression %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global k
n     =  16;
sigma =  10; eps   =  1;
k    =  4; gdim = n^2*(k+1)^2;



%%%%%%% densities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TC0 %%%%%%%%%%%%
nu0 = 0.4; nu1 = 0.42;
rho_0 = @(x,y)  1/0.03141592 * exp(-100*((x-nu0)^2 + (y-nu0)^2));
rho_1 = @(x,y)  1/0.03141592 * exp(-100*((x-nu1)^2 + (y-nu1)^2)); 

rho0_proj = computeDirectProjection(n,k,rho_0);
rho1_proj = computeDirectProjection(n,k,rho_1);


%%%%%% initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta     = @(x,y) 0.5*(x.^2 + y.^2) ;
beta_x   = computeDirectProjection(n,k,beta);

gamma   = @(x,y)  x.^3/3 - x.^2/2 + y.^3/3 - y.^2/2 ; 
gamma_x = computeDirectProjection(n,k,gamma);

lambda_x = 0.6*(gamma_x + 1/6);



%%% compute dlambda_x
global ss0
ss0 = 1/2^8; %ss0 = 0.005; 
ss0old = ss0;
dlambda_x =  ss0*(gamma_x + 1/6);  %% !! coeff à changer dans OT_Xxi2() !!


lambda_x_eps = lambda_x + dlambda_x;



Alpha_x = beta_x - lambda_x;
Alpha_x_eps = beta_x - lambda_x_eps;



%%% compute rho0_tilde_xi_1
det_x_eps = computeDet_ddl(n,k, Alpha_x_eps); 
rho0_tilde_x = rho1_proj./det_x_eps;   

[rho0_tilde_xi_1, ~] = OT_compute_rho0Tilde_xi_TEST(n,k, rho0_tilde_x);   %%ici ss0 non nulle
rho0_tilde_xi_1(rho0_tilde_xi_1<0) = 0; % positive density
vol =  integrate(rho0_tilde_xi_1,n); rho0_tilde_xi_1 = 1/vol * rho0_tilde_xi_1;




%%% compute rho0_tilde_xi_0
det_x = computeDet_ddl(n,k, Alpha_x); 
rho0_tilde_x = rho1_proj./det_x; 
ss0 = 0;
[rho0_tilde_xi_0, ~] = OT_compute_rho0Tilde_xi_TEST(n,k, rho0_tilde_x);    

rho0_tilde_xi_0(rho0_tilde_xi_0<0) = 0; % positive density
vol =  integrate(rho0_tilde_xi_0,n); rho0_tilde_xi_0 = 1/vol * rho0_tilde_xi_0;


%%% calcul de \delta\Tilde{\rho_0}
d_rho0Tilde_x = OT_compute_d_rho0Tilde_x(n,k, rho1_proj, dlambda_x)./det_x;
    %%ici ss0=0
[d_rho0Tilde_xi, ~] = OT_compute_rho0Tilde_xi_TEST(n,k, d_rho0Tilde_x); 


%%% calcul de \delta\Tilde{\rho_0}(xi) avec \delta U = \nabla_xi\delta\Phi
ss0 = ss0old;
d_rho0Tilde_xi_2 = OT_compute_d_rho0Tilde_xi(n,k, rho1_proj);

%%% ploting %%%%%%%

U1 = rho0_tilde_xi_1-rho0_tilde_xi_0;
U2_1 = d_rho0Tilde_xi;
U2_2 = d_rho0Tilde_xi_2;



% OT_plot_sol(n,k, abs(U1-U2_1));
% title('avec $\delta U = \nabla_x\delta\lambda$', 'interpreter', 'latex' )
% OT_plot_sol(n,k, abs(U1-U2_2),3);
% title('avec $\delta U = \nabla_\xi\delta\Phi$', 'interpreter', 'latex')
% 

% 
% OT_plot_sol(n,k, U1);
% hold on
% OT_plot_sol(n,k, d_rho0Tilde_x);




normL2 =  sqrt(integrate((U1-U2_1).^2,n));
normLinf = max(abs(U1-U2_1));
fprintf("L2-norm=%f | Linf-norm=%f\n", normL2, normLinf);


normL2 =  sqrt(integrate((U1-U2_2).^2,n));
normLinf = max(abs(U1-U2_2));
fprintf("L2-norm=%f | Linf-norm=%f\n", normL2, normLinf);

return




%%%%%%%%%%%%%%% test det[I-\nabla^2\lambda]=det[I+\nabla^2\Phi] %%%%%%%%%%%%
% ss= 0.6;
% func_Phi_xi = @(xi1, xi2) 1/(2*ss)*(  (1+ss)*(xi1+xi2) - 1/(6*ss) *( ((1+ss)^2 - 4*ss*xi1 )^1.5 ...
%                                                                     +((1+ss)^2 - 4*ss*xi2 )^1.5) ) ...
%                          ;% -0.5*(xi1^2 + xi2^2);
% 
% Alpha_xi =  computeDirectProjection(n,k,func_Phi_xi);
% 
% det_xi_Phi  =  computeDet_ddl(n,k, Alpha_xi);
% 
% ss0 = 0;
% [det_xi_lambda, ~] = OT_compute_rho0Tilde_xi_TEST(n,k, det_x); 



%%%%%%%%%%% test if \nabla_x \delta\lambda = \nabla_xi \delta\Phi %%%%%%%%%
%%% not really !!
%GRAD_dlambda_x = computeGrad_ddl(n,k,dlambda_x); 

%%%added
GRAD_dlambda_x = -computeGrad_ddl(n,k,lambda_x);  
%%%


ss0 = ss0old;
ss= 0.6;
% func_dPhi_xi1 = @(xi, dump) OT_Xxi3(xi, ss+ss0) - OT_Xxi3(xi, ss) ;
% func_dPhi_xi2 = @(dump, xi) OT_Xxi3(xi, ss+ss0) - OT_Xxi3(xi, ss) ; 

%%%added
func_dPhi_xi1 = @(xi, dump)  -OT_Xxi3(xi, ss) +xi;
func_dPhi_xi2 = @(dump, xi)  -OT_Xxi3(xi, ss) +xi;
%%%




GRAD_dPhi_xi1 =  computeDirectProjection(n,k, func_dPhi_xi1);  
GRAD_dPhi_xi2 =  computeDirectProjection(n,k, func_dPhi_xi2);  


ss0 = 0;
[GRAD_dlambda_xi1, ~] = OT_compute_rho0Tilde_xi_TEST(n,k, GRAD_dlambda_x(:,1)); 
[GRAD_dlambda_xi2, ~] = OT_compute_rho0Tilde_xi_TEST(n,k, GRAD_dlambda_x(:,2)); 


OT_plot_sol(n,k, GRAD_dlambda_xi1);
hold on
OT_plot_sol(n,k, GRAD_dPhi_xi1);



