clear
clc
close all
all_marks = {'o','+','*','s','^','p'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Verifier l'ordre de convergence sur la solution du probleme de l'obstacle
%%% -> Algorithme avec equation linearisée
%%% -> Laplacien + descriminant hessien au sens fort
%%% -> on résount min(u-g,-Δu-f) = 0 dans Ω
%%% -> u = u_exct sur ∂Ω                    
%%% -> c = 1e3

%%% -> la solution exacte est Δu = - f 
%%% -> La base de Lagrange sur les points de quadrature 
%%% -> NIPG, sigma =10

%%% -> u_exct(x,y) = sin(pi*x) + sin(pi*y)
%%% -> g(x,y) = sin(pi*x).*cos(2*pi*x).^2  + sin(pi*y).*cos(2*pi*y).^2 sur Ω


%%% -> k= 2 
n_k2 = [2 4 8 16 32];
%err2_k2 = [0.02582239 0.006056285 0.00148937 0.0003724169 9.308232e-05];  
%err2_k2_det =[9.45997 2.790365   0.7249961  0.1834325 0.06286849]; 

err2_k2 = [0.02582239  0.006056285 0.00148937 0.000372417 9.309112e-05]; %w.r.t (u) 
err2_k2_det = [9.45997 2.790365 0.7249961 0.1834325 0.04597416]; %w.r.t det(u) 

coeff_k2 = polyfit(log(1./n_k2),log(err2_k2_det),1)



%%% -> k= 3
n_k3 = [2 4 8 16 32];
err2_k3 = [0.006995482 0.0003897485 2.35365e-05 1.457272e-06 9.083458e-08]; % %w.r.t (u) 
err2_k3_det = [6.401128 1.532836 0.3774553 0.09400774 0.02347954]; %w.r.t det(u) 

coeff_k3 = polyfit(log(1./n_k3),log(err2_k3_det),1)



%%% -> k= 4
n_k4 = [2 4 8 16];
err2_k4 = [0.0002795525 1.805764e-05 1.129623e-06 7.052905e-08];  %w.r.t (u) 
err2_k4_det =[0.1592159 0.01306694 0.0008948391 5.725017e-05];    %w.r.t det(u) 

coeff_k4 = polyfit(log(1./n_k4),log(err2_k4_det),1)


%%% -> k= 5
n_k5 = [2 4 8];
err2_k5 = [7.380058e-06 1.283426e-07 2.111433e-09];  %w.r.t (u) 
err2_k5_det =[0.07357469 0.004608114 0.0002876344]; %w.r.t det(u) 

coeff_k5 = polyfit(log(1./n_k5),log(err2_k5_det),1)



figure
loglog(1./n_k2,err2_k2_det,'marker',all_marks{1}, 'LineWidth',1)
hold on
loglog(1./n_k3,err2_k3_det,'marker',all_marks{2}, 'LineWidth',1)
hold on
loglog(1./n_k4,err2_k4_det,'marker',all_marks{3}, 'LineWidth',1)
hold on
loglog(1./n_k5,err2_k5_det,'marker',all_marks{4}, 'LineWidth',1)
grid on


text(0.6*1/n_k2(end),err2_k2_det(end) ,'$\mathcal{O}(h^{1.93})$','interpreter','latex');
text(0.6*1/n_k3(end),err2_k3_det(end) ,'$\mathcal{O}(h^{2.02})$','interpreter','latex');
text(0.6*1/n_k4(end),err2_k4_det(end) ,'$\mathcal{O}(h^{3.82})$','interpreter','latex');
text(1.1*1/n_k5(end),err2_k5_det(end) ,'$\mathcal{O}(h^{4.00})$','interpreter','latex');


legend('$\bf{P}_2$','$\bf{P}_3$','$\bf{P}_4$','$\bf{P}_5$' ,'interpreter','latex')

xlabel('$h$', 'interpreter','latex')
ylabel('$|| \det u - \det  u_{exct} ||_{L^2}$', 'interpreter','latex')

