clear
close all
clc
% 
% 
% n = [4 8 16 32 64];
% err2 = [0.001058  0.000544 0.000278 0.000131 0.000061];
% err2 = [0.021155  0.010912 0.005576 0.002639 0.001251];
% err2 = [0.055140  0.025836 0.012682 0.006290 0.003110];  %(avec dirichlet non homegene codé k=1)
% 
% 
% %%% sinusoidale k=1 sigma = 100 %%%%%%%%%
% n = [4 8 16 32 64 70];
% err2 = [1.548555 0.936879 0.499916 0.255493 0.128826 0.120421];
% 
% %%% polynome k=2 sigma = 100 %%%%
% n = [49 32 16 8 4];
% err2 = [0.004133 0.006301 0.012721 0.025952 0.055766];
% 
% %%% sinusoidale k=2 sigma = 100 %%%%%
% n = [49 32 16 8 4];
% err2 = [0.162624 0.255366 0.499999 0.943422 1.637507];
% 
% 
% %%% polynome k=2 sigma = 100 Dirichlet non homegene non codé %%
% n = [4 8 16 32 49];
% err2 = [0.015103 0.008007 0.004125 0.002126 0.001509];



%%% exponentielle (dirichlet=0) k=2 sigma=100
% n= [4 8 16 32];
% err2 = [0.000463 0.000059 0.000007 0.000001];




%%% sinusoidale (dirichlet=0) k=1 sigma =100
n_k1 = [4 8 16 32 64];
err2_k1 =[0.119994 0.030294 0.007595 0.001900 0.000475];



%%% sinusoidale (direchlet=0) k=2 sigma=100
n_k2 = [4 8 16 32 49];
err2_k2 = [0.013948 0.001858 0.000236 0.000030 0.000006];


%%% sinusoidale (dirichlet=0) k=3 sigma =100
n_k3 = [2 4 6 8 12];
err2_k3= [0.025827 0.001322  0.000272  0.000088 0.000017];

coeff_k1 = polyfit(log(1./n_k1),log(err2_k1),1);
coeff_k2 = polyfit(log(1./n_k2),log(err2_k2),1);
coeff_k3 = polyfit(log(1./n_k3),log(err2_k3),1);

figure
loglog(1./n_k1,err2_k1, 'marker','s');
hold on
loglog(1./n_k2,err2_k2, 'marker','s');
hold on
loglog(1./n_k3,err2_k3, 'marker','s');
%loglog(1./n,1./n.^2)
grid on;
legend('$\bf{Q}_1$','$\bf{Q}_2$','$\bf{Q}_3$','interpreter','latex')
xlabel('$h$', 'interpreter','latex')
ylabel('$|| e ||_{L^2}$', 'interpreter','latex')
text(0.65*1/n_k1(end),err2_k1(end) ,'$\mathcal{O}(h^2)$','interpreter','latex');
text(0.65*1/n_k2(end),err2_k2(end) ,'$\mathcal{O}(h^3)$','interpreter','latex');
text(0.65*1/n_k3(end),err2_k3(end) ,'$\mathcal{O}(h^4)$','interpreter','latex');


% legend(['y= ',num2str(coeff_k1(1)),'x + ',num2str(coeff_k1(2))], ...
%        ['y= ',num2str(coeff_k2(1)),'x + ',num2str(coeff_k2(2))],...
%        ['y= ',num2str(coeff_k3(1)),'x + ',num2str(coeff_k3(2))],...
% 'interpreter','latex')
