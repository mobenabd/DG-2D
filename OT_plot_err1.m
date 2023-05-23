clear; clc; close all;

all_marks = {'o','+','*','s','^','p'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Verifier l'ordre de convergence sur \delta rho0_tilde(xi)  
%%%%% lambda = @(x,y) ss * ( x.^3/3 - x.^2/2 + y.^3/3 - y.^2/2 +1/6 ); ss=0.6
%%%%% dlambda = ss0 * ( x.^3/3 - x.^2/2 + y.^3/3 - y.^2/2 +1/6 ); ss0 <<1
%%%%% On fixe k=4, n=16, et on fait varier ss0


n_k4 = [2^4 2^5 2^6 2^7 2^8];


%%% avec \delta U = \nabla_x\delta\lambda  | -> ordre 2 non vérifié !
%err2_k4_2   = [0.219822 0.079486 0.032545 0.014605 0.006910];
%err2_k4_inf = [2.097403 0.779593 0.321345 0.144072 0.067949];


%%% avec $\delta U = \nabla_\xi\delta\Phi$ | -> ordre 2 vérifié !
err2_k4_2   = [0.123898 0.031477 0.007915 0.001989 0.000501];
err2_k4_inf = [1.030675 0.265319 0.067141 0.016932 0.004273];




coeff_k4_2 = polyfit(log(1./n_k4),log(err2_k4_2),1);
coeff_k4_inf = polyfit(log(1./n_k4),log(err2_k4_inf),1);


figure
loglog(1./n_k4,err2_k4_2,'marker',all_marks{1}, 'LineStyle', '-.', 'LineWidth',1)
hold on
loglog(1./n_k4,err2_k4_inf,'marker',all_marks{2},'LineStyle', '-','LineWidth',1)
grid on
text(0.55*1/n_k4(end),err2_k4_2(end) ,['$\mathcal{O}(h^{' num2str(coeff_k4_2(1),'%.2f') '})$'], ...
        'interpreter','latex', 'fontsize', 12);
text(0.55*1/n_k4(end),err2_k4_inf(end) ,['$\mathcal{O}(h^{' num2str(coeff_k4_inf(1),'%.2f') '})$'], ...
        'interpreter','latex', 'fontsize', 12);
xlabel('$\sigma_0$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('error-norm', 'interpreter', 'latex', 'fontsize', 12)
legend('$|| \tilde{\rho}_0^{(1)} - \tilde{\rho}_0^{(0)} - \delta\tilde{\rho}_0 ||_{L^2}$', ...
       '$|| \tilde{\rho}_0^{(1)} - \tilde{\rho}_0^{(0)} - \delta\tilde{\rho}_0 ||_{L^\infty}$', ...
    'interpreter', 'latex', 'fontsize', 12, 'location', 'best')
%title('taking $\delta U(\xi) = \nabla_\xi\delta\Phi(\xi)$', 'interpreter','latex', 'fontsize', 12)


title('taking $\delta U(\xi) = \nabla_x\delta\lambda(x)$', 'interpreter','latex', 'fontsize', 12)

