clear; clc; close all;
%'+' | 'o' | '*' | '.' | 'x' | 'square' | 'diamond' | 'v' | '^' | '>' | '<' |
%'pentagram' | 'hexagram' | 'none'.
all_marks = {'o','+','*','s','^','p'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Verifier l'ordre de convergence sur rho0_tilde(xi) = rho1(X(xi))*det(JacX(xi)) 
%%%%% lambda = @(x,y) ss * ( (x.^3/3 - x.^2/2) + (y.^3/3 - y.^2/2) ); ss=0.6


test = 2;  %% 1: X(xi) exact | 2: X(xi) par Gauss-Newton


%%% -> k= 2
if (test==1)
    %%% X(xi) exact
    n_k2 = [4 8 16 32 64];
    err2_k2_2   = [0.9086966939 0.3650241213 0.1310624008 0.0455164694 0.0347565571];  
    err2_k2_inf = [9.1136942554 3.7015801160 2.0650794158 0.4939637257 0.4894694585];  
elseif (test==2)
    %%% X(xi) by Guass-Newton
    n_k2 = [4 8 16 32];
    err2_k2_2   = [1.0231448526  0.3440726202 0.1317839867 0.0447644679];  
    err2_k2_inf = [10.6790420041 3.2704132536 2.0927822601 0.4888465514];
end 
    
coeff_k2_2 = polyfit(log(1./n_k2),log(err2_k2_2),1);
coeff_k2_inf = polyfit(log(1./n_k2),log(err2_k2_inf),1);


%%% -> k= 3
if (test==1)
    %%% X(xi) exact
    n_k3 = [4 8 16 32 64];
    err2_k3_2   = [0.4214670561 0.0321726589 0.0020414407 0.0002294004 0.0000100396];  
    err2_k3_inf = [3.7790305557 0.4187258619 0.0218184903 0.0034753467  0.0001587751];  
elseif (test==2)
    %%% X(xi) by Guass-Newton
    n_k3 = [4 8 16 32];
    err2_k3_2   = [0.4214675828 0.0321727765 0.0020414127 0.0002298099];  
    err2_k3_inf = [3.7790396581 0.4187245404 0.0218216863 0.0034820098];
end 
    
coeff_k3_2 = polyfit(log(1./n_k3),log(err2_k3_2),1);
coeff_k3_inf = polyfit(log(1./n_k3),log(err2_k3_inf),1);




%%% -> k= 4
if (test==1)
    %%% X(xi) exact
    n_k4 = [4 8 16 32 64];
    err2_k4_2 = [0.1416043893 0.0073614078 0.0002757561  0.0000060839 0.0000001781];  
    err2_k4_inf = [1.4852860999 0.1136622813 0.0034227387 0.0000658518 0.0000024111];
elseif (test==2)
    %%% X(xi) by Guass-Newton
    n_k4 = [4 8 16 32];
    err2_k4_2   = [0.1416044825 0.0073612888 0.0002758355 0.0000060218];
    err2_k4_inf = [1.4852847454 0.1136591991 0.0034215700 0.0000697708];
end 
    
coeff_k4_2 = polyfit(log(1./n_k4),log(err2_k4_2),1);
coeff_k4_inf = polyfit(log(1./n_k4),log(err2_k4_inf),1);




%%% -> k= 5
if (test==1)
    %%% X(xi) exact
    n_k5 = [4 8 16 32];
    err2_k5_2 = [0.0465146618 0.0012508642 0.0000267736 0.0000003479];
    err2_k5_inf = [0.4176913400 0.0130845204 0.0004253213 0.0000044971];
elseif (test==2)
    %%% X(xi) by Guass-Newton
    n_k5 = [4 8 16 32];
    err2_k5_2   = [0.0465145012 0.0012509541 0.0000264697 0.0000009991];
    err2_k5_inf = [0.4176984020 0.0130915657 0.0004177000 0.0000112574];
end 
    
coeff_k5_2 = polyfit(log(1./n_k5),log(err2_k5_2),1);
coeff_k5_inf = polyfit(log(1./n_k5),log(err2_k5_inf),1);



figure(1)
loglog(1./n_k2,err2_k2_2,'marker',all_marks{1}, 'LineWidth',1)
hold on
loglog(1./n_k3,err2_k3_2,'marker',all_marks{2}, 'LineWidth',1)
hold on
loglog(1./n_k4,err2_k4_2,'marker',all_marks{3}, 'LineWidth',1)
hold on
loglog(1./n_k5,err2_k5_2,'marker',all_marks{4}, 'LineWidth',1)

grid on
legend('$\bf{P}_2$','$\bf{P}_3$', '$\bf{P}_4$', '$\bf{P}_5$' ,'interpreter','latex', 'Location','best', 'fontsize', 12)
text(0.6*1/n_k2(end),err2_k2_2(end) ,['$\mathcal{O}(h^{' num2str(coeff_k2_2(1),'%.2f') '})$'],'interpreter','latex', 'fontsize', 12);
text(0.6*1/n_k3(end),err2_k3_2(end) ,['$\mathcal{O}(h^{' num2str(coeff_k3_2(1),'%.2f') '})$'],'interpreter','latex', 'fontsize', 12);
text(0.6*1/n_k4(end),err2_k4_2(end) ,['$\mathcal{O}(h^{' num2str(coeff_k4_2(1),'%.2f') '})$'],'interpreter','latex', 'fontsize', 12);
text(0.6*1/n_k5(end),err2_k5_2(end) ,['$\mathcal{O}(h^{' num2str(coeff_k5_2(1),'%.2f') '})$'],'interpreter','latex', 'fontsize', 12);
xlabel('$h$', 'interpreter','latex', 'fontsize', 12)
ylabel('$|| \rho_0 - \tilde{\rho}_0 ||_{L^2}$', 'interpreter','latex', 'fontsize', 12)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
loglog(1./n_k2,err2_k2_inf,'marker',all_marks{1}, 'LineWidth',1)
hold on
loglog(1./n_k3,err2_k3_inf,'marker',all_marks{2}, 'LineWidth',1)
hold on
loglog(1./n_k4,err2_k4_inf,'marker',all_marks{3}, 'LineWidth',1)
hold on
loglog(1./n_k5,err2_k5_inf,'marker',all_marks{4}, 'LineWidth',1)

grid on
legend('$\bf{P}_2$', '$\bf{P}_3$', '$\bf{P}_4$', '$\bf{P}_5$' ,'interpreter','latex', 'Location','best', 'fontsize', 12)
text(0.6*1/n_k2(end),err2_k2_inf(end) ,['$\mathcal{O}(h^{' num2str(coeff_k2_inf(1),'%.2f') '})$'],'interpreter','latex', 'fontsize', 12);
text(0.6*1/n_k3(end),err2_k3_inf(end) ,['$\mathcal{O}(h^{' num2str(coeff_k3_inf(1),'%.2f') '})$'],'interpreter','latex', 'fontsize', 12);
text(0.6*1/n_k4(end),err2_k4_inf(end) ,['$\mathcal{O}(h^{' num2str(coeff_k4_inf(1),'%.2f') '})$'],'interpreter','latex', 'fontsize', 12);
text(0.6*1/n_k5(end),err2_k5_inf(end) ,['$\mathcal{O}(h^{' num2str(coeff_k5_inf(1),'%.2f') '})$'],'interpreter','latex', 'fontsize', 12);
xlabel('$h$', 'interpreter','latex', 'fontsize', 12)
ylabel('$|| \rho_0 - \tilde{\rho}_0 ||_{L^\infty}$', 'interpreter','latex', 'fontsize', 12)






