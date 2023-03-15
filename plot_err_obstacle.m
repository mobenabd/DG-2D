clear
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Verifier l'ordre de convergence sur la solution du probleme de l'obstacle
%%%%% et le descriminant hessien au sens faible %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% on résount min(u-g,-Δu-f) = 0 dans Ω
%%%%% u = u_exct sur ∂Ω %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% puis on calcule le laplacien projeté 

%%%%% la solution exacte est Δu = - f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% La base de Lagrange sur les points de quadrature %%%%%%%%%%%%%%%%%%%

%%% u_exct(x,y) = sin(pi*x) + sin(pi*y)
%%% g(x,y) = sin(pi*x).*cos(2*pi*x).^2  + sin(pi*y).*cos(2*pi*y).^2 sur Ω
%%% sigma =100


%%%  k=1, SIPG, 100
n_k1 = [4 8 16 32];
%err2_k1 =[0.0399956327 0.0100243402 0.0025068410 0.0006267601]; %w.r.t u
%err2_k1 =[0.2288769558 0.0574018811 0.0143611956 0.0035909560]; %w.r.t lap(u)
%err2_k1 =[8.0969826536 7.8803791463 7.8631494415 7.8607023945]; %w.r.t det(u)

%%% k=1, NIPG,100
%err2_k1 =[0.0387239619 0.0097145667 0.0024341458 0.0006094034]; %w.r.t u
%err2_k1 =[0.2288769558 0.0574018811 0.0143611956 0.0035909560]; %w.r.t lap(u)
err2_k1 =[8.0988871304 7.8804903208 7.8631515709 7.8607019374]; %w.r.t det(u)

coeff_k1 = polyfit(log(1./n_k1),log(err2_k1),1)



%%%  k=2, SIPG, 100
n_k2 = [4 8 16 32];
%err2_k2 =[0.0025703029 0.0003306729  0.0000417993 0.0000052515]; %w.r.t u
%err2_k2 =[0.0162505177 0.0020532983 0.0002572653 0.0000321674]; %w.r.t lap(u)
%err2_k2 =[5.1325753939 2.8295581501 1.4855162870 0.7326677304]; %w.r.t det(u)

%%%  k=2, NIPG, 100
% err2_k2 = [0.0028395970 0.0003983824 0.0000628460 0.0000120590]; %w.r.t u
% err2_k2 = [0.0162505177 0.0020532983 0.0002572653 0.0000321674]; %w.r.t lap(u)
 err2_k2 = [5.1059719268 2.8213500448 1.4835230434 0.7303423986]; %w.r.t det(u)


coeff_k2 = polyfit(log(1./n_k2),log(err2_k2),1)



%%%  k=3, SIPG, 100
n_k3 = [2 4 8 16];
%err2_k3 =[0.0018023359 0.0001210731 0.0000077642 0.0000004902];  %w.r.t u
%err2_k3 =[0.0122457352 0.0007861769 0.0000490176 0.0000030613];  %w.r.t lap(u)
%err2_k3 =[2.0792112069 0.4813267073 0.1261322491 0.0335723810];  %w.r.t det(u)


%%%  k=3, NIPG, 100
% err2_k3 =[0.0019816046 0.0001253819 0.0000079038 0.0000004964];  %w.r.t u
% err2_k3 =[0.0122457352 0.0007861769 0.0000490176 0.0000030613];  %w.r.t lap(u)
 err2_k3 =[2.0774705297 0.4816060846 0.1263318557 0.0335850310];  %w.r.t det(u)

%%% k=3 SIPG, 100, determiant coputed directyly without weak laplacian
%err2_k3 =[7.6827203538 1.8952380288 0.4725148057 0.1180654761];  %w.r.t det(u)


coeff_k3 = polyfit(log(1./n_k3),log(err2_k3),1)



%%%  k=4, SIPG, 100
n_k4 = [2 4 8 16];
%err2_k4 =[0.0001356318 0.0000043606 0.0000001391 0.0000000044];  %w.r.t u
% err2_k4 =[0.0010086327 0.0000344064 0.0000010866 0.0000000340];  %w.r.t lap(u)
% err2_k4 =[0.5692793811 0.0676240809 0.0082802347 0.0010876482];  %w.r.t det(u)


%%%  k=4, NIPG, 100
% err2_k4 = [0.0001692295 0.0000063438 0.0000002618 0.0000000128]; %w.r.t u
% err2_k4 = [0.0010086327 0.0000344064 0.0000010866 0.0000000340]; %w.r.t lap(u)
 err2_k4 = [0.5178826550 0.0619560102 0.0076940711 0.0010202861]; %w.r.t det(u)


coeff_k4 = polyfit(log(1./n_k4),log(err2_k4),1)



figure
% loglog(1./n_k1,err2_k1,'marker','s')
% hold on
loglog(1./n_k2,err2_k2,'marker','s')
hold on
loglog(1./n_k3,err2_k3,'marker','s')
hold on
loglog(1./n_k4,err2_k4,'marker','s')
grid on

%text(0.6*1/n_k1(end),err2_k1(end) ,'$\mathcal{O}(h^{1.99})$','interpreter','latex');
text(0.6*1/n_k2(end),err2_k2(end) ,'$\mathcal{O}(h^{0.93})$','interpreter','latex');
text(0.6*1/n_k3(end),err2_k3(end) ,'$\mathcal{O}(h^{1.98})$','interpreter','latex');
text(0.6*1/n_k4(end),err2_k4(end) ,'$\mathcal{O}(h^{2.99})$','interpreter','latex');

legend('$\bf{P}_2$','$\bf{P}_3$','$\bf{P}_4$','interpreter','latex')

xlabel('$h$', 'interpreter','latex')
ylabel('$|| \det u - \det  u_{exct} ||_{L^2}$', 'interpreter','latex')


