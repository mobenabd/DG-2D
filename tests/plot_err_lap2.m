clear
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Verifier l'ordre de convergence sur le laplacien au sens faible %%%%
%%%%% on résount u_t-Δu = f puis on calcule le laplacien projeté  %%%%%%%%
%%%%% la solution exacte est Δu = u_t - f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% La base de Lagrange sur les points de quadrature %%%%%%%%%%%%%%%%%%%

%%% u_exct(x,y,t) = sin(2pix)*sin(2piy)*cos(pi/2*t)
%%% dt = 0.001
%%% T_final = 100*dt
%%% eps=-1, sigma =100


%%%  k=1, SIPG, 100
n_k1 = [4 8 16 32];
err2_k1 =[0.1159657536 0.0297703436  0.0075080196 0.0018858308]; %w.r.t u
err2_k1 =[5.0344174531 1.2737923170 0.3192593978 0.0798646470]; %w.r.t lap(u)

%%% k=1, NIPG,100
err2_k1 =[0.1134104987 0.0288720069 0.0072799070 0.0018302460]; %w.r.t u
err2_k1 =[5.0344845413 1.2737961324 0.3192595019 0.0798645679]; %w.r.t lap(u)

coeff_k1 = polyfit(log(1./n_k1),log(err2_k1),1)


%%%  k=2, SIPG, 100
n_k2 = [4 8 16 24];
err2_k2 =[0.0138114078 0.0018370021 0.0002337705 0.0000701028]; %w.r.t u
err2_k2 =[0.7012633612 0.0905404813 0.0114290172 0.0034411661]; %w.r.t lap(u)

%%%  k=2, NIPG, 100
%err2_k2 = [0.0144259704 0.7012445873 0.0002848207 0.0000994703];
%err2_k2 = [0.7012445873 0.0905400575 0.0114300260 0.0034428576];


coeff_k2 = polyfit(log(1./n_k2),log(err2_k2),1)
