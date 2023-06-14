function PROJ = OT_compute_d_rho0Tilde_xi(n,k, rho1_proj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Compute  \Tilde{\rho_0}(xi) directly using \Phi, \delta\Phi && X(xi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloc = (k+1)^2;
gdim = n^2*Nloc;
PROJ = zeros(gdim,1);
[nodes, ~] = getWeightsNodes(k+1);

global ss0
ss= 0.6; ss0old = ss0;
ss_eps = ss + ss0;

%%% compute \nabla_xi\delta\Phi(xi)
func_dPhi_xi1 = @(xi, dump) OT_Xxi3(xi, ss_eps) - OT_Xxi3(xi, ss);
func_dPhi_xi2 = @(dump, xi) OT_Xxi3(xi, ss_eps) - OT_Xxi3(xi, ss);
GRAD_dPhi_xi1 =  computeDirectProjection(n,k, func_dPhi_xi1);  
GRAD_dPhi_xi2 =  computeDirectProjection(n,k, func_dPhi_xi2);  

%%% compute \nabla_x\rho_1(xi) && rho_1(xi)
GRAD_rho1_x = computeGrad_ddl(n,k,rho1_proj);
ss0 = 0;
[GRAD_rho1_xi1, ~] = OT_compute_rho0Tilde_xi_TEST(n,k, GRAD_rho1_x(:,1)); 
[GRAD_rho1_xi2, ~] = OT_compute_rho0Tilde_xi_TEST(n,k, GRAD_rho1_x(:,2)); 
[rho1_xi, ~] = OT_compute_rho0Tilde_xi_TEST(n,k, rho1_proj); 
ss0 = ss0old;



kk = 0;
for elem=1:n^2
    for i=1:Nloc
        ie = i + kk;
        
        if mod(i, k+1) ==0
            ii = k;
        else
            ii =  mod(i, k+1) -1;
        end
        jj = i/(k+1) - (ii+1)/(k+1);
        
        [xi1, xi2] = mapp_xy(nodes(ii+1), nodes(jj+1), elem, n);
        
        %%% compute tr((I+\nabla^2\Phi)^-1 * \nabla^2\delta\Phi)
        tr_xi = OT_Jac_Xxi(xi1, ss_eps)*OT_Jac_Xxi(xi1, ss)^(-1) ...
              + OT_Jac_Xxi(xi2, ss_eps)*OT_Jac_Xxi(xi2, ss)^(-1) ...
              - 2;
        
       
        %%% compute det(I+\nabla^2 Phi)
        det_xi = OT_Jac_Xxi(xi1, ss) * OT_Jac_Xxi(xi2, ss);
        
        
        %%% compute \delta\rho_0 (xi)
        PROJ(ie) = det_xi * ( GRAD_rho1_xi1(ie)*GRAD_dPhi_xi1(ie) ...
                            + GRAD_rho1_xi2(ie)*GRAD_dPhi_xi2(ie) ...
                            + rho1_xi(ie)*tr_xi);
           
    end
    kk = kk + Nloc;
end

end
