function Bk = ACM_bk(k,c,U,LAP, f_projected,g_projected,num, DET, Eig_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%return local b^k(x) on quadrature points  for ACM method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloc = (k+1)^2;
Bk = zeros(Nloc,1);
idx = (num-1)*Nloc;

if  (Eig_m == 0) %%solve obstacle problem/concavization
    for i=1:Nloc
        if ( LAP(idx+i) + f_projected(idx+i) + c*( U(idx+i) - g_projected(idx+i) ) < 0 ...%&& DET(idx+i) > 0)
                && -DET(idx+i) + f_projected(idx+i) + c*( U(idx+i) - g_projected(idx+i))  < 0 )
            
            Bk(i) = c*g_projected(idx+i) - LAP(idx+i);
        else
            Bk(i) = sqrt( LAP(idx+i)^2 + 2*( f_projected(idx+i)-DET(idx+i)) );
            %Bk(i) =  f_projected(idx+i);
        end
    end

elseif (Eig_m == 1) %%regularize solution
    for i=1:Nloc
        if ( -DET(idx+i) + f_projected(idx+i) + c*( U(idx+i) - g_projected(idx+i))  < 0 )
            
            Bk(i) = c*g_projected(idx+i) - LAP(idx+i);
        else
            Bk(i) = sqrt( LAP(idx+i)^2 + 2*( f_projected(idx+i)-DET(idx+i)) );
            %Bk(i) =  f_projected(idx+i);
        end
    end
    
end


%%%%%%%%% For Monge-Ampère Equation %%%%%%%%%%
% for i=1:Nloc
% 
%     Bk(i) = sqrt( LAP(idx+i)^2 +2*( f_projected(idx+i)-DET(idx+i) ) );
% 
% end


end