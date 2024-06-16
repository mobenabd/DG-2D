function Ak = ACM_ak(k,c,U,LAP, f_projected,g_projected,num, DET, Eig_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%return local a^k(x) on quadrature points  for ACM method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloc = (k+1)^2;

Ak = zeros(Nloc,1);

idx = (num-1)*Nloc;


if  (Eig_m == 0) %%solve obstacle problem/concavization
    for i=1:Nloc
        if ( LAP(idx+i) + f_projected(idx+i) + c*( U(idx+i) - g_projected(idx+i) ) < 0 ...%&& DET(idx+i) > 0)
               )% && -DET(idx+i) + f_projected(idx+i) + c*( U(idx+i) - g_projected(idx+i))  < 0 )
            
            Ak(i) = c;
            
        end
    end
    

elseif (Eig_m == 1) %%regularize solution
    for i=1:Nloc
        if ( -DET(idx+i) + f_projected(idx+i) + c*( U(idx+i) - g_projected(idx+i))  < 0 )
            
            Ak(i) = c;
            
        end
    end

end


end