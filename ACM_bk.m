function Bk = ACM_bk(k,c,U,LAP, f_projected,g_projected,num, DET, Eig_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%return local b^k(x) on quadrature points  for ACM method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloc = (k+1)^2;
Bk = zeros(Nloc,1);
idx = (num-1)*Nloc;

for i=1:Nloc
    if ( LAP(idx+i) + f_projected(idx+i) + c*( U(idx+i) - g_projected(idx+i) ) < 0 )
        Bk(i) = c*g_projected(idx+i) - LAP(idx+i);
        %Bk(i) = c*g_projected(idx+i) ;
    else
        Bk(i) = sqrt( LAP(idx+i)^2 + 2*( f_projected(idx+i)-DET(idx+i)) );
        %Bk(i) =  f_projected(idx+i);
    end
end


% for i=1:Nloc
% 
%     Bk(i) = sqrt( LAP(idx+i)^2 +2*( f_projected(idx+i)-DET(idx+i) ) );
% 
% end


end