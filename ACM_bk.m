function Bk = ACM_bk(k,c,U,LAP, f_projected,g_projected,num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%return local b^k(x) on quadrature points  for ACM method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloc = (k+1)^2;
Bk = zeros(Nloc,1);
idx = (num-1)*Nloc;

for i=1:Nloc
    if ( LAP(idx+i) + f_projected(idx+i) + c*( U(idx+i) - g_projected(idx+i) ) < 0 )
        Bk(i) = c*g_projected(idx+i) - LAP(idx+i);
    else
        Bk(i) =  f_projected(idx+i);
        %Bk(i) =  10;%DET(idx+i);
    end
end


end