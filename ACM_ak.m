function Ak = ACM_ak(k,c,U,LAP, f_projected,g_projected,num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%return local a^k(x) on quadrature points  for ACM method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloc = (k+1)^2;

Ak = zeros(Nloc,1);

idx = (num-1)*Nloc;

for i=1:Nloc
    if ( LAP(idx+i) + f_projected(idx+i) + c*( U(idx+i) - g_projected(idx+i) ) < 0 )
        Ak(i) = c;
    end
end


end