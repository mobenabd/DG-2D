function  val =  ACM_quadrature(Ak,f)
%%%% Compute integrals over [-1,1]*[-1,1] %%%%%
%%%% Ak: b^k or a^k evaluated in quadrature points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = getGlobal_k();
[s, w] = getWeightsNodes(k+1);
l = k+1;

val = 0;
for i=1:l
    for j=1:l
        val = val + w(i)*w(j)*f(s(i),s(j)) *Ak(l*(j-1)+i);
    end
end


end