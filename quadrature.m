function  val =  quadrature(f, dim,k)
%%%% Compute integrals over [-1,1]*[-1,1] for dim=2 or [-1,1] for dim =1 %%
%%%% k here is the number of quadrature nodes. (used to integer source/initial condition function)
%%%% otherwise we take k+1 nodes where k is the degree of polynomials.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if  ~exist('k','var')
    k = getGlobal_k();
    [s, w] = getWeightsNodes(k+1);
else
    %k = getGlobal_k();
    %[s, w] = getWeightsNodes(k+1);
    [s, w] = getWeightsNodes(k);
end


l = length(s);

val = 0;
if (dim==2)
    for i=1:l
        for j=1:l
            val = val + w(i)*w(j)*f(s(i),s(j));
        end
    end
elseif (dim==1)
    for i=1:l
        val = val + w(i)*f(s(i));
    end
end

end