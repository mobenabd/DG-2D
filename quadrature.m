function  val =  quadrature(f, dim)
%%%% Compute integrals over [-1,1]*[-1,1] for dim=2 or [-1,1] for dim =1 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%gauss-legendre nodes and weights (n=5) -> excat for poylnomials up to degree 9
% w = [0.5688888888888889, 0.4786286704993665, 0.4786286704993665, ...
%     0.2369268850561891, 0.2369268850561891];
% s = [0, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, ...
%     0.9061798459386640];


s = [0, 1/21*sqrt(245-14*sqrt(70)), -1/21*sqrt(245-14*sqrt(70)),...
        1/21*sqrt(245+14*sqrt(70)), -1/21*sqrt(245+14*sqrt(70))];

w = [128/225, 1/900*(322+13*sqrt(70)), 1/900*(322+13*sqrt(70)),...
        1/900*(322-13*sqrt(70)), 1/900*(322-13*sqrt(70))];    
   

val = 0;
if (dim==2)
    for i=1:5
        for j=1:5
            val = val + w(i)*w(j)*f(s(i),s(j));
        end
    end
elseif (dim==1)
    for i=1:5
        val = val + w(i)*f(s(i));
    end
end

end