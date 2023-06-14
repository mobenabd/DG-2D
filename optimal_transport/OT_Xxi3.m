function x = OT_Xxi3(xi, ss)
%%%%% compute X(xi) known (the symetric test case) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta = (1+ss)^2 - 4*ss*xi;

x = ((1+ss) - sqrt(delta))/(2*ss);

end