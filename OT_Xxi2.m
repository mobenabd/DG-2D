function [x, y] = OT_Xxi2(xi1, xi2)
%global ss
global ss0
ss = 0.6   + ss0;

delta1 = (1+ss)^2 - 4*ss*xi1;
delta2 = (1+ss)^2 - 4*ss*xi2;


x = ((1+ss) - sqrt(delta1))/(2*ss);
y = ((1+ss) - sqrt(delta2))/(2*ss);

end