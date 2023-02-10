clear
clc
close all





syms q ;
eqn = q.^2 .*(1-log(q/2)) == 1;
r = double(solve(eqn,q));
r = solve(eqn,q)



return
syms a b c d;

A = sym([a b; c d]);

[~, D] = eig(A);

D(1,1)*D(2,2) - det(A)