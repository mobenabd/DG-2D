function val = integrate2(func, n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Integrate a function f(x,y) defined on [0,1]^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[~,~, detJac] = jacobian_elem(1,n);    % (uniforme grid)

val = 0;
for ielm=1:n^2
    source  =  func_mapping(func,n,ielm);
    
    val = val + quadrature(@(x,y) source(x,y),2);
    
end

val = val * detJac;

end
