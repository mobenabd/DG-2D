function val = Lagrange(x,i,k,derv)
%%%%%%%%%%%% Lagrange polynomials %%%%%%%%%%%%%%%%
%%%%%%%%%%%% derv=0: polynomial
%%%%%%%%%%%% derv=1: first deriviative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[nodes, ~] = getWeightsNodes(k+1);
l = length(nodes);


if derv ==0
    val = 1;
    for j=1:l
        if (j ~= i)
            val = val*(x-nodes(j))/(nodes(i) - nodes(j));
        end
    end
end


if derv==1
    val = 0;
    for j=1:l
        if (j ~= i)
            tmp = 1.0/(nodes(i) - nodes(j));
            for m=1:l
                if (m ~= i && m ~= j)
                    tmp = tmp*(x-nodes(m))/(nodes(i) - nodes(m));
                end
            end
            val = val + tmp;
        end
        
    end 
end

if derv==2
    val = 0;
    for j=1:l %somme
        if (j ~= i)
            tmp = 1/(nodes(i) - nodes(j));
            
            tmp2 = 0;
            for m=1:l  %%somme
                if (m ~= i && m ~= j)
                    tmp1 = 1/(nodes(i) - nodes(m)); 
                     
                    for nn=1:l  %%produit
                        if (nn ~= i && nn ~= j && nn ~= m)
                             tmp1 = tmp1.*(x-nodes(nn))/(nodes(i) - nodes(nn));
                        end
                    end
                   tmp2 = tmp2 + tmp1;
                end 
            end
            tmp = tmp*tmp2;
            
            val = val + tmp;
        end
        
    end
end

end