function nodeType = boundary_position(N1,n)
if (N1<=n+1)
    nodeType = 'bottom';
    if (N1==n+1)
        nodeType = 'bottom-right';
    elseif(N1==1)
        nodeType = 'bottom-left';
    end
elseif (N1>=(n+1)^2-n)
    nodeType = 'top';
    if (N1==(n+1)^2)
        nodeType = 'top-right';
    elseif(N1==(n+1)^2-n)
        nodeType = 'top-left';
    end
elseif (mod(N1,n+1)==0 && N1~=n && N1~=(n+1)^2)
    nodeType = 'right';
elseif (mod(N1-1,n+1)==0 && N1~=1 && N1~=(n+1)^2-n)
    nodeType = 'left';
else
    nodeType = 'interior';
end
end