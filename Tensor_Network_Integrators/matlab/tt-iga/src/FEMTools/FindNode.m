function [index,x,y] = FindNode(x,y,ctpx,ctpy)
M = max(max(ctpx(:)),max(ctpy(:)));
tor = M*10^-6;
n = length(ctpx);
    for i=1:n
        if abs(ctpx(i)-x)< tor && abs(ctpy(i)-y)< tor
            index=i;
        end
    end
end

