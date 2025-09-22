function [index,x,y] = FindNodeAX(x0,x1,y0,ctpx,ctpy)
tor = 10^-3*max(ctpy(:));
index=[];
x=[];
y=[];
n = length(ctpx);
    for i=1:n
         if (abs(ctpy(i)-y0)<tor && ctpx(i)>=x0 && ctpx(i)<=x1) 
            index = [index i];
            x     = [x ctpx(i)];
            y     = [y ctpy(i)];
        end
    end
end