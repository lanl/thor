function [index,x,y] = FindNodeAY(x0,y0,y1,ctpx,ctpy)
tor= 10^-3*max(ctpx(:));
index=[];
x=[];
y=[];
n = length(ctpx);
    for i=1:n
         if (abs(ctpx(i)-x0)<tor && ctpy(i)>=y0 && ctpy(i)<=y1) 
            index = [index i];
            x     = [x ctpx(i)];
            y     = [y ctpy(i)];
        end
    end
end