function [index,x,y] = FindNodeY(y0,ctpx,ctpy)
tor= 10^-8;
index=[];
x=[];
y=[];
n = length(ctpx);
    for i=1:n
         if abs(ctpy(i)-y0)<tor
            index = [index i];
            x     = [x ctpx(i)];
            y     = [y ctpy(i)];
        end
    end
end
