function k = sref(step,a)
if  length(step) ~= length(a)+1
    fprintf('Wrong Input Data \n')
end
a = sort([a 0 1]);
L = length(step);
k = [];
for i=1:L
k = [k (a(i+1)-a(i))*href1(step(i))+a(i)];
end

end

