function t = NonUniformKnot(nn,p)
n = nn - 1;
m = n + p + 1;
t = zeros(1,m+1);

for j=2:(n - p + 1)
 t(j + p) = (j - 1)/(n - p + 1);
end;

for k=(m - p + 1):(m + 1)
 t(k) = 1;
end

end

