function CC = KnotInsert(knot,t)


k      = Findk(knot,t);
[m, p] = EvaluateKnot(knot);
Alpha  = zeros(m + 1,1);
CC     = zeros(m,m + 1);


for i=1:(k - p)
 Alpha(i) = 1;
end

for i=(k - p + 1):k
 Alpha(i) = (t - knot(i))/(knot(i + p) - knot(i));
end

for i=(k + 1):(m + 1)
 Alpha(i) = 0;
end

for i=1:m
 CC(i, i)     = Alpha(i);
 CC(i, i + 1) = 1 - Alpha(i + 1);
end;

CC=CC';

end

