function f = cross_fun_nD(a,ftrue,nout)
% x is an input array of tensor train
if nargin<3
  nout = 1;
end
n = size(a,1); %number of samples
f = zeros(n,nout);
for i = 1:n
  c = num2cell(a(i,:));
  f(i,:) = ftrue(c{:});
end