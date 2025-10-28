function [tt] = tt_read_ascii(fname)
% tt_read_ascii: reads ASCII file
%   - fname: file name to read
%
% Output:
%   - tt : tt-tensor
  % determine the tensor rank (number of dimensions)
  fid = fopen(fname,'r');
  d = size(strsplit(fgetl(fid),','),2);
  fclose(fid);

  % read it
  fid = fopen(fname,'r');
  nn = fscanf(fid, '%d,', d);
  rr = fscanf(fid, '%d,', d+1);
  tt = tt_rand(nn, d, rr);
  tt.r = rr;
  tt.ps = cumsum([1;tt.n.*tt.r(1:d).*tt.r(2:d+1)]);
  tt.core = fscanf(fid, '%e', mem(tt));
end
