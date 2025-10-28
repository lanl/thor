function tt_write_ascii(tt, fname)
% tt_write_ascii: writes ASCII file
%   - tt   : tensor train
%   - fname: file name to write
%

  fid = fopen(fname,'w');
  d = size(tt.n, 1);
  for k = 1:d-1
      fprintf(fid, '%d,', tt.n(k));
  end
  fprintf(fid, '%d\n', tt.n(d));
  for k = 1:d
      fprintf(fid, '%d,', tt.r(k));
  end
  fprintf(fid, '%d\n', tt.r(d+1));
  for k = 1:mem(tt)
      fprintf(fid, '%22.14e\n', tt.core(k));
  end
  fclose(fid);
end
