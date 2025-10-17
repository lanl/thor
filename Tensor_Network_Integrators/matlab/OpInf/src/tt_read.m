function tt = tt_read(filename)
  %
  fid = fopen(filename,"r");
  %
  tt = tt_tensor;
  %
  tt.d    = fread(fid,1,"int");
  tt.n    = fread(fid,tt.d,"int");
  tt.r    = fread(fid,tt.d+1,"int");
  tt.ps   = fread(fid,tt.d+1,"int");
  %
  tt.core = fread(fid,"double");
  %
  fclose(fid);
  %
end