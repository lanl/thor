function tt_write(filename,tt)
  %
  fid = fopen(filename,"w");
  %
  fwrite(fid,tt.d,"int");
  fwrite(fid,tt.n,"int");
  fwrite(fid,tt.r,"int");
  fwrite(fid,tt.ps,"int");
  fwrite(fid,tt.core,"double");
  %
  fclose(fid);
  %
end