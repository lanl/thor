function [x,Tpod] = get_snapshots_from_file_with_cross(cfg,eps_tt,j2)
  %
  fid     = fopen(cfg.filename,"r");
  cfg.N   = fread(fid,1,"int");
  cfg.Neq = fread(fid,1,"int");
  cfg.Nt  = fread(fid,1,"int");
  fclose(fid);
  %
  mmap = prepare_memmap(cfg);
  %
  m = cfg.N*cfg.Neq;
  n = j2;
  %
  t1=toc;
  if(cfg.crosstype == "greedy2")
    x = cross_interpolation([m,n],@(xx)funmap(xx,mmap),eps_tt,1,true);
  elseif(cfg.crosstype == "cross2d")
    %
    if(~isfield(cfg,"r0"))
      cfg.r0 = 2;
    end
    %
    x = cross2d_vec(@(xx) funmap(xx,mmap), m, n, eps_tt, "r0",cfg.r0);
  else
    error("Invalid cross interpolation option %s",cfg.crosstype)
  end
  %
  x = round(x,eps_tt);
  x = tt_orthogonolize(x);
  %
  t2 = toc;
  %
  Tpod = t2 - t1;
  %
end