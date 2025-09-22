function print_stats(fid,opinf,type,xp,x_sim,i,n)
  %
  if(i==1)
    fprintf(fid,"|-------|------------|--------------|----------------|-------------|----------------------------------|\n");
    fprintf(fid,"| OpInf |  Tpod (s)  |  t_learn (s) |  t_predict (s) |  rel error  | tt: rank/POD: n, retained energy |\n");
    fprintf(fid,"|-------|------------|--------------|----------------|-------------|----------------------------------|\n");
  end
  %
  fmt = "|";
  fmt = fmt +  "%6s | %10s |   %.2e   |    %.2e    |   %.2e  |";
  %
  if(type=="tt")
    ranks = "          %3d, %3d, %3d          ";
    fmt = fmt + ranks + " \n";
    ranks = opinf.xp.r(2:4)';
  elseif(type=="qtt")
    ranks = "";
    for d=2:opinf.dims
      ranks = ranks + "%3d";
      if(d<opinf.dims)
        ranks = ranks + ",";
      end
    end
    fmt = fmt + ranks + " \n";
    ranks = opinf.ranks';
  elseif(type=="rom" || type=="ttrom")
    ranks = "         %3d, %8.4f %%         ";
    fmt = fmt + ranks + " \n";
    ranks = [opinf.n 100*opinf.retained_energy];
  else
    fmt = fmt + "                 -                \n";
    ranks = [];
  end
  %
  relErr = norm(x_sim(:)-xp(:)) / norm(x_sim(:));
  %
  if(type~="ft")
    tpod_str = sprintf('%.2e',opinf.Tpod);
  else
    tpod_str = '-';
  end
  %
  if(isempty(ranks))
    fprintf(fid,fmt, upper(type), tpod_str, opinf.Tlearn, opinf.Tpredict, relErr);
  else
    fprintf(fid,fmt, upper(type), tpod_str, opinf.Tlearn, opinf.Tpredict, relErr, ranks);
  end
  %
  if(i==n)
    fprintf(fid,"|-------|------------|--------------|----------------|-------------|----------------------------------|\n");
  end
  %
end
