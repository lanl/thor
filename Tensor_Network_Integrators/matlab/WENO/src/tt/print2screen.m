function print2screen(data,rki,is_final_step)
    %
    if(mod(data.time_step,1)==0 || data.time_step==1)
      %
      if(rki==0)
        %
        fmt = "\n%10s";
        %
        for n=1:3
          fmt = fmt + "%15s";
        end
        fmt = fmt + "%12s";
        fmt = fmt + "%12s";
        %
        for n=1:data.tt.Neq
          fmt = fmt + "%9s";
        end
        fmt = fmt + "\n";
        %
        if(data.tt.Neq==1)
          fprintf("%s",repelem("-",88));
          fprintf(fmt,"Time Step","Time (s)","dt (s)","max_eig (m/s)","eps", ...
                  " | TT_info: ","u");
          fprintf("%s",repelem("-",88));
        elseif(data.tt.Neq==5)
            fprintf("%s",repelem("-",124));
            fprintf(fmt,"Time Step","Time (s)","dt (s)","max_eig (m/s)","eps", ...
                    " | TT_info: ","rho","rhoU","rhoV","rhoW","rhoE");
            fprintf("%s",repelem("-",124));
        end
        %
        fprintf("\n");
        %
      elseif(rki==1)
        if(data.tt.Neq==1)
          fprintf("%s",repelem("..",44));
        elseif(data.tt.Neq==5)
          fprintf("%s",repelem("..",62));
        end
        fprintf("\n");
      end
      %
      if(rki==0)
        h32 = sqrt(data.vol);
        %
        fmt = "%10d";
        for n=1:3
          fmt = fmt + "%15.5e";
        end
        fmt = fmt + "%12.4e";
        fmt = fmt + "%12s";
        %
        fprintf(fmt,data.time_step,data.curr_time,data.dt,data.max_eig,data.tt.eps," | F-norm : ");
        for n=1:data.tt.Neq
            fprintf("%9.2e",norm(data.tt.Q{n})*h32);
        end
        fprintf("\n");
        %
      end
      %
      if(rki>0 || is_final_step)
        %
        if(rki==0)
          if(data.tt.Neq==1)
            fprintf("%s",repelem("..",44));
          elseif(data.tt.Neq==5)
            fprintf("%s",repelem("..",62));
          end
          fprintf("\n");
          fprintf("%55s%12s%12s","",""," | (r1,r2): ");
        else
          fprintf("%22s%9s%2d%22s%12s%12s","","rk-stage:",rki,"",""," | (r1,r2) : ");
        end
        %
        for n=1:data.tt.Neq
          if(data.tt.d==2)
            str = sprintf("%1s%1d%s","(",data.tt.Q{n}.r(2),")");
          else
            str = sprintf("%1s%1d%1s%1d%s","(",data.tt.Q{n}.r(2),",",data.tt.Q{n}.r(3),")");
          end
          fprintf("%9s",str);
        end
        %
        fprintf("\n");
        %
      end
      %
      if(rki==data.rkmax)
        fprintf("\n");
      end
      %
    end
    %
end