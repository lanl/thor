function print2screen(data,rki,is_final_step)
    %
    if(mod(data.time_step,10)==0 || data.time_step==1) 
        %
        if(rki==0)
            %
            %
            fmt = "\n%10s";
            for n=1:3
                fmt = fmt + "%15s";
            end
            fmt = fmt + "%12s";
            fmt = fmt + "%12s";
            for n=1:data.tt.Neq
                fmt = fmt + "%9s";
            end
            fmt = fmt + "\n";
            %
            fprintf("%s",repelem("-",106));
            if(data.SWEtype == "linear")
                fprintf(fmt,"Time Step","Time (s)", "dt (s)", "max_eig (m/s)", "eps", ...
                        " | TT info: ", "eta", "u", "v");
            else
                fprintf(fmt,"Time Step","Time (s)", "dt (s)", "max_eig (m/s)", "eps", ...
                        " | TT info: ", "h", "hu", "hv");
            end
            fprintf("%s",repelem("-",106));
            %
            fprintf("\n"); 
        elseif(rki==1) 
            fprintf("%s",repelem("..",53));
            fprintf("\n");  
        end

        if(rki==0)
            h32 = sqrt(data.vol);
            %
            fmt = "%10d";
            fmt = fmt + "%15.2f";
            for n=1:2
                fmt = fmt + "%15.5e";
            end
            fmt = fmt + "%12.4e";
            fmt = fmt + "%12s";

            fprintf(fmt,data.time_step,data.ref.t*data.curr_time,data.ref.t*data.dt,data.ref.U*data.max_eig,data.tt.eps," | F-norm : ");
            for n=1:data.tt.Neq
                fprintf("%9.2e",norm(data.tt.Q{n})*h32);
            end
            fprintf("\n");  
            %
        end

        if(rki>0 || is_final_step)
            %
            fprintf("%22s%9s%2d%22s%12s%12s","","rk-stage:", rki,"","", " |       r: ");
            %
            for n=1:data.tt.Neq
                str=sprintf("%1d",data.tt.Q{n}.r(2));
                fprintf("%9s",str);
            end
            %
            fprintf("\n"); 
            %
        end
        %
        if(rki==3) 
            fprintf("\n");  
        end
    end
  %
end