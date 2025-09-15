function print_errors(fid,error,error_kind,N)
    %
    fprintf(fid,"\n\n%s errors:\n",error_kind);
    %
    maxline=10;
    %
    for cnt=1:size(error,2)
        if(size(N,2)==1) 
            Nx = N(cnt,1);
            grid=sprintf("Nx=%d",Nx);
            maxline=max(maxline,max(ceil(log10(N)),[],"all")+5);
        elseif(size(N,2)==2)
            Nx = N(cnt,1);
            Ny = N(cnt,2);
            grid=sprintf("%dx%d",Nx,Ny);
            maxline=max(maxline,2*max(ceil(log10(N)),[],"all")+3);
        else
            Nx = N(cnt,1);
            Ny = N(cnt,2);
            Nz = N(cnt,3);
            grid=sprintf("%dx%dx%d",Nx,Ny,Nz);
            maxline=max(maxline,3*max(ceil(log10(N)),[],"all")+4);
        end
        a=length(grid{1});
        b=floor((maxline-a)/2);
        c=maxline-a-b;
        %
        text = sprintf("%s%ds%s%ds%s%ds","%",b,"%",a,"%",c);
        text = sprintf(" | %s",text);
        %
        fprintf(fid,text,"",grid,"");
        %
        if(cnt>1)
            fprintf(fid," | ORDER");
        end
        %
    end
    %
    fprintf(fid," |\n ");
    %
    for cnt=1:size(error,2)
        fprintf(fid,"|-");
        for jj=1:maxline
            fprintf(fid,"-");
        end
        fprintf(fid,"-");
        if(cnt>1)
            fprintf(fid,"|-------");
        end
    end
    fprintf(fid,"|\n");
    %
    a=8;
    b=floor((maxline-a)/2);
    c=maxline-a-b;
    %
    format=sprintf(" |%s%ds%s%s%ds","%",b,"%.3e","%",c);
    %
    for eq=1:size(error,1)
        for cnt=1:size(error,2)
            fprintf(fid,format,"",error(eq,cnt),"");
            if(cnt>1)
                order=log(error(eq,cnt-1)/error(eq,cnt))/log(N(cnt,1)/N(cnt-1,1));
                fprintf(fid," | %.3f",order);
            end
        end
        fprintf(fid," |\n");
    end
    %
end