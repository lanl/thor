function save_errors(filename,Lerror)
    %
    fid=fopen(filename,"w");
    %
    for eq=1:size(Lerror,1)
        for lvl=1:size(Lerror,2)
            fprintf(fid,"%e\t",Lerror(eq,lvl));
        end
        fprintf(fid,"\n");
    end
    %
    fclose(fid);
    %
end