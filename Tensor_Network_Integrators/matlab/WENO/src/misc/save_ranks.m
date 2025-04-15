function save_ranks(filename,rank)
    fid=fopen(filename,"w");
    Neq = (size(rank,1)-1)/2;
    for time_step = 1:size(rank,2)
        fprintf(fid,"%e\t",rank(1:2,time_step));
        for eq=1:Neq
            fprintf(fid,"%d\t",rank(2*eq+1:2*eq+2,time_step));
        end
        fprintf(fid,"\n");
    end
    fclose(fid);
end