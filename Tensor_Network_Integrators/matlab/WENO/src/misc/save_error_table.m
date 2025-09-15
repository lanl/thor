function save_error_table(filename,error,error_kind,N)
    %
    fid = fopen(filename,"w");
    %
    print_errors(fid,error,error_kind,N);
    %
    fclose(fid);
    %
end