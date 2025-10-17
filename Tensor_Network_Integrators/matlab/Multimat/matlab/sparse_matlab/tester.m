% Tester script 

filenames = {'input_advection','input_blowoff','input_piston','input_sod'};

for i = 1:length(filenames)
    diaryname = [filenames{i},'_output.txt'];
    diary(diaryname)
    main(filenames{i})
    diary off
end 