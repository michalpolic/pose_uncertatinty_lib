function write_corresp_to_file( file, corresp, model )

    fileID = fopen(file,'w');
    fprintf(fileID,'%f ',model);
    fprintf(fileID,'\n');
    for i = 1:size(corresp,1)
        fprintf(fileID,'%f ',corresp(i,:));
        fprintf(fileID,'\n');
        
        % TODO: output if we know the covariance matrices
    end
    fclose(fileID);
end

