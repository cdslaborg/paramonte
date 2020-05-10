function bashrcContents = getBashrcContents()
    % returns the contents of the .bashrc file as a character vector
    bashrcPath = fullfile(getHomePath(), ".bashrc");
    if  exist(bashrcPath,"file")
        bashrcContents = fileread(bashrcPath);
    else
        fid = fopen(bashrcPath, 'w'); fclose(fid);
        bashrcContents = '';
    end
end
