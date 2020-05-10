function bashProfileContents = getBashProfileContents()
    % returns the contents of the .bash_profile file as a character vector.
    % also adds the command to read the .bashrc file if it does not exist already.
    bashProfilePath = fullfile(getHomePath(), ".bash_profile");
    bashProfileContents = '';
    if exist(bashProfilePath,"file")
        bashProfileContents = fileread(bashProfilePath);
    end
    if ~contains(bashProfileContents,'.bashrc')
        fid = fopen(bashProfilePath, 'a+');
        fprintf(fid, '%s', '\n[ -f $HOME/.bashrc ] && . $HOME/.bashrc\n');
        fclose(fid);
        bashProfileContents = fileread(bashProfilePath);
    end
end