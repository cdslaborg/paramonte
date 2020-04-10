function Err = copyFile(pathOld, pathNew, isWindows)
    
    FUNCTION_NAME   = "@copyFile()";
    
    Err             = Err_class();
    Err.occurred    = false;

    if length(strtrim(pathOld)) == 0, return; end

    % First check whether file exists:
    
    fileExists      = isfile(pathNew);   % check if the file already exists
    
    if fileExists
        Err.occurred    = true;
        Err.msg         = FUNCTION_NAME + ": The requested copy file = '" + pathNew + "' already exists.";
        return
    end

    if isWindows
        cmd = 'copy "'  + pathOld + '" "'   + pathNew + '" > nul';
    else
        cmd = "cp "     + pathOld + " "     + pathNew;
    end
    
    counter = 0;
    while true
        counter = counter + 1;
        Err = SystemInfo_class.executeCmd(cmd);
        if Err.occurred
            Err.msg = FUNCTION_NAME + ": Error occurred while executing command " + cmd + "'.\n";
            return
        end
        % ensure file is copied
        fileExists      = isfile(pathNew);   % check if the file already exists
        if ~fileExists && counter < 100, continue; end
        break
    end
    
    if ~fileExists
        Err.occurred   = true;
        Err.msg        = FUNCTION_NAME + ": Failed to copy file from '" + pathOld + "' to '" + pathNew + "' after " + num2str(counter) + " attempts.";
        return
    end
    
end