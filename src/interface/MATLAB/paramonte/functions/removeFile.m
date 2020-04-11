function Err = removeFile(path,isWindows)
    
    FUNCTION_NAME   = "@removeFile()";
    
    Err             = Err_class();
    Err.occurred    = false;

    % First check whether file exists:
    
    fileExists      = isfile(path); % check if the file already exists
    
    if ~fileExists
        Err.occurred    = true;
        Err.msg         = FUNCTION_NAME + ": The requested file = '" + path + "' does not exist.";
        return
    end
    
    if isWindows
        cmd = "del " + path + " > nul";
    else
        cmd = "rm " + path;
    end
    
    counter = 0;
    while true
        counter = counter + 1;
        Err = executeCmd(cmd);
        if Err.occurred
            Err.msg = FUNCTION_NAME + ": Error occurred while executing command " + cmd + "'.\n";
            return
        end
        % ensure file is removed
        
        fileExists  = isfile(path); % check if the file already exists
        if fileExists
            Err.occurred    = true;
            Err.msg         = FUNCTION_NAME + ": Failed to remove file = '" + path + "' after " + num2str(counter) + " attempts.";
            return
        end
        if fileExists && counter < 100, continue; end
        break
    end
    
    if fileExists
        Err.occurred = true;
        Errmsg = FUNCTION_NAME + ": Failed to remove file = '" + path + "' after " + num2str(counter) + " attempts.";
        return
    end

end