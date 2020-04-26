function Err = executeCmd(command)

    FUNCTION_NAME   = SystemInfo_class.CLASS_NAME + "@executeCmd()";
    
    command         = convertStringsToChars(command);
    Err             = Err_class();

    [Err.stat, cmdout] = system(command);

    if Err.stat == 0
        Err.occurred    = false;
        Err.msg         = cmdout;
        return
    elseif Err.stat == 1
        Err.occurred    = true;
        Err.msg         = FUNCTION_NAME + ": Error occurred. The processor does not support command execution of the command: "...
                        + command + ". " + newline + strtrim(cmdout);
        return
    else
        Err.occurred = true;
        Err.msg = FUNCTION_NAME + ": Unknown error occurred while attempting to execute the command: "...
                + command + ". The compiler/processor's explanatory message: " + newline + strtrim(cmdout);
        return
    end

end