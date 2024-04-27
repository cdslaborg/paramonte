function str = arch()
    %
    %   Return the processor architecture name of
    %   the current system on which MATLAB is running,
    %   which is also supported by the ParaMonte MATLAB library.
    %
    %   Parameters
    %   ----------
    %
    %       None
    %
    %   Returns
    %   -------
    %
    %       str
    %
    %           The output scalar MATLAB string containing the processor architecture.
    %           The output value can be either:
    %
    %               arm64
    %
    %                   This can be currently true only if the OS is Darwin.
    %
    %               amd64
    %
    %                   This can be currently true for either Darwin, Linux, or Windows.
    %
    %               If the function cannot detect any of the above,
    %               the system-returned value will be output as is.
    %
    %   Interface
    %   ---------
    %
    %       str = pm.sys.arch()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    if isunix
        [failed, str] = system("uname -m");
        if failed
            error   ( newline ...
                    + "Failed to detect the processor architecture via" + newline ...
                    + "the MATLAB intrinsic call `system('uname -m'): " + newline ...
                    + newline ...
                    + pm.io.arch + strtrim(str) + newline ...
                    + newline ...
                    )
        end
    elseif ispc
         [~, str] = system("echo %PROCESSOR_ARCHITECTURE%");
    else
        error   ( newline ...
                + "Unrecognized operating system for MATLAB." + newline ...
                + "This is impossible, unless we live in the distant future." + newline ...
                + newline ...
                );
    end
    strLower = lower(strtrim(str));
    if contains(strLower, "arm64")
        str = "arm64";
    elseif contains(strLower, "amd64") || contains(strLower, "x86_64")
        str = "amd64";
    else
        warning ( newline ...
                + "Unrecognized unsupported processor architecture" + newline ...
                + "for the ParaMonte MATLAB library: " + newline ...
                + newline ...
                + pm.io.tab() + strtrim(str) + newline ...
                + newline ...
                );
    end
    %str = computer('arch');
end