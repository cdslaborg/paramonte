%>  \brief
%>  Return the processor architecture name of
%>  the current system on which MATLAB is running,
%>  which is also supported by the ParaMonte MATLAB library.<br>
%>
%>  \return
%>  ``str`` :   The output scalar MATLAB string containing the processor architecture.<br>
%>              The output value can be either:<br>
%>              <ol>
%>                  <li>    ``"arm64"`` :   This can be currently true only if the OS is Darwin.<br>
%>                  <li>    ``"amd64"`` :   This can be currently true for either Darwin, Linux, or Windows.<br>
%>              </ol>
%>
%>  \interface{arch}
%>  \code{.m}
%>
%>      str = pm.sys.arch()
%>
%>  \endcode
%>
%>  \note
%>  If the function cannot detect any of the above,
%>  the system-returned value will be output as is.<br>
%>
%>  \example{arch}
%>  \include{lineno} example/sys/arch/main.m
%>  \output{arch}
%>  \include{lineno} example/sys/arch/main.out.m
%>
%>  \final{arch}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 5:34 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function str = arch()
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