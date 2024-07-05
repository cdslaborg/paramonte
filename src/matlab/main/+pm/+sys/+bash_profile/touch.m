%>  \brief
%>  Return the contents of the ``.bash_profile`` file
%>  in the system home folder as a scalar MATLAB string.<br>
%>
%>  \details
%>  If the file does not exist, create it.<br>
%>  If the file does not contain the Bash command to
%>  source the contents of the ``.bashrc`` file, then add
%>  it to the ``.bash_profile`` and return its contents.<br>
%>
%>  \return
%>  ``str`` :   The output scalar MATLAB string containing
%>              the contents of the ``.bash_profile``.<br>
%>
%>  \interface{touch}
%>  \code{.m}
%>
%>      str = pm.sys.bash_profile.touch()
%>
%>  \endcode
%>
%>  \example{touch}
%>  \include{lineno} example/sys/bash_profile/touch/main.m
%>  \output{touch}
%>  \include{lineno} example/sys/bash_profile/touch/main.out.m
%>
%>  \final{touch}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 4:58 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function str = touch()
    path = fullfile(pm.sys.path.home(), ".bash_profile");
    if  isfile(path)
        str = string(fileread(path));
    else
        str = "";
    end
    if ~contains(str, ".bashrc")
        fid = fopen(path, 'a+');
        fprintf(fid, '%s', [newline, '[ -f $HOME/.bashrc ] && . $HOME/.bashrc', newline]);
        fclose(fid);
        str = string(fileread(path));
    end
end