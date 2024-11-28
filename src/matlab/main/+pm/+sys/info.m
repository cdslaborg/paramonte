%>  \brief
%>  Return a MATLAB string  and the corresponding cache
%>  file path containing the current system information.<br>
%>
%>  \return
%>  ``txt``     :   The output scalar MATLAB string containing the current system information.<br>
%>  ``cache``   :   The output scalar MATLAB string representing the path to the
%>                  cache file containing the current system information.<br>
%>                  The returned cache file path has the form:<br>
%>                  \code{.m}
%>                      <DIR>/.info<YEAR><MONTH>.cache<br>
%>                  \endcode
%>                  where,<br>
%>                  <ol>
%>                      <li>    ``<DIR>`` is replaced by the directory containing the function [pm.sys.info()](@ref info),
%>                      <li>    ``<YEAR>`` is replaced by the current year,
%>                      <li>    ``<MONTH>`` is replaced by the current month.
%>                  </ol>
%>                  This means that the cache file contents are supposed to be updated only every month if necessary.<br>
%>                  A cache file is generated because retrieving system information is an expensive shell command-line operation.<br>
%>                  The time expense is particularly notable on Windows systems.<br>
%>
%>  \interface{info}
%>  \code{.m}
%>
%>      [txt, cache] = pm.sys.info()
%>
%>  \endcode
%>
%>  \example{info}
%>  \include{lineno} example/sys/info/main.m
%>  \output{info}
%>  \include{lineno} example/sys/info/main.out.m
%>
%>  \final{info}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 5:31 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function [txt, cache] = info()
    prefix = ".info";
    suffix = ".cache";
    [dirname, ~, ~] = fileparts(mfilename('fullpath'));
    cache = fullfile(string(dirname), prefix + string(year(date)) + string(month(date)) + suffix); % update cache once a month.
    if ~isfile(cache)
        delete(prefix + "*" + suffix);
        if ispc
            cmd = "systeminfo";
        end
        if ismac
            cmd = "uname -a; sysctl -a | grep machdep.cpu";
        end
        if isunix && ~ismac
            cmd = "uname -a; lscpu";
        end
        [failed, txt] = system(cmd);
        failed = failed ~= 0;
        if  failed
            warning ( newline ...
                    + "Failed to fetch the system information on the current system. skipping..." ...
                    + newline ...
                    );
        end
        fid = fopen(cache, 'wt');
        fprintf(fid, "%s", txt);
        fclose(fid);
    else
        txt = fileread(cache);
    end
    txt = string(txt);
end