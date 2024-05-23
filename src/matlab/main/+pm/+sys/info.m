%>  \brief
%>  Return a MATLAB string containing the current system information.
%>
%>  \param[in]  `None`
%>
%>  \return
%>  `str`   :   The output scalar MATLAB string containing the current system information.<br>
%>  `cache` :   The output scalar MATLAB string representing the path to
%>              the cache file containing the current system information.
%>              The returned cache file path has the form:<br>
%>                  <DIR>/.info<YEAR><MONTH>.cache<br>
%>          where,<br>
%>              <DIR> is replaced by the directory containing the function ``pm.sys.info.cache()``,
%>              <YEAR> is replaced by the current year,
%>              <MONTH> is replaced by the current month.<br>
%>          This means that the cache file contents are supposed to be updated only every month if necessary.
%>          A cache file is generated because retrieving system information is an expensive shell command-line operation.
%>          The time expense is particularly notable on Windows systems.
%>
%>  \interface{info}
%>  \code{.m}
%>
%>      [str, cache] = pm.sys.info()
%>
%>  \endcode
%>  \final{info}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 5:31 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function [str, cache] = info()
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
        [errorOccurred, str] = system(cmd);
        if errorOccurred
            warning ( newline ...
                    + "Failed to fetch the system information on the current system. skipping..." ...
                    + newline ...
                    );
        end
        fid = fopen(cache, 'wt');
        fprintf(fid, "%s", str);
        fclose(fid);
    else
        str = fileread(cache);
    end
    str = string(str);
end