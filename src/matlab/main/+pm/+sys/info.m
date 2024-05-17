function [str, cache] = info()
    %
    %   Return a MATLAB string containing the current system information.
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
    %           The output scalar MATLAB string containing the current system information.
    %
    %       cache
    %
    %           The output scalar MATLAB string representing the path to
    %           the cache file containing the current system information.
    %           The returned cache file path has the form:
    %
    %               <DIR>/.info<YEAR><MONTH>.cache
    %
    %           where,
    %
    %               <DIR> is replaced by the directory containing the function ``pm.sys.info.cache()``,
    %               <YEAR> is replaced by the current year,
    %               <MONTH> is replaced by the current month.
    %
    %           This means that the cache file contents are supposed to be updated only every month if necessary.
    %           A cache file is generated because retrieving system information is an expensive shell command-line operation.
    %           The time expense is particularly notable on Windows systems.
    %
    %   Interface
    %   ---------
    %
    %       [str, cache] = pm.sys.info()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
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