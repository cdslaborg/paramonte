function str = contents()
    %
    %   Return the contents of the ``.bash_profile`` file 
    %   in the system home folder as a scalar MATLAB string.
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
    %           The output scalar MATLAB string containing the contents
    %           of the ``.bash_profile`` file if it exists or an empty string.
    %
    %   Interface
    %   ---------
    %
    %       str = pm.sys.bash_profile.contents()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    path = fullfile(pm.sys.path.home(), ".bash_profile");
    if  isfile(path)
        str = string(fileread(path));
    else
        str = "";
    end
end