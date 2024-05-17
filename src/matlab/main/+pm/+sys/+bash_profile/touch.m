function str = touch()
    %
    %   Return the contents of the ``.bash_profile`` file
    %   in the system home folder as a scalar MATLAB string.
    %
    %   If the file does not exist, create it. and add to it
    %   If the file does not contain the Bash to source the
    %   contents of the ``.bashrc`` file, then add it to
    %   the ``.bash_profile`` and return its contents.
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
    %           The output scalar MATLAB string containing
    %           the contents of the ``.bash_profile``.
    %
    %   Interface
    %   ---------
    %
    %       str = pm.sys.bash_profile.touch()
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
    if ~contains(str, ".bashrc")
        fid = fopen(path, 'a+');
        fprintf(fid, '%s', [newline, '[ -f $HOME/.bashrc ] && . $HOME/.bashrc', newline]);
        fclose(fid);
        str = string(fileread(path));
    end
end