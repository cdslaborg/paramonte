function matlabRelease = release(type)
    %
    %   Return a scalar MATLAB string containing the
    %   MATLAB release version, year, or season as requested.
    %
    %   Parameters
    %   ----------
    %
    %       type
    %
    %           The input scalar MATLAB string that can be either:
    %
    %               1.  ``"year"``, implying that the
    %                   MATLAB release year must be returned.
    %
    %               2.  ``"season"`` or ``"minor"``, implying
    %                   that the MATLAB release year must be returned.
    %
    %           (**optional**. default = ``""`` implying the full version.)
    %
    %   Returns
    %   -------
    %
    %       matlabRelease
    %
    %           The output scalar MATLAB string containing the
    %           MATLAB release version, year, or season as requested.
    %
    %   Interface
    %   ---------
    %
    %       matlabRelease = pm.matlab.release()
    %       matlabRelease = pm.matlab.release(type)
    %
    %   Example
    %   -------
    %
    %       matlabReleaseSeason = pm.matlab.release("minor");
    %       matlabReleaseSeason = pm.matlab.release("season");
    %       matlabReleaseYear = str2double(pm.matlab.release("year"));
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    if nargin < 1
        type = "";
    end
    % WARNING: Do NOT convert the input char '-release' to string or the output will be wrong.
    matlabRelease = version('-release');
    if string(type) == "year"
        matlabRelease = matlabRelease(1:4);
    elseif string(type) == "season" || string(type) == "minor"
        matlabRelease = matlabRelease(5:5);
    elseif type ~= ""
        disp("type = ");
        disp(type);
        error("Invalid input value for the argument ``type`` of getMatlabRelease().");
    end
    matlabRelease = string(matlabRelease);
end