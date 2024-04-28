function hasit = parallel()
    %
    %   Return a scalar MATLAB logical that is ``true`` if and
    %   only if the current installation of MATLAB contains
    %   the MATLAB Parallel Computing Toolbox.
    %
    %   This function searches the MATLAB license
    %   for an installation of the Parallel Computing Toolbox.
    %   If the search fails, a parallel code section will be tested.
    %
    %   Parameters
    %   ----------
    %
    %       None
    %
    %   Returns
    %   -------
    %
    %       hasit
    %
    %           The output scalar MATLAB logical that is ``true`` if and
    %           only if the current installation of MATLAB contains
    %           the required MATLAB Toolbox.
    %
    %   Interface
    %   ---------
    %
    %       hasit = pm.matlab.has.parallelism();
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    hasit = license('test', 'Distrib_Computing_Toolbox');
    if ~hasit
        try
            delete(gcp('nocreate'));
            hasit = true;
        end
    end
end