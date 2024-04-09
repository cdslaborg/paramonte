function hasit = riskman()
    %
    %   Return a scalar MATLAB logical that is ``true`` if and
    %   only if the current installation of MATLAB contains
    %   the MATLAB Risk Management Toolbox.
    %
    %   This function searches the MATLAB license
    %   for an installation of the Toolbox.
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
    %       hasit = pm.matlab.has.riskman();
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    hasit = license('test', 'Risk_Management_Toolbox');
end