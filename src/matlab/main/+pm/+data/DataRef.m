classdef DataRef < pm.matlab.Handle
    %
    %   This is the abstract class for generating instances of objects
    %   that can contain basic attributes required for storing reference
    %   to a given data.
    %
    %   This class is merely a convenience read-only
    %   wrapper to reference external data.
    %
    %   This class primarily exist to facilitate bypassing
    %   the lack of references and pointers in MATLAB.
    %
    %   Parameters
    %   ----------
    %
    %       dfref
    %
    %           The input object containing the target dataset
    %           or function handle that takes no arguments and returns the dataset.
    %           Specifying a function handle is superior to specifying the dataset
    %           directly, because the function handle will always allow the use of
    %           the most updated version of the user table or matrix.
    %           (**optional**. default = ``[]``)
    %
    %   Returns
    %   -------
    %
    %       self
    %
    %           The output scalar object of class ``pm.data.DataRef``.
    %
    %   Interface
    %   ---------
    %
    %       df = pm.data.DataRef(dfref);
    %
    %   Attributes
    %   ----------
    %
    %       See the class attributes descriptions below.
    %
    properties(Access = protected, Hidden)
        %
        %   dfref
        %
        %       A ``protected`` and ``Hidden`` copy of the user-specified
        %       input data (or its handle) ``dfref`` to the class constructor.
        %       This is an exact copy of the user-specified function handle or data.
        %       Users are supposed to access this component only
        %       via the class method ``pm.data.DataRef.copy()``.
        %
        dfref;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = DataRef(dfref)
            if  nargin < 1
                dfref = [];
            end
            self.dfref = dfref;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function df = copy(self)
            %
            %   Generate and return a copy of the data (or what is returned by the function handle)
            %   given the constructor of the parent object or the data return by the input.
            %
            %   This class method offer the only way to access the user-specified data (reference).
            %   The underlying logic behind the use of function to access the data (reference)
            %   originates from the lack of the concept of references (pointers)
            %   in the MATLAB computing language.
            %
            %   Parameters
            %   ----------
            %
            %       None
            %
            %   Returns
            %   -------
            %
            %       df
            %
            %           The output scalar MATLAB table a full copy of the data (reference)
            %           contained in the user-specified input ``dfref`` passed
            %           to the constructor of the parent object.
            %
            %   Interface
            %   ---------
            %
            %       df = pm.data.DataRef.copy()
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            if  isa(self.dfref, 'function_handle')
                df = self.dfref();
            else
                df = self.dfref;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end