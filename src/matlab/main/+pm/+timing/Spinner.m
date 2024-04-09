classdef Spinner < pm.matlab.Handle
    %
    %   This is the base class for generating objects
    %   that can display the time spinner on the console.
    %
    %   \devnote
    %
    %       The ``handle`` superclass is essential to allow
    %       object modification by the object methods.
    %
    %   Parameters
    %   ----------
    %
    %       See the documentation of the class constructor.
    %
    %   Attributes
    %   ----------
    %
    %       See below for information on the attributes (properties).
    %
    %   Methods
    %   -------
    %
    %       See below for information on the methods.
    %
    %   Returns
    %   -------
    %
    %       An object of class ``pm.timing.Spinner``.
    %
    %   Interface
    %   ---------
    %
    %       self = pm.timing.Spinner()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    properties(Access = protected)
        %
        %   tickmarks
        %
        %       The MATLAB ``char`` vector containing the set of
        %       characters that represent the passage of time in the spinner.
        %
        tickmarks = '|/-\';
    end

    properties(Hidden)
        %
        %   tickCount
        %
        %       The MATLAB integer containing the length of ``tickmarks``.
        %
        tickCount = 4;
        %
        %   format
        %
        %       The MATLAB ``char`` vector containing the spinner display format.
        %
        format = [repmat('\b', 1, 4 + 1), '%s'];
        %
        %   clock
        %
        %       The scalar integer representing the index of the ``tickmarks`` attribute.
        %
        clock = 0;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public)

        function self = Spinner()
            %
            %   Return a scalar object of class ``pm.timing.Spinner``.
            %
            %   This is the constructor of the class ``pm.timing.Spinner``.
            %
            %   Parameters
            %   ----------
            %
            %       tickmarks
            %
            %       The input MATLAB ``char`` vector containing the set of
            %       characters that represent the passage of time in the spinner.
            %       (**optional**, default = ``'|/-\'``)
            %
            %   Returns
            %   -------
            %
            %       self
            %
            %           The output scalar object of class ``pm.timing.Spinner``.
            %
            %   Interface
            %   ---------
            %
            %       self = pm.timing.Spinner()
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            self.tickmarks = '|/-\';
            self.tickCount = length(self.tickmarks);
            self.format = [repmat('\b', 1, self.tickCount + 3), '%s'];
            self.clock = 0;
        end

        function spin(self, fraction)
            %
            %   Rotate the tick mark of the spinner
            %   and display the percentage value of the input fraction.
            %
            %   This is a dynamic method of the class ``pm.timing.Spinner()``.
            %
            %   Parameters
            %   ----------
            %
            %       fraction
            %
            %           The input scalar MATLAB fractional real (``0 <= fraction <= 1``)
            %           representing the fraction of work so far accomplished.
            %
            %   Returns
            %   -------
            %
            %       None
            %
            %   Interface
            %   ---------
            %
            %       spinner = pm.timing.Spinner()
            %       spinner.spin(fraction)
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            if  0 < self.clock
                if  fraction < 1
                    self.clock = mod(self.clock, self.tickCount) + 1;
                else
                    self.clock = 1;
                end
                fprintf(self.format, [self.tickmarks(self.clock), sprintf(' %3.0f', 100 * fraction), '% ']);
            else
                self.clock = self.clock + 1;
                fprintf('%s', [self.tickmarks(self.clock), sprintf(' %3.0f', 100 * fraction), '% ']);
            end
            %if  1 <= fraction
            %    %fprintf('\b\b\b\b\b');
            %    fprintf('\n');
            %end
        end

    end

end