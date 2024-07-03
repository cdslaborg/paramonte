%>  \brief
%>  This is the base class for generating objects
%>  that can display the time spinner on the console.
%>
%>  \devnote
%>  The ``handle`` superclass is essential to allow
%>  object modification by the object methods.
%>  See the documentation of the class constructor.
%>
%>  \note
%>  See below for information on the attributes (properties).
%>
%>  \note
%>  See below for information on the methods.
%>
%>  \return
%>  An object of class ``pm.timing.Spinner``.
%>
%>  \interface{Spinner}
%>  \code{.m}
%>
%>      self = pm.timing.Spinner()
%>
%>  \endcode
%>
%>  \final{Spinner}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 5:39 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef Spinner < pm.matlab.Handle

    properties(Access = protected)
        %>
        %>  \param[in]  tickmarks   :   The MATLAB ``char`` vector containing the set of
        %>                              characters that represent the passage of time in the spinner.
        %>
        tickmarks = '|/-\';
    end

    properties(Hidden)
        %>
        %>  \param  tickCount   :   The MATLAB integer containing the length of ``tickmarks``.
        %>
        tickCount = 4;
        %>
        %>  \param  format      :   The MATLAB ``char`` vector containing the spinner display format.
        %>
        format = [repmat('\b', 1, 4 + 1), '%s'];
        %>
        %>  \param  clock       :   The scalar integer representing the index of the ``tickmarks`` attribute.
        %>
        clock = 0;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %>  \brief
        %>  Return a scalar object of class ``pm.timing.Spinner``.<br>
        %>  This is the constructor of the class ``pm.timing.Spinner``.
        %>
        %>  \param[in]  tickmarks   :   The input MATLAB ``char`` vector containing the set of
        %>                              characters that represent the passage of time in the spinner.
        %>                              (**optional**, default = ``'|/-\'``)
        %>
        %>  \return
        %>  ``self``:   The output scalar object of class ``pm.timing.Spinner``.
        %>
        %>  \interface{Spinner}
        %>  \code{.m}
        %>
        %>      self = pm.timing.Spinner()
        %>
        %>  \endcode
        %>
        %>  \final{Spinner}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 5:43 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Spinner()
            self.tickmarks = '|/-\';
            self.tickCount = length(self.tickmarks);
            self.format = [repmat('\b', 1, self.tickCount + 3), '%s'];
            self.clock = 0;
        end

        %>  \brief
        %>  Rotate the tick mark of the spinner
        %>  and display the percentage value of the input fraction.<br>
        %>  This is a dynamic method of the class ``pm.timing.Spinner()``.
        %>
        %>  \param[in]  fraction    :   The input scalar MATLAB fractional real (``0 <= fraction <= 1``)
        %>                              representing the fraction of work so far accomplished.
        %>
        %>  \interface{spin}
        %>  \code{.m}
        %>
        %>      spinner = pm.timing.Spinner()
        %>      spinner.spin(fraction)
        %>
        %>  \endcode
        %>
        %>  \final{spin}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 5:45 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function spin(self, fraction)
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