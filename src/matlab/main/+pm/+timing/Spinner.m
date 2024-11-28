%>  \brief
%>  This is the base class for generating objects
%>  that can display the time spinner on the console.<br>
%>
%>  \devnote
%>  The ``handle`` superclass is essential to allow
%>  object modification by the object methods.<br>
%>
%>  \note
%>  See the documentation of the class constructor.<br>
%>
%>  \note
%>  See below for information on the attributes (properties).<br>
%>
%>  \note
%>  See below for information on the methods.<br>
%>
%>  \final
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 5:39 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef Spinner < pm.matlab.Handle

    properties(Access = protected)
        %>
        %>  ``tickmarks``
        %>
        %>  The MATLAB ``char`` vector containing the set of
        %>  characters that represent the passage of time in the spinner.
        %>
        tickmarks = '|/-\';
    end

    properties(Hidden)
        %>
        %>  ``tickCount``
        %>
        %>  The MATLAB integer containing the length of ``tickmarks``.
        %>
        %>  \warning
        %>  This is an internal ``Hidden`` class attribute
        %>  that is inaccessible to the end users.<br>
        %>
        tickCount = 4;
        %>
        %>  ``format``
        %>
        %>  The MATLAB ``char`` vector containing the spinner display format.
        %>
        %>  \warning
        %>  This is an internal ``Hidden`` class attribute
        %>  that is inaccessible to the end users.<br>
        %>
        format = [repmat('\b', 1, 4 + 1), '%s'];
        %>
        %>  ``clock``
        %>
        %>  The scalar integer representing the index of the ``tickmarks`` attribute.
        %>
        %>  \warning
        %>  This is an internal ``Hidden`` class attribute
        %>  that is inaccessible to the end users.<br>
        %>
        clock = 0;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %>  \brief
        %>  Return a scalar object of class [pm.timing.Spinner](@ref Spinner).<br>
        %>
        %>  \details
        %>  This is the constructor of the class [pm.timing.Spinner](@ref Spinner).<br>
        %>
        %>  \param[in]  tickmarks   :   The input MATLAB ``char`` vector containing the set of
        %>                              characters that represent the passage of time in the spinner.<br>
        %>                              (**optional**, default = ``'|/-\'``)
        %>
        %>  \return
        %>  ``self``                :   The output scalar object of class [pm.timing.Spinner](@ref Spinner).<br>
        %>
        %>  \interface{Spinner}
        %>  \code{.m}
        %>
        %>      self = pm.timing.Spinner()
        %>      self = pm.timing.Spinner(tickmarks)
        %>
        %>  \endcode
        %>
        %>  \example{Spinner}
        %>  \include{lineno} example/timing/Spinner/main.m
        %>  \output{Spinner}
        %>  \include{lineno} example/timing/Spinner/main.out.m
        %>
        %>  \final{Spinner}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 5:43 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Spinner(tickmarks)
            if  nargin < 1
                tickmarks = '|/-\';
            end
            self.tickmarks = tickmarks;
            self.tickCount = length(self.tickmarks);
            self.format = [repmat('\b', 1, self.tickCount + 3), '%s'];
            self.clock = 0;
        end

        %>  \brief
        %>  Rotate the tick mark of the spinner
        %>  and display the percentage value of the input fraction.<br>
        %>  This is a dynamic method of the class [pm.timing.Spinner](@ref Spinner).<br>
        %>
        %>  \param[inout]   self        :   The **implicitly-passed** input argument representing the parent object of the method.<br>
        %>  \param[in]      fraction    :   The input scalar MATLAB fractional real (``0 <= fraction <= 1``)
        %>                                  representing the fraction of work so far accomplished.<br>
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
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
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