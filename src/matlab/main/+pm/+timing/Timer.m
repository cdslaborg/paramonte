%
%   This is the base class for generating objects
%   that can time interval consecutively.
%
%   The main utility of this timer class
%   is its dynamic ``del()`` method which
%   can compute the time elapsed since the
%   last measurement in one function call.
%
%   \devnote
%
%       The ``handle`` superclass is essential to allow
%       object modification by the object methods.
%
%       See the documentation of the class constructor.
%
%   \note
%       See below for information on the attributes (properties).
%
%   \note
%       See below for information on the methods.
%
%>  \return
%       An object of class ``pm.timing.Timer``.
%>
%>  \interface{}
%>  \code{.m}
%>  \endcode
%>
%       timer = pm.timing.Timer()
%
%   Example
%   -------
%
%       timer = pm.timing.Timer()
%       timer.tic();
%       pause(1);
%       timer.toc()
%       pause(.5);
%       timer.toc()
%       pause(.2);
%       timer.del()
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef Timer < pm.matlab.Handle

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = protected)
        %
        %   clock
        %
        %       The scalar MATLAB real containing the most recent
        %       timing since the construction of the timer.
        %
        clock = 0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Hidden)
        %
        %   start
        %
        %       The scalar MATLAB real containing the timer start.
        %
        start;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %
        %   Return a scalar object of class ``pm.timing.Timer``.
        %
        %   This is the constructor of the class ``pm.timing.Timer``.
        %
        %       None
        %>
        %>  \return
        %>
        %       self
        %
        %           The output scalar object of class ``pm.timing.Timer``.
        %
        %>  \interface{copy}
        %>  \code{.m}
        %>  \endcode
        %
        %       self = pm.timing.Timer()
        %
        %>  \final{}
        %>
        %>  \author
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Timer()
            self.tic();
        end

        %
        %   Reset the timer, equivalent to reconstructing the timer object.
        %
        %   This is a dynamic method of the class ``pm.timing.Timer()``.
        %
        %       None
        %>
        %>  \return
        %>
        %       None
        %
        %>  \interface{copy}
        %>  \code{.m}
        %>  \endcode
        %
        %       timer = pm.timing.Timer()
        %       timer.tic()
        %
        %>  \final{}
        %>
        %>  \author
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function tic(self)
            self.start = tic();
            self.clock = 0;
        end

        %
        %   Return a scalar MATLAB ``real`` containing the time
        %   past since the (re)construction of the timer object.
        %   Also, set the ``clock`` attribute of the parent object.
        %
        %   This is a dynamic method of the class ``pm.timing.Timer()``.
        %
        %       None
        %>
        %>  \return
        %>
        %       clock
        %
        %           The output scalar MATLAB ``real`` containing the time
        %           past since the (re)construction of the timer object.
        %
        %>  \interface{copy}
        %>  \code{.m}
        %>  \endcode
        %
        %       timer = pm.timing.Timer()
        %       clock = timer.toc()
        %
        %>  \final{}
        %>
        %>  \author
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function clock = toc(self)
            self.clock = toc(self.start);
            clock = self.clock;
        end

        %
        %   Return a scalar MATLAB ``real`` containing the time
        %   past since the last time measurement by the timer object.
        %   Also, set the ``clock`` attribute of the parent object.
        %
        %   This is a dynamic method of the class ``pm.timing.Timer()``.
        %
        %       None
        %>
        %>  \return
        %>
        %       delta
        %
        %           The output scalar MATLAB ``real`` containing the time
        %           past since the last time measurement by the timer object.
        %
        %>  \interface{copy}
        %>  \code{.m}
        %>  \endcode
        %
        %       timer = pm.timing.Timer()
        %       delta = timer.del()
        %
        %>  \final{}
        %>
        %>  \author
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function delta = del(self)
            delta = self.clock;
            delta = self.toc() - delta;
        end

    end

end