%>  \brief
%>  This is the base class for generating objects
%>  that can time interval consecutively.<br>
%>
%>  \details
%>  The main utility of this timer class
%>  is its dynamic ``del()`` method which
%>  can compute the time elapsed since the
%>  last measurement in one function call.
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
%>  \JoshuaOsborne, May 21 2024, 5:47 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef Timer < pm.matlab.Handle

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = protected)
        %>
        %>  ``clock``
        %>
        %>  The scalar MATLAB real containing the most recent
        %>  timing since the construction of the timer.<br>
        %>
        %>  \warning
        %>  This is an internal ``protected`` class attribute
        %>  that is inaccessible to the end users.<br>
        %>
        clock = 0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Hidden)
        %>
        %>  ``start``
        %>
        %>  The scalar MATLAB real containing the timer start.<br>
        %>
        %>  \warning
        %>  This is an internal ``Hidden`` class attribute
        %>  that is inaccessible to the end users.<br>
        %>
        start;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %>  \brief
        %>  Return a scalar object of class [pm.timing.Timer](@ref Timer).<br>
        %>
        %>  \details
        %>  This is the constructor of the class [pm.timing.Timer](@ref Timer).<br>
        %>
        %>  \return
        %>  ``self``    :   The output scalar object of class [pm.timing.Timer](@ref Timer).<br>
        %>
        %>  \interface{Timer}
        %>  \code{.m}
        %>
        %>      self = pm.timing.Timer()
        %>
        %>  \endcode
        %>
        %>  \example{Timer}
        %>  \include{lineno} example/timing/Timer/main.m
        %>  \output{Timer}
        %>  \include{lineno} example/timing/Timer/main.out.m
        %>
        %>  \final{Timer}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 5:49 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Timer()
            self.tic();
        end

        %>  \brief
        %>  Reset the timer, equivalent to reconstructing the timer object.<br>
        %>  This is a dynamic method of the class [pm.timing.Timer](@ref Timer).
        %>
        %>  \interface{tic}
        %>  \code{.m}
        %>
        %>      timer = pm.timing.Timer()
        %>      timer.tic()
        %>
        %>  \endcode
        %>
        %>  \final{tic}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 5:50 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function tic(self)
            self.start = tic();
            self.clock = 0;
        end

        %>  \brief
        %>  Return a scalar MATLAB ``real`` containing the time
        %>  past since the (re)construction of the timer object.<br>
        %>
        %>  \details
        %>  Also, set the ``clock`` attribute of the parent object.<br>
        %>  This is a dynamic method of the class [pm.timing.Timer](@ref Timer).<br>
        %>
        %>  \return
        %>  ``clock``   :   The output scalar MATLAB ``real`` containing the time
        %>                  past since the (re)construction of the timer object.<br>
        %>
        %>  \interface{toc}
        %>  \code{.m}
        %>
        %>      timer = pm.timing.Timer()
        %>      clock = timer.toc()
        %>
        %>  \endcode
        %>
        %>  \final{toc}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 5:58 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function clock = toc(self)
            self.clock = toc(self.start);
            clock = self.clock;
        end

        %>  \brief
        %>  Return a scalar MATLAB ``real`` containing the time past
        %>  since the last time measurement by the timer object.<br>
        %>  Also, set the ``clock`` attribute of the parent object.<br>
        %>  This is a dynamic method of the class [pm.timing.Timer](@ref Timer).<br>
        %>
        %>  \return
        %>  ``delta``   :   The output scalar MATLAB ``real`` containing the time
        %>                  past since the last time measurement by the timer object.<br>
        %>
        %>  \interface{del}
        %>  \code{.m}
        %>
        %>      timer = pm.timing.Timer()
        %>      delta = timer.del()
        %>  \endcode
        %>
        %>  \final{del}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 5:59 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function delta = del(self)
            delta = self.clock;
            delta = self.toc() - delta;
        end

    end

end