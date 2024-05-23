%>  \brief
%>  This is the abstract class for generating instances of objects
%>  that contain the specifications of various types of figures.<br>
%>  This is a generic class for generating figures
%>  containing a single subplot (axes).
classdef Plot < pm.vis.figure.Figure

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = public)
        %>
        %>  \param[in]  subplot :   The scalar object of superclass ``pm.vis.subplot.Subplot``
        %>                          representing the set of subplots to display in the figure.
        %>  
        subplot = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>
        %>  \param[in]  subplot     :   The input scalar object of superclass ``pm.vis.subplot.Subplot``.
        %>                              The input ``subplot`` object must minimally have the ``make()`` and ``reset()`` methods.
        %>  
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.
        %>
        %>  \note
        %>  The input ``varargin`` can also contain the components
        %>  of the ``subplot`` component of the parent object.
        %>
        %>  \return
        %>  `self`                  :   The output scalar object of class ``pm.vis.plot.Plot``.
        %>
        %>  \interface{Plot}
        %>  \code{.m}
        %>
        %>      plot = pm.vis.plot.Plot(subplot);
        %>
        %>  \endcode
        %>
        %>  \note
        %>  See the list of class attributes below,
        %>  also those of the superclass ``pm.vis.figure.Figure``.
        %>
        %>  \final{Plot}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 9:32 AM, University of Texas at Arlington<br>
        function self = Plot(subplot, varargin)
            varargin = {"subplot", subplot, varargin{:}};
            self = self@pm.vis.figure.Figure(varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Reset the properties of the plot to the original default settings.
        %>  Use this method when you change many attributes of the plot and
        %>  you want to clean up and go back to the default settings.
        %>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.
        %>
        %>  \note
        %>  The input ``varargin`` can also contain the components
        %>  of the ``subplot`` component of the parent object.
        %>
        %>  \return
        %>  `None`
        %>
        %>  \interface{reset}
        %>  \code{.m}
        %>
        %>      pm.vis.plot.Plot.reset() # reset the plot to the default settings.
        %>
        %>  \endcode
        %>
        %>  \final{reset}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 9:35 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function reset(self, varargin)

            [varfig, varsub] = pm.matlab.hashmap.popKeyVal(["figure", "subplot"], varargin);
            reset@pm.vis.figure.Figure(self, varfig{:});
            self.subplot.reset(varsub{:});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: Any non-MATLAB-default setting must be preferably set in premake() method to override user null values.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Configure the plot settings and specifications,
        %>  make the plot, and return nothing.
        %>
        %>  \details
        %>  In making the figure, this method we call the ``make()``
        %>  methods of each of the subplot objects stored in the
        %>  ``subplot`` component.
        %>
        %>  \warning
        %>  This method has side-effects by manipulating
        %>  the existing attributes of the parent object.
        %>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.
        %>
        %>  \note
        %>  The input ``varargin`` can also contain the components
        %>  of the ``subplot`` component of the parent object.
        %>
        %>  \return
        %>  `None`
        %>
        %>  \interface{make}
        %>  \code{.m}
        %>
        %>      p = pm.vis.plot.Plot.make(varargin);
        %>
        %>  \endcode
        %>
        %>  \example{make}
        %>
        %>      p = pm.vis.plot.Plot(pm.vis.subplot.Line());
        %>      p.make()
        %>
        %>  \final{make}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 9:37 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function make(self, varargin)

            [varfig, varsub] = pm.matlab.hashmap.popKeyVal(["figure", "subplot"], varargin);
            make@pm.vis.figure.Figure(self, varfig{:});
            self.subplot.make(varsub{:});
            hold off;

        end % function

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Configure the plot settings and specifications and return nothing.
        %>
        %>  \warning
        %>  This method has side-effects by manipulating
        %>  the existing attributes of the parent object.
        %>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``premake()`` method.
        %>
        %>  \return
        %>  `None`
        %>
        %>  \interface{premake}
        %>  \code{.m}
        %>
        %>      f = pm.vis.plot.Plot.premake(varargin);
        %>
        %>  \endcode
        %>
        %>  \example{getBorder}
        %>
        %>      f = pm.vis.plot.Plot(pm.vis.Line());
        %>      f.premake("figure", {"color", "none"})
        %>
        %>  \final{premake}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 9:39 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function premake(self, varargin)

            [varfig, varsub] = pm.matlab.hashmap.popKeyVal(["figure", "subplot"], varargin);
            premake@pm.vis.figure.Figure(self, varfig{:});
            if ~isempty(varsub)
                self.subplot.hash2comp(varsub);
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% These settings must happen here so that they can be reset every time user nullifies the values.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end