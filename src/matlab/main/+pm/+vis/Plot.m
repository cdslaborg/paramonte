%>  \brief
%>  This is the base class for generating instances of objects
%>  that contain the specifications of various types of figures.<br>
%>
%>  \details
%>  This is a generic class for generating
%>  figures containing a **single** subplot (axes).<br>
%>
%>  \note
%>  See the documentation of the constructor of the class
%>  [pm.vis.Plot](@ref Plot::Plot) for example usage.<br>
%>
%>  \note
%>  See the documentation of the attributes
%>  of the superclass [pm.vis.figure.Figure](@ref Figure).<br>
%>
%>  \see
%>  [pm.vis.Cascade](@ref Cascade)<br>
%>  [pm.vis.Subplot](@ref Subplot)<br>
%>  [pm.vis.Figure](@ref Figure)<br>
%>  [pm.vis.Corner](@ref Corner)<br>
%>  [pm.vis.Plot](@ref Plot)<br>
%>  [pm.vis.Tile](@ref Tile)<br>
%>
%>  \final
%>
%>  \AmirShahmoradi, July 5 2024, 1:07 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
classdef Plot < pm.vis.figure.Figure

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = public)
        %>
        %>  ``subplot``
        %>
        %>  The scalar object of superclass [pm.vis.Subplot](@ref Subplot)
        %>  representing the set of subplots to display in the figure.<br>
        %>
        subplot = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>
        %>  \param[in]  subplot     :   The input scalar object of superclass [pm.vis.Subplot](@ref Subplot).<br>
        %>                              The input ``subplot`` object must minimally have the ``make()`` and ``reset()`` methods.<br>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.<br>
        %>                              The input ``varargin`` can also contain the components
        %>                              of the ``subplot`` component of the parent object.<br>
        %>
        %>  \return
        %>  ``self``                :   The output scalar object of class [pm.vis.Plot](@ref Plot).<br>
        %>
        %>  \interface{Plot}
        %>  \code{.m}
        %>
        %>      plot = pm.vis.Plot(subplot);
        %>      plot = pm.vis.Plot(subplot, varargin);
        %>
        %>  \endcode
        %>
        %>  \note
        %>  The input ``varargin`` can also contain the components
        %>  of the ``subplot`` component of the parent object.<br>
        %>
        %>  \note
        %>  See the list of class attributes below,
        %>  also those of the superclass [pm.vis.figure.Figure](@ref Figure).<br>
        %>
        %>  \example{Plot}
        %>  \include{lineno} example/vis/Plot/main.m
        %>  \vis{Plot}
        %>  <br>\image html example/vis/Plot/PlotLine.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotLine3.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotScatter.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotScatter3.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotLineScatter.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotLineScatter3.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotHistfit.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotHistogram.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotHistogram2.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotContour.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotContourf.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotContour3.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotHeatmap.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotEllipse3.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotEllipse.1.png width=700
        %>
        %>  \final{Plot}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 9:32 AM, University of Texas at Arlington<br>
        %>  \AmirShahmoradi, July 5 2024, 1:07 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        function self = Plot(subplot, varargin)
            varargin = {"subplot", subplot, varargin{:}};
            self = self@pm.vis.figure.Figure(varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Reset the properties of the plot to the original default settings.<br>
        %>
        %>  \details
        %>  Use this method when you change many attributes of the plot and
        %>  you want to clean up and go back to the default settings.<br>
        %>
        %>  \param[inout]   self        :   The input/output parent object of class [pm.vis.Subplot](@ref Subplot)
        %>                                  which is **implicitly** passed to this dynamic method (not by the user).<br>
        %>  \param[in]      varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                                  If the property is a ``struct()``, then its value must be given as a cell array,
        %>                                  with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                                  Note that all of these property-value pairs can be also directly set via the
        %>                                  parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \interface{reset}
        %>  \code{.m}
        %>
        %>      pm.vis.Plot.reset() % reset the plot to the default settings.
        %>
        %>  \endcode
        %>
        %>  \note
        %>  The input ``varargin`` can also contain the components
        %>  of the ``subplot`` component of the parent object.<br>
        %>
        %>  \final{reset}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 9:35 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
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
        %>  make the plot, and return nothing.<br>
        %>
        %>  \details
        %>  In making the figure, this method we call the ``make()`` methods of
        %>  each of the subplot objects stored in the ``subplot`` component.<br>
        %>
        %>  \warning
        %>  This method has side-effects by manipulating
        %>  the existing attributes of the parent object.<br>
        %>
        %>  \param[inout]   self        :   The input/output parent object of class [pm.vis.Subplot](@ref Subplot)
        %>                                  which is **implicitly** passed to this dynamic method (not by the user).<br>
        %>  \param[in]      varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                                  If the property is a ``struct()``, then its value must be given as a cell array,
        %>                                  with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                                  Note that all of these property-value pairs can be also directly set via the
        %>                                  parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \interface{make}
        %>  \code{.m}
        %>
        %>      p = pm.vis.Plot.make(varargin);
        %>
        %>  \endcode
        %>
        %>  \note
        %>  The input ``varargin`` can also contain the components
        %>  of the ``subplot`` component of the parent object.<br>
        %>
        %>  \example{Plot}
        %>  \include{lineno} example/vis/Plot/main.m
        %>  \vis{Plot}
        %>  <br>\image html example/vis/Plot/PlotLine.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotLine3.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotScatter.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotScatter3.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotLineScatter.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotLineScatter3.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotHistfit.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotHistogram.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotHistogram2.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotContour.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotContourf.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotContour3.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotHeatmap.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotEllipse3.1.png width=700
        %>  <br>\image html example/vis/Plot/PlotEllipse.1.png width=700
        %>
        %>  \final{make}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 9:37 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function make(self, varargin)

            [varfig, varsub] = pm.matlab.hashmap.popKeyVal(["figure", "subplot"], varargin);
            make@pm.vis.figure.Figure(self, varfig{:});
            self.subplot.make(varsub{:});
            if  isprop(self.subplot, "axes") % Heatmap does not accept ``hold``.
                hold off;
            end

        end % function

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Configure the plot settings and specifications and return nothing.<br>
        %>
        %>  \warning
        %>  This method has side-effects by manipulating
        %>  the existing attributes of the parent object.<br>
        %>
        %>  \param[inout]   self        :   The input/output parent object of class [pm.vis.Subplot](@ref Subplot)
        %>                                  which is **implicitly** passed to this dynamic method (not by the user).<br>
        %>  \param[in]      varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                                  If the property is a ``struct()``, then its value must be given as a cell array,
        %>                                  with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                                  Note that all of these property-value pairs can be also directly set via the
        %>                                  parent object attributes, before calling the ``premake()`` method.<br>
        %>
        %>  \interface{premake}
        %>  \code{.m}
        %>
        %>      f = pm.vis.Plot.premake(varargin);
        %>
        %>  \endcode
        %>
        %>  \example{getBorder}
        %>
        %>      f = pm.vis.Plot(pm.vis.Line());
        %>      f.premake("figure", {"color", "none"})
        %>
        %>  \final{premake}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 9:39 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
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