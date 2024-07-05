%>  \brief
%>  This is the class for generating instances of objects
%>  that contain the specifications of various types of plots.<br>
%>
%>  \details
%>  This class primarily serves as the superclass for
%>  the visualization-ready subclass [pm.vis.subplot.Subplot](@ref Subplot)
%>  and its subclasses, all accessible to the end users.<br>
%>
%>  Dynamic class attributes
%>  ------------------------
%>
%>  This class contains a set of attributes that are defined dynamically at runtime
%>  for the output object depending on its subclass (plot type it represents).<br>
%>  The following is the list of all class attributes that are dynamically added
%>  to the instantiated class objects based on the specified input plot type.<br>
%>  See also the explicit class and superclass attributes not listed below.<br>
%>
%>  <ol>
%>      <li>    ``axes`` (available for all subplots except [pm.vis.subplot.Heatmap](@ref Heatmap))<br>
%>
%>              A MATLAB ``struct`` whose fields and values are passed as
%>              keyword arguments to the MATLAB intrinsic ``set()`` for
%>              the current active axes object in the plot ``gca()``.<br>
%>
%>      <li>    ``colorbar`` (available for all axes types that allow color-mapping)<br>
%>
%>              A MATLAB ``struct`` whose fields and their values will
%>              be passed as keyword arguments to the MATLAB intrinsic ``colorbar``.<br>
%>              The following are the default components of ``colorbar``:<br>
%>
%>              <ol>
%>                  <li>    ``enabled``
%>
%>                          A ``logical`` value. If ``true``, the
%>                          color bar will be applied to the axes.<br>
%>
%>                  <li>    others
%>
%>                          See the acceptable keyword arguments
%>                          of the MATLAB intrinsic ``colorbar()``.<br>
%>              </ol>
%>
%>              \warning
%>              Keep in mind that MATLAB keyword arguments are case-INsensitive.<br>
%>              Hence, ensure you do not add the same keyword as multiple different fields.<br>
%>              For example, ``colorbar.color`` and ``colorbar.Color`` are the same,
%>              and only one of the two will be processed.<br>
%>
%>              \example{colorbar}
%>              \code{.m}
%>
%>                  self.colorbar.enabled = true;
%>                  self.colorbar.location = "west";
%>
%>              \endcode
%>              <br>
%>
%>      <li>    ``colormap`` (available for all axes types that allow color-mapping)
%>
%>              A MATLAB ``struct`` whose fields and their values will be passed
%>              as keyword arguments to the MATLAB intrinsic ``colormap``.<br>
%>              The following are the default components of ``colormap``:<br>
%>
%>              <ol>
%>                  <li>    ``enabled``
%>
%>                          A logical value. If ``true``, the
%>                          colormap will be applied to the axes.<br>
%>
%>                  <li>    ``map``
%>
%>                          A string or a vector of color triplets or any other value
%>                          that the intrinsic MATLAB ``colormap`` accepts as input.<br>
%>
%>                          This option is relevant only to visualizations that allow color-mapping.<br>
%>              </ol>
%>
%>              \warning
%>              Keep in mind that MATLAB keyword arguments are case-INsensitive.<br>
%>              Hence, ensure you do not add the same keyword as multiple different fields.<br>
%>              For example, ``colormap.map`` and ``colormap.Map`` are the same,
%>              and only one of the two will be processed.<br>
%>
%>              \example{colormap}
%>              \code{.m}
%>
%>                  self.colormap.enabled = true;
%>                  self.colormap.map = "winter";
%>                  self.colormap.map = "winter";
%>                  self.colormap.map = 'default';
%>                  self.colormap.map = pm.vis.cmap.redblue();
%>
%>              \endcode
%>              <br>
%>
%>      <li>    ``contour`` (available only for [pm.vis.subplot.Contour](@ref Contour) axes types)
%>
%>              A MATLAB ``struct`` whose fields and their values will be passed
%>              as keyword arguments to the MATLAB intrinsic ``contour``.<br>
%>              The following are the default components of ``contour``:<br>
%>
%>              <ol>
%>                  <li>    ``enabled``
%>
%>                          A logical value. If ``true``, the
%>                          contour will be added to the axes.
%>
%>                  <li>    ``levels``
%>
%>                          See the corresponding positional argument of the MATLAB intrinsic ``contourf()``.<br>
%>
%>                  <li>    ``lineSpec``
%>
%>                          See the corresponding positional argument of the MATLAB intrinsic ``contourf()``.<br>
%>
%>                  <li>    others
%>
%>                          See the acceptable keyword arguments of the MATLAB intrinsic ``contour()``.<br>
%>              </ol>
%>
%>              \warning
%>              Keep in mind that MATLAB keyword arguments are case-INsensitive.<br>
%>              Hence, ensure you do not add the same keyword as multiple different fields.<br>
%>
%>              \example{colormap}
%>              \code{.m}
%>
%>                  self.contour.enabled = true;
%>                  self.contour.lineWidth = "none";
%>
%>              \endcode
%>              <br>
%>
%>      <li>    ``contour3`` (available only for [pm.vis.subplot.Contour3](@ref Contour3) axes types)
%>
%>              A MATLAB ``struct`` whose fields and their values will be passed
%>              as keyword arguments to the MATLAB intrinsic ``contour3``.<br>
%>              The following are the default components of ``contour3``:<br>
%>
%>              <ol>
%>                  <li>    ``enabled``
%>
%>                          A logical value. If ``true``, the
%>                          contour3 will be added to the axes.<br>
%>
%>                  <li>    ``levels``
%>
%>                          See the corresponding positional argument of the MATLAB intrinsic ``contourf()``.<br>
%>
%>                  <li>    ``lineSpec``
%>
%>                          See the corresponding positional argument of the MATLAB intrinsic ``contourf()``.<br>
%>
%>                  <li>    others
%>
%>                          See the acceptable keyword arguments of the MATLAB intrinsic ``contour3()``.<br>
%>              </ol>
%>
%>              \warning
%>              Keep in mind that MATLAB keyword arguments are case-INsensitive.<br>
%>              Hence, ensure you do not add the same keyword as multiple different fields.<br>
%>
%>              \example{colormap}
%>              \code{.m}
%>
%>                  self.contour3.enabled = true;
%>                  self.contour3.lineWidth = "none";
%>
%>              \endcode
%>              <br>
%>
%>      <li>    ``contourf`` (available only for [pm.vis.subplot.Contourf](@ref Contourf) axes types)
%>
%>              A MATLAB ``struct`` whose fields and their values will be passed
%>              as keyword arguments to the MATLAB intrinsic ``contourf``.<br>
%>              The following are the default components of ``contourf``:<br>
%>
%>              <ol>
%>                  <li>    ``enabled``
%>
%>                          A logical value. If ``true``, the
%>                          contourf will be added to the axes.<br>
%>
%>                  <li>    ``levels``
%>
%>                          See the corresponding positional argument of the MATLAB intrinsic ``contourf()``.<br>
%>
%>                  <li>    ``lineSpec``
%>
%>                          See the corresponding positional argument of the MATLAB intrinsic ``contourf()``.<br>
%>
%>                  <li>    others
%>
%>                          See the acceptable keyword arguments of the MATLAB intrinsic ``contourf()``.<br>
%>              </ol>
%>
%>              \warning
%>              Keep in mind that MATLAB keyword arguments are case-INsensitive.<br>
%>              Hence, ensure you do not add the same keyword as multiple different fields.<br>
%>
%>              \example{contourf}
%>              \code{.m}
%>
%>                  self.contourf.enabled = true;
%>                  self.contourf.lineWidth = "none";
%>
%>              \endcode
%>              <br>
%>
%>      <li>    ``histfit`` (available only for [pm.vis.subplot.Histfit](@ref Histfit) axes types)
%>
%>              A MATLAB ``struct`` whose fields and their values will be passed
%>              as keyword arguments to the MATLAB intrinsic ``histfit``.<br>
%>              The following are the default components of ``histfit``:<br>
%>
%>              <ol>
%>                  <li>    ``enabled``
%>
%>                          A logical value. If ``true``, the
%>                          histfit will be added to the axes.<br>
%>
%>                  <li>    ``nbins``
%>
%>                          See the corresponding positional argument of the MATLAB intrinsic ``histfit()``.<br>
%>
%>                  <li>    ``dist``
%>
%>                          See the corresponding positional argument of the MATLAB intrinsic ``histfit()``.<br>
%>              </ol>
%>
%>              \warning
%>              Keep in mind that MATLAB keyword arguments are case-INsensitive.<br>
%>              Hence, ensure you do not add the same keyword as multiple different fields.<br>
%>
%>              \example{histfit}
%>              \code{.m}
%>
%>                  self.histfit.enabled = true;
%>                  self.histfit.nbins = 20;
%>
%>              \endcode
%>              <br>
%>
%>      <li>    ``histogram`` (available only for [pm.vis.subplot.Histogram](@ref Histogram) axes types)
%>
%>              A MATLAB ``struct`` whose fields and their values will be passed
%>              as keyword arguments to the MATLAB intrinsic ``histogram``.<br>
%>              The following are the default components of ``histogram``:<br>
%>
%>              <ol>
%>                  <li>    ``enabled``
%>
%>                          A logical value. If ``true``, the
%>                          histogram will be added to the axes.<br>
%>
%>                  <li>    ``nbins``
%>
%>                          See the corresponding positional argument of the MATLAB intrinsic ``histogram()``.<br>
%>
%>                  <li>    ``edges``
%>
%>                          See the corresponding positional argument of the MATLAB intrinsic ``histogram()``.<br>
%>
%>                  <li>    others
%>
%>                          See the corresponding acceptable keyword arguments of the MATLAB intrinsic ``histogram()``.<br>
%>              </ol>
%>
%>              \warning
%>              Keep in mind that MATLAB keyword arguments are case-INsensitive.<br>
%>              Hence, ensure you do not add the same keyword as multiple different fields.<br>
%>
%>              \example{histogram}
%>              \code{.m}
%>
%>                  self.histogram.enabled = true;
%>                  self.histogram.edgeColor = "none";
%>                  self.histogram.nbins = 20;
%>
%>              \endcode
%>              <br>
%>
%>      <li>    ``histogram2`` (available only for [pm.vis.subplot.Histogram2](@ref Histogram2) axes types)
%>
%>              A MATLAB ``struct`` whose fields and their values will be passed
%>              as keyword arguments to the MATLAB intrinsic ``histogram2``.<br>
%>              The following are the default components of ``histogram2``:<br>
%>
%>              <ol>
%>                  <li>    ``enabled``
%>
%>                          A logical value. If ``true``, the
%>                          histogram2 will be added to the axes.<br>
%>
%>                  <li>    ``nbins``
%>
%>                          See the corresponding positional argument of the MATLAB intrinsic ``histogram2()``.<br>
%>
%>                  <li>    ``xedges``
%>
%>                          See the corresponding positional argument of the MATLAB intrinsic ``histogram2()``.<br>
%>
%>                  <li>    ``yedges``
%>
%>                          See the corresponding positional argument of the MATLAB intrinsic ``histogram2()``.<br>
%>
%>                  <li>    others
%>
%>                          See the corresponding acceptable keyword arguments of the MATLAB intrinsic ``histogram2()``.<br>
%>              </ol>
%>
%>              \warning
%>              Keep in mind that MATLAB keyword arguments are case-INsensitive.<br>
%>              Hence, ensure you do not add the same keyword as multiple different fields.<br>
%>
%>              \example{histogram2}
%>              \code{.m}
%>
%>                  self.histogram2.enabled = true;
%>                  self.histogram2.edgeColor = "none";
%>                  self.histogram2.nbins = 20;
%>
%>              \endcode
%>              <br>
%>
%>      <li>    ``legend`` (available for all axes types except [pm.vis.subplot.Heatmap](@ref Heatmap))
%>
%>              A MATLAB ``struct`` whose fields and values are passed
%>              as keyword arguments to the MATLAB intrinsic ``title()``.<br>
%>
%>      <li>    ``maxnoise`` (available only for [Contour](@ref Contour)/[Contourf](@ref Contourf)/[Contour3](@ref Contour3) axes types)
%>
%>              A float indicating the threshold below which the kernel density
%>              estimate is considered to be noise and is rounded to zero.<br>
%>              The higher this value is, the less noise will be
%>              visible in the resulting contour plots.<br>
%>              If empty, the default value is ``0.001``.<br>
%>
%>      <li>    ``plot`` (available only for [pm.vis.subplot.Line](@ref Line)/[pm.vis.subplot.LineScatter](@ref LineScatter) axes types)
%>
%>              A MATLAB ``struct`` whose fields and their values will be passed
%>              as keyword arguments to the MATLAB intrinsic ``plot``.<br>
%>              The following are the default components of ``plot``:<br>
%>
%>              <ol>
%>                  <li>    ``enabled``
%>
%>                          A logical value. If ``true``, the
%>                          plot will be added to the axes.<br>
%>
%>                  <li>    ``lineSpec``
%>
%>                          See the corresponding positional argument of the MATLAB intrinsic ``plot()``.<br>
%>
%>                  <li>    others
%>
%>                          See the corresponding acceptable keyword arguments of the MATLAB intrinsic ``plot()``.<br>
%>              </ol>
%>
%>              \warning
%>              Keep in mind that MATLAB keyword arguments are case-INsensitive.<br>
%>              Hence, ensure you do not add the same keyword as multiple different fields.<br>
%>
%>              \example{plot}
%>              \code{.m}
%>
%>                  self.plot.enabled = true;
%>                  self.plot.lineWidth = 1;
%>
%>              \endcode
%>              <br>
%>
%>      <li>    ``plot3`` (available only for [Line3](@ref Line3)/[LineScatter3](@ref LineScatter3) axes types)
%>
%>              A MATLAB ``struct`` whose fields and their values will be passed
%>              as keyword arguments to the MATLAB intrinsic ``plot3``.<br>
%>              The following are the default components of ``plot3``:<br>
%>
%>              <ol>
%>                  <li>    ``enabled``
%>
%>                          A logical value. If ``true``, the
%>                          ``plot3`` will be added to the axes.<br>
%>
%>                  <li>    ``lineSpec``
%>
%>                          See the corresponding positional argument of the MATLAB intrinsic ``plot3()``.<br>
%>
%>                  <li>    others
%>
%>                          See the corresponding acceptable keyword arguments of the MATLAB intrinsic ``plot3()``.<br>
%>              </ol>
%>
%>              \warning
%>              Keep in mind that MATLAB keyword arguments are case-INsensitive.<br>
%>              Hence, ensure you do not add the same keyword as multiple different fields.<br>
%>
%>              \example{plot3}
%>              \code{.m}
%>
%>                  self.plot3.enabled = true;
%>                  self.plot3.lineWidth = 1;
%>
%>              \endcode
%>              <br>
%>
%>      <li>    ``precision`` (available only for [pm.vis.subplot.Heatmap](@ref Heatmap) axes types)
%>
%>              A scalar integer representing the number of digits after
%>              the decimal point for the values that appear in each cell
%>              of the heatmap. The default value is set by MATLAB.<br>
%>
%>      <li>    ``resolution`` (available only for [Contour](@ref Contour)/[Contourf](@ref Contourf)/[Contour3](@ref Contour3) axes types)
%>
%>              A scalar integer indicating the grid resolution for discretization of
%>              the data during the kernel density estimation. It must be a power of
%>              two, otherwise it will be changed to the next power of two at the
%>              time of using it. If empty, the default value is ``2^9``.<br>
%>
%>      <li>    ``scatter`` (available only for [Scatter](@ref Scatter)/[LineScatter](@ref LineScatter) axes types)
%>
%>              A MATLAB ``struct`` whose fields and their values will be passed
%>              as keyword arguments to the MATLAB intrinsic ``scatter``.<br>
%>              The following are the default components of ``scatter``:<br>
%>
%>              <ol>
%>                  <li>    ``enabled``
%>
%>                          A logical value. If ``true``, the
%>                          scatter will be added to the axes.<br>
%>
%>                  <li>    ``size``
%>
%>                          See the corresponding ``sz`` positional argument of the MATLAB intrinsic ``plot3()``.<br>
%>
%>                  <li>    ``color``
%>
%>                          See the corresponding ``C`` positional argument of the MATLAB intrinsic ``plot3()``.<br>
%>
%>                  <li>    ``filled``
%>
%>                          See the corresponding ``filled`` positional argument of the MATLAB intrinsic ``plot3()``.<br>
%>
%>                  <li>    ``marker``
%>
%>                          See the corresponding ``mkr`` positional argument of the MATLAB intrinsic ``plot3()``.<br>
%>
%>                  <li>    others
%>
%>                          See the corresponding acceptable keyword arguments of the MATLAB intrinsic ``scatter()``.<br>
%>              </ol>
%>
%>              \warning
%>              Keep in mind that MATLAB keyword arguments are case-INsensitive.<br>
%>              Hence, ensure you do not add the same keyword as multiple different fields.<br>
%>
%>              \example{scatter}
%>              \code{.m}
%>
%>                  self.scatter.enabled = true; % add scatter()
%>                  self.scatter.color = "red"; % set the points color
%>                  self.scatter.marker = "."; % set the marker type
%>                  self.scatter.size = 10; % set the point size
%>                  self.scatter.lineWidth = 0;
%>
%>              \endcode
%>              <br>
%>
%>      <li>    ``scatter3`` (available only for [Scatter3](@ref Scatter3)/[LineScatter3](@ref LineScatter3) axes types)
%>
%>              A MATLAB ``struct`` whose fields and their values will be passed
%>              as keyword arguments to the MATLAB intrinsic ``scatter3``.<br>
%>              The following are the default components of ``scatter3``:<br>
%>
%>              <ol>
%>                  <li>    ``enabled``
%>
%>                          A logical value. If ``true``, the
%>                          scatter3 will be added to the axes.<br>
%>
%>                  <li>    ``size``
%>
%>                          See the corresponding ``sz`` positional argument of the MATLAB intrinsic ``plot3()``.<br>
%>
%>                  <li>    ``color``
%>
%>                          See the corresponding ``C`` positional argument of the MATLAB intrinsic ``plot3()``.<br>
%>
%>                  <li>    ``filled``
%>
%>                          See the corresponding ``filled`` positional argument of the MATLAB intrinsic ``plot3()``.<br>
%>
%>                  <li>    ``marker``
%>
%>                          See the corresponding ``mkr`` positional argument of the MATLAB intrinsic ``plot3()``.<br>
%>
%>                  <li>    others
%>
%>                          See the corresponding acceptable keyword arguments of the MATLAB intrinsic ``scatter3()``.<br>
%>              </ol>
%>
%>              \warning
%>              Keep in mind that MATLAB keyword arguments are case-INsensitive.<br>
%>              Hence, ensure you do not add the same keyword as multiple different fields.<br>
%>
%>              \example{scatter3}
%>              \code{.m}
%>
%>                  self.scatter3.enabled = true; % add scatter3()
%>                  self.scatter3.color = "red"; % set the points color
%>                  self.scatter3.marker = "."; % set the marker type
%>                  self.scatter3.size = 10; % set the point size
%>                  self.scatter3.lineWidth = 0;
%>
%>              \endcode
%>              <br>
%>
%>      <li>    ``surface`` (available only for [Line](@ref Line)/[LineScatter](@ref LineScatter)/[Line3](@ref Line3)/[LineScatter3](@ref LineScatter3) axes types)
%>
%>              A MATLAB ``struct`` whose fields and their values will be passed
%>              as keyword arguments to the MATLAB intrinsic ``surface``.<br>
%>              The following are the default components of ``surface``:<br>
%>
%>              <ol>
%>                  <li>    ``enabled``
%>
%>                          A logical value. If ``true``, the
%>                          surface will be added to the axes.<br>
%>
%>                  <li>    others
%>
%>                          See the corresponding acceptable keyword arguments of the MATLAB intrinsic ``surface()``.<br>
%>              </ol>
%>
%>              \warning
%>              Keep in mind that MATLAB keyword arguments are case-INsensitive.<br>
%>              Hence, ensure you do not add the same keyword as multiple different fields.<br>
%>
%>              \example{surface}
%>              \code{.m}
%>
%>                  self.surface.enabled = true;
%>                  self.surface.edgeColor = 1;
%>
%>              \endcode
%>              <br>
%>
%>      <li>    ``title`` (available for all axes types)
%>
%>              A MATLAB ``struct`` whose fields and values are passed
%>              as keyword arguments to the MATLAB intrinsic ``title()``.<br>
%>
%>      <li>    ``xlabel`` (available for all axes types)
%>
%>              A MATLAB ``struct`` whose fields and values are passed
%>              as keyword arguments to the MATLAB intrinsic ``xlabel()``.<br>
%>
%>      <li>    ``ylabel`` (available for all axes types)
%>
%>              A MATLAB ``struct`` whose fields and values are passed
%>              as keyword arguments to the MATLAB intrinsic ``ylabel()``.<br>
%>
%>      <li>    ``zlabel`` (available only for all tri-axes axes types)
%>
%>              A MATLAB ``struct`` whose fields and values are passed
%>              as keyword arguments to the MATLAB intrinsic ``zlabel()``.<br>
%>
%>      <li>    ``xlim`` (available for all axes types)
%>
%>              A MATLAB vector of length ``2`` whose fields and values are
%>              passed as keyword arguments to the MATLAB intrinsic ``xlim()``.<br>
%>
%>      <li>    ``ylim`` (available for all axes types)
%>
%>              A MATLAB vector of length ``2`` whose fields and values are
%>              passed as keyword arguments to the MATLAB intrinsic ``ylim()``.<br>
%>
%>      <li>    ``zlim`` (available only for all tri-axes axes types)
%>
%>              A MATLAB vector of length ``2`` whose fields and values are
%>              passed as keyword arguments to the MATLAB intrinsic ``zlim()``.<br>
%>
%>      <li>    ``xscale`` (available for all axes types)
%>
%>              A MATLAB string whose value is passed directly to the MATLAB intrinsic
%>              ``xscale()`` to set the axis scale to either logarithmic or linear.<br>
%>              Possible values are: ``"log"``, ``"linear"``.<br>
%>              The default behavior is set by MATLAB.<br>
%>
%>      <li>    ``yscale`` (available for all axes types)
%>
%>              A MATLAB string whose value is passed directly to the MATLAB intrinsic
%>              ``yscale()`` to set the axis scale to either logarithmic or linear.<br>
%>              Possible values are: ``"log"``, ``"linear"``.<br>
%>              The default behavior is set by MATLAB.<br>
%>
%>      <li>    ``zscale`` (available only for all tri-axes axes types)
%>
%>              A MATLAB string whose value is passed directly to the MATLAB intrinsic
%>              ``zscale()`` to set the axis scale to either logarithmic or linear.<br>
%>              Possible values are: ``"log"``, ``"linear"``.<br>
%>              The default behavior is set by MATLAB.<br>
%>  </ol>
%>
%>  \devnote
%>  While dynamic addition of class attributes is not ideal, the current
%>  design was deemed unavoidable and best, given the constraints of the
%>  MATLAB language and visualization tools.<br>
%>
%>  \see
%>  [pm.matlab.Handle](@ref Handle)<br>
%>
%>  \final
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef Axes < pm.matlab.Handle

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = protected, Hidden)
        %>
        %>  ``type``        :   auxiliary struct containing plot type information.
        %>
        type = struct();
        %>
        %>  ``cenabled``    :   auxiliary logical scalar that is true if plot is color-mapped.
        %>
        cenabled = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Construct and return an object of class [pm.vis.axes.Axes](@ref Axes).<br>
        %>
        %>  \details
        %>  This function is the constructor of the class [pm.vis.axes.Axes](@ref Axes).<br>
        %>  For more information, see the documentation of the class [pm.vis.axes.Axes](@ref Axes).<br>
        %>
        %>  \param[in]  ptype       :   The input scalar MATLAB string containing the name of the
        %>                              subclass that whose parent is Axes (e.g., "heatmap").<br>
        %>                              Supported plot names are:<br>
        %>                              <ol>
        %>                                  <li>    ``"line"``
        %>                                  <li>    ``"line3"``
        %>                                  <li>    ``"scatter"``
        %>                                  <li>    ``"scatter3"``
        %>                                  <li>    ``"lineScatter"``
        %>                                  <li>    ``"lineScatter3"``
        %>                                  <li>    ``"histogram2"``
        %>                                  <li>    ``"histogram"``
        %>                                  <li>    ``"contour3"``
        %>                                  <li>    ``"contourf"``
        %>                                  <li>    ``"contour"``
        %>                                  <li>    ``"histfit"``
        %>                                  <li>    ``"heatmap"``
        %>                              </ol>
        %>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``premake()`` method.<br>
        %>
        %>  \return
        %>  ``self``                :   The output scalar object of class ``pm.vis.axes.Axes``.<br>
        %>
        %>  \interface{Axes}
        %>  \code{.m}
        %>
        %>      axes = pm.vis.axes.Axes(ptype);
        %>
        %>  \endcode
        %>
        %>  \final{Axes}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 5:59 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Axes(ptype, varargin)

            if  nargin < 1 || ~pm.introspection.istype(ptype, "string", 1) || ~pm.array.len(ptype)
                help("pm.vis.axes.Axes");
                error   ( newline ...
                        + "The input argument ``ptype`` is missing." + newline ...
                        + "For more information, see the class documentation displayed above." + newline ...
                        + newline ...
                        );
            else
                % lower the first character.
                ptype = string(ptype);
                ptype{1} = lower(ptype{1});
            end

            self.type = struct();
            self.type.name = ptype;
            pnamel = lower(ptype);
            self.type.is.line           = strcmpi(pnamel, "line"        );
            self.type.is.line3          = strcmpi(pnamel, "line3"       );
            self.type.is.scatter        = strcmpi(pnamel, "scatter"     );
            self.type.is.scatter3       = strcmpi(pnamel, "scatter3"    );
            self.type.is.lineScatter    = strcmpi(pnamel, "lineScatter" );
            self.type.is.lineScatter3   = strcmpi(pnamel, "lineScatter3");
            self.type.is.histogram2     = strcmpi(pnamel, "histogram2"  );
            self.type.is.histogram      = strcmpi(pnamel, "histogram"   );
            self.type.is.contour3       = strcmpi(pnamel, "contour3"    );
            self.type.is.contourf       = strcmpi(pnamel, "contourf"    );
            self.type.is.contour        = strcmpi(pnamel, "contour"     );
            self.type.is.heatmap        = strcmpi(pnamel, "heatmap"     );
            self.type.is.histfit        = strcmpi(pnamel, "histfit"     );
            self.type.is.diffusion      = self.type.is.contour || self.type.is.contourf || self.type.is.contour3;
            self.type.is.density        = self.type.is.diffusion || self.type.is.histfit || self.type.is.histogram || self.type.is.histogram2;
            self.type.is.d3             = self.type.is.line3 || self.type.is.scatter3 || self.type.is.lineScatter3;
            self.type.is.d1             = self.type.is.histfit || self.type.is.histogram || self.type.is.heatmap;
            self.type.is.d2             = ~(self.type.is.d1 || self.type.is.d3);
            self.type.is.triaxes        = self.type.is.d3 || self.type.is.histogram2 || self.type.is.contour3;
           %self.type.is.targetable     = false; %xx ~self.type.is.heatmap && ~ self.type.is.d3;
            self.type.has.line          = self.type.is.line || self.type.is.line3 || self.type.is.lineScatter || self.type.is.lineScatter3;
            self.type.has.scatter       = self.type.is.scatter || self.type.is.scatter3 || self.type.is.lineScatter || self.type.is.lineScatter3;

            self.reset(varargin{:}); % This is the subclass method!

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Reset the properties of the plot to the original default settings and return nothing.<br>
        %>
        %>  \details
        %>  Use this method when you change many attributes of the plot and
        %>  you want to clean up and go back to the default settings.<br>
        %>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes.<br>
        %>
        %>  \interface{reset}
        %>  \code{.m}
        %>
        %>      a = pm.vis.axes.Axes(ptype)
        %>
        %>      a.reset(varargin) % reset the plot to settings in ``varargin`` and the rest to default.
        %>      a.reset() % reset the plot to the default settings.
        %>
        %>  \endcode
        %>
        %>  \final{reset}
        %>
        %>  \author
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function reset(self, varargin)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: Any non-MATLAB-default setting must be preferably set in the premake() method to override user null values.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            fontSize_def = 12;

            %%%% axes

            if ~self.type.is.heatmap
                self.newprop("axes", struct());
                self.axes.box = [];
                self.axes.color = [];
                self.axes.colorScale = [];
                self.axes.fontName = [];
                self.axes.fontSize = fontSize_def;
                self.axes.fontSizeMode = [];
                self.axes.fontSmoothing = [];
                self.axes.fontWeight = [];
                self.axes.xgrid = [];
                self.axes.ygrid = [];
                if  self.type.is.triaxes
                    self.axes.zgrid = [];
                end
                self.axes.enabled = [];
            end

            %%%% title

            self.newprop("title", struct());
            if ~self.type.is.heatmap
                self.title.fontSize = fontSize_def;
                self.title.interpreter = [];
                self.title.fontWeight = [];
                self.title.color = []; %[0, 0, 0];
            end
            self.title.enabled = [];
            self.title.titletext = [];
            if ~self.type.is.heatmap
                self.title.subtitletext = [];
            end

            %%%% xlabel, xlim

            self.newprop("xlabel", struct());
            if ~self.type.is.heatmap
                self.xlabel.color = []; %[0.15, 0.15, 0.15];
                self.xlabel.fontAngle = [];
                self.xlabel.fontSize = fontSize_def;
                self.xlabel.fontWeight = [];
                self.xlabel.interpreter = [];
                self.xlabel.rotation = [];
            end
            self.xlabel.enabled = [];
            self.xlabel.txt = [];

            self.newprop("xlim", []);
            self.newprop("xscale", []);

            %%%% ylabel, ylim

            self.newprop("ylim", []);
            self.newprop("ylabel", self.xlabel);
            self.newprop("yscale", []);

            %%%% zlabel, zlim

            if  self.type.is.triaxes
                self.newprop("zlabel", self.xlabel);
                self.newprop("zlim", []);
                self.newprop("zscale", []);
            end

            %%%% colc, colorbar, colormap

            if  self.type.is.heatmap || ~self.type.is.d1

                if  self.type.has.line || self.type.has.scatter
                    self.newprop("colc", {});
                end

                if ~self.type.is.heatmap
                    self.newprop("colorbar", struct());
                    self.colorbar.fontSize = fontSize_def;
                    self.colorbar.direction = 'normal';
                    self.colorbar.limits = [];
                    self.colorbar.location = 'eastoutside';
                    self.colorbar.ticks = [];
                    self.colorbar.tickLabels = [];
                    self.colorbar.tickLabelInterpreter = 'tex';
                    self.colorbar.enabled = [];
                end

                self.newprop("colormap", struct());
                self.colormap.enabled = [];
                self.colormap.map = []; % 'default'

            end

            %%%% legend

            if ~self.type.is.heatmap
                self.newprop("legend", struct());
                self.legend.box = 'off';
                self.legend.color = 'none';
                self.legend.fontSize = fontSize_def;
                self.legend.interpreter = 'none';
                self.legend.location = [];
                self.legend.numColumns = [];
                self.legend.textColor = [];
                self.legend.enabled = [];
                self.legend.labels = {};
            end

            %%%% target

            %if ~self.type.is.targetable
            %    self.newprop("target");
            %end

            %%%%%%%%%%%%%%%%%%%%%%%
            %%%% heatmap attributes
            %%%%%%%%%%%%%%%%%%%%%%%

            if  self.type.is.heatmap
                self.newprop("heatmap", struct());
                self.newprop("precision", []);
                self.heatmap.colorbarVisible = [];
                self.heatmap.colorLimits = [];
                self.heatmap.missingDataColor = [];
                self.heatmap.fontName = [];
                self.heatmap.fontSize = fontSize_def;
                self.heatmap.enabled = [];
            end

            %%%%%%%%%%%%%%%%%%%%%%%
            %%%% density attributes
            %%%%%%%%%%%%%%%%%%%%%%%

            if  self.type.is.histfit
                self.newprop("histfit", struct());
                self.histfit.enabled = [];
                self.histfit.dist = [];
                self.histfit.nbins = [];
            end

            if  self.type.is.histogram
                self.newprop("histogram", struct());
                self.histogram.binLimitsMode = [];
                self.histogram.binMethod = [];
                self.histogram.binWidth = [];
                self.histogram.displayStyle = [];
                self.histogram.edgeAlpha = [];
                self.histogram.edgeColor = [];
                self.histogram.faceAlpha = [];
                self.histogram.faceColor = [];
                self.histogram.lineStyle = [];
                self.histogram.lineWidth = [];
                self.histogram.normalization = [];
                self.histogram.orientation = [];
                self.histogram.enabled = [];
                self.histogram.edges = [];
                self.histogram.nbins = [];
            end

            if  self.type.is.histogram2
                self.newprop("histogram2", struct());
                self.histogram2.binWidth = [];
                self.histogram2.xbinLimits = [];
                self.histogram2.xbinLimitsMode = [];
                self.histogram2.ybinLimits = [];
                self.histogram2.ybinLimitsMode = [];
                self.histogram2.binMethod = [];
                self.histogram2.showEmptyBins = [];
                self.histogram2.normalization = [];
                self.histogram2.displayStyle = [];
                self.histogram2.edgeAlpha = [];
                self.histogram2.edgeColor = [];
                self.histogram2.faceAlpha = [];
                self.histogram2.faceColor = [];
                self.histogram2.faceLighting = [];
                self.histogram2.lineStyle = [];
                self.histogram2.lineWidth = [];
                self.histogram2.enabled = [];
                self.histogram2.nbins = [];
                self.histogram2.xedges = [];
                self.histogram2.yedges = [];
            end

            if self.type.is.contour || self.type.is.contourf || self.type.is.contour3
                self.newprop("resolution", 2^9);
                self.newprop("maxnoise", 0.001);
            end

            if  self.type.is.contour
                self.newprop("contour", struct());
                self.contour.color = [];
                self.contour.labelSpacing = [];
                self.contour.lineWidth = [];
                self.contour.showText = [];
                self.contour.enabled = [];
                self.contour.levels = [];
                self.contour.lineSpec = [];
            end

            if  self.type.is.contourf
                self.newprop("contourf", struct());
                self.contourf.edgecolor = [];
                self.contourf.color = [];
                self.contourf.labelSpacing = [];
                self.contourf.lineWidth = [];
                self.contourf.showText = [];
                self.contourf.enabled = [];
                self.contourf.levels = [];
                self.contourf.lineSpec = [];
            end

            if  self.type.is.contour3
                self.newprop("contour3", struct());
                self.contour3.color = [];
                self.contour3.labelSpacing = [];
                self.contour3.lineWidth = [];
                self.contour3.showText = [];
                self.contour3.enabled = [];
                self.contour3.levels = [];
                self.contour3.lineSpec = [];
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% line/scatter attributes
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if  self.type.is.line || self.type.is.lineScatter
                self.newprop("surface", struct());
                self.newprop("plot", struct());
                self.plot.color = [];
                self.plot.lineWidth = [];
                self.plot.lineStyle = [];
                self.plot.markeredgeColor = [];
                self.plot.markerFaceColor = [];
                self.plot.markerIndices = [];
                self.plot.markerSize = [];
                self.plot.enabled = [];
            end

            if  self.type.is.line3 || self.type.is.lineScatter3
                self.newprop("surface", struct());
                self.newprop("plot3", struct());
                self.plot3.color = {};
                self.plot3.lineWidth = {};
                self.plot3.lineStyle = '-';
                self.plot3.markeredgeColor = [];
                self.plot3.markerFaceColor = [];
                self.plot3.markerIndices = [];
                self.plot3.markerSize = [];
                self.plot3.enabled = [];
            end

            if  isprop(self, "surface")
                self.surface.alignVertexCenters = {};
                self.surface.backFaceLighting = {};
                self.surface.edgeAlpha = {};
                self.surface.edgeColor = {};
                self.surface.edgeLighting = {};
                self.surface.faceAlpha = {};
                self.surface.faceColor = {};
                self.surface.faceLighting = {};
                self.surface.lineStyle = {};
                self.surface.lineWidth = {};
                self.surface.marker = {};
                self.surface.markerEdgeColor = {};
                self.surface.markerFaceColor = {};
                self.surface.markerSize = {};
                self.surface.meshStyle = {};
                self.surface.enabled = [];
            end

            if  self.type.is.scatter || self.type.is.lineScatter
                self.newprop("scatter", struct());
                self.scatter.colorVariable = [];
                self.scatter.lineWidth = [];
                self.scatter.markeredgeColor = [];
                self.scatter.markerFaceColor = [];
                self.scatter.enabled = [];
                self.scatter.color = [];
                self.scatter.filled = [];
                self.scatter.marker = [];
                self.scatter.size = [];
            end

            if  self.type.is.scatter3 || self.type.is.lineScatter3
                self.newprop("scatter3", struct());
                self.scatter3.colorVariable = [];
                self.scatter3.lineWidth = [];
                self.scatter3.markeredgeColor = [];
                self.scatter3.markerFaceColor = [];
                self.scatter3.enabled = [];
                self.scatter3.color = [];
                self.scatter3.filled = [];
                self.scatter3.marker = [];
                self.scatter3.size = [];
            end

            %if  self.type.is.targetable
            %    self.target = pm.vis.Target();
            %end

            if ~isempty(varargin)
                self.hash2comp(varargin); % parse arguments
            end
            %self.premake(varargin{:}); % This is the subclass method!

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Configure the plot settings and specifications and return nothing.<br>
        %>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``premake()`` method.<br>
        %>
        %>  \interface{premake}
        %>  \code{.m}
        %>
        %>      a = pm.vis.axes.Axes(ptype);
        %>
        %>      a.premake(varargin);
        %>      a.premake();
        %>
        %>  \endcode
        %>
        %>  \warning
        %>  This method has side-effects by manipulating
        %>  the existing attributes of the parent object.<br>
        %>
        %>  \example{premake}
        %>  \code{.m}
        %>
        %>      a = pm.vis.axes.Axes("line");
        %>      a.premake("xlim", [0, 1])
        %>
        %>  \endcode
        %>
        %>  \final{premake}
        %>
        %>  \author
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function premake(self, varargin)

            if ~isempty(varargin)
                self.hash2comp(varargin); % parse arguments
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% These settings must happen here so that they can be reset every time user nullifies the values.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%
            %%%% Set the enabled visualization components.
            %%%%

            if  isprop(self, "axes") && isempty(self.axes.enabled)
                self.axes.enabled = true;
            end

            if  isprop(self, "title") && isempty(self.title.enabled)
                self.title.enabled = ~isempty(self.title.titletext);
            end

            if  isprop(self, "xlabel") && isempty(self.xlabel.enabled)
                self.xlabel.enabled = true;
            end

            if  isprop(self, "ylabel") && isempty(self.ylabel.enabled)
                self.ylabel.enabled = true;
            end

            if  isprop(self, "zlabel") && isempty(self.zlabel.enabled)
                self.zlabel.enabled = true;
            end

            if  isprop(self, "colorbar") && isempty(self.colorbar.enabled)
                self.colorbar.enabled = true;
            end

            if  isprop(self, "colormap") && isempty(self.colormap.enabled)
                self.colormap.enabled = true;
            end

            if  isprop(self, "contour") && isempty(self.contour.enabled)
                self.contour.enabled = true;
            end

            if  isprop(self, "contourf") && isempty(self.contourf.enabled)
                self.contourf.enabled = true;
            end

            if  isprop(self, "contour3") && isempty(self.contour3.enabled)
                self.contour3.enabled = true;
            end

            if  isprop(self, "heatmap") && isempty(self.heatmap.enabled)
                self.heatmap.enabled = true;
            end

            if  isprop(self, "histfit") && isempty(self.histfit.enabled)
                self.histfit.enabled = true;
            end

            if  isprop(self, "histogram") && isempty(self.histogram.enabled)
                self.histogram.enabled = true;
            end

            if  isprop(self, "histogram2") && isempty(self.histogram2.enabled)
                self.histogram2.enabled = true;
            end

            if  isprop(self, "legend") && isempty(self.legend.enabled)
                self.legend.enabled = false;
            end

            if  isprop(self, "plot") && isempty(self.plot.enabled)
                self.plot.enabled = true;
            end

            if  isprop(self, "plot3") && isempty(self.plot3.enabled)
                self.plot3.enabled = true;
            end

            if  isprop(self, "scatter") && isempty(self.scatter.enabled)
                self.scatter.enabled = true;
            end

            if  isprop(self, "scatter3") && isempty(self.scatter3.enabled)
                self.scatter3.enabled = true;
            end

            %%%%
            %%%% ensure line plots are monochromatic when there is scatter plot.
            %%%%

            if  self.type.has.line
                if  self.type.has.scatter && isempty(self.surface.enabled)
                    self.surface.enabled = false;
                    % if  self.type.is.d2
                    %     self.plot.color = [];
                    % else
                    %     self.plot3.color = [];
                    % end
                elseif  isempty(self.surface.enabled)
                    self.surface.enabled = true;
                end
                proplist = ["plot", "plot3"];
                for prop = proplist
                    if  isprop(self, prop) && isempty(self.(prop).enabled)
                        self.(prop).enabled = ~self.surface.enabled;
                    end
                end
            end

            %%%%
            %%%% Set heatmap settings.
            %%%%

            if  self.type.is.heatmap
                if ~isempty(self.precision) && ~isnumeric(self.precision)
                    help("pm.vis.axes.Axes");
                    disp("self.precision");
                    disp( self.precision );
                    error   ( newline ...
                            + "The specified precision must be a positive whole number." + newline ...
                            + "For more information, see the documentation displayed above." + newline ...
                            + newline ...
                            );
                end
                self.setKeyVal(self.type.name, "colorbarVisible", "on");
                self.setKeyVal(self.type.name, "missingDataColor", "[0.1500 0.1500 0.1500]");
            end

            %%%%
            %%%% Set histfit/histogram/histogram2/contour/contour3/contourf settings.
            %%%%

            if  self.type.is.histogram
                self.setKeyVal(self.type.name, "edgeColor", "none");
            end

            if  self.type.is.histfit
                self.setKeyVal(self.type.name, "dist", "Normal");
                self.setKeyVal(self.type.name, "nbins", []);
            end

            if  self.type.is.histogram2
                self.setKeyVal(self.type.name, "edgeColor", "none");
                self.setKeyVal(self.type.name, "displayStyle", "bar3");
                self.setKeyVal(self.type.name, "showEmptyBins", "off");
               %self.setKeyVal(self.type.name, "numbins", [100 100]);
                if ~self.colormap.enabled
                    if ~pm.introspection.istype(self.histogram2.faceColor, "string", 1)
                        self.setKeyVal(self.type.name, "faceColor", "auto");
                    elseif self.histogram2.faceColor == "flat"
                        % enforce monochrome by removing the colormapping.
                        self.histogram2.faceColor = "auto";
                    end
                else
                    self.histogram2.faceColor = "flat";
                   %self.setKeyVal(self.type.name, "faceColor", "flat");
                end
            end

            if  self.type.is.diffusion
                if  self.type.is.contour3
                    self.setKeyVal(self.type.name, "color", []);
                    self.setKeyVal(self.type.name, "levels", 50);
                else
                    if  self.type.is.contourf
                        self.setKeyVal(self.type.name, "color", "none");
                        self.setKeyVal(self.type.name, "edgecolor", "none");
                        self.setKeyVal(self.type.name, "levels", 50); % 8
                        %if ~self.colormap.enabled && self.contourf.enabled
                        %    self.colormap.map = flipud(gray);
                        %end
                    else
                        self.setKeyVal(self.type.name, "levels", 10); % 8
                        self.setKeyVal(self.type.name, "color", []);
                    end
                end
                self.setKeyVal("resolution", 2^9);
                self.setKeyVal("maxnoise", 0.001);
                self.setKeyVal(self.type.name, "labelSpacing", 144);
                self.setKeyVal(self.type.name, "showText", "off");
                self.setKeyVal(self.type.name, "lineStyle", "-");
                self.setKeyVal(self.type.name, "lineWidth", 0.5);
            end

            %%%%
            %%%% Set line/scatter settings.
            %%%%

            if  self.type.is.scatter || self.type.is.lineScatter
                self.setKeyVal("scatter", "size", 5);
                self.setKeyVal("scatter", "marker", "o");
                self.setKeyVal("scatter", "filled", true);
                if  self.type.is.lineScatter && (self.scatter.enabled || self.colormap.enabled)
                    self.setKeyVal("plot", "color", uint8([200 200 200 150]));
                end
            end

            if  self.type.is.scatter3 || self.type.is.lineScatter3
                self.setKeyVal("scatter3", "size", 5);
                self.setKeyVal("scatter3", "marker", "o");
                self.setKeyVal("scatter3", "filled", true);
                if  self.type.is.lineScatter3 && (self.scatter3.enabled || self.colormap.enabled)
                    self.setKeyVal("plot3", "color", uint8([200 200 200 150]));
                end
            end

            if  self.type.has.line
                if self.type.is.d2; self.setKeyVal("plot", "lineWidth", 1); end
                if self.type.is.d3; self.setKeyVal("plot3", "lineWidth", 1); end
                self.setKeyVal("surface", "faceColor", "none");
                self.setKeyVal("surface", "edgeColor", "flat");
                self.setKeyVal("surface", "edgeAlpha", 0.5);
                self.setKeyVal("surface", "lineStyle", "-");
                self.setKeyVal("surface", "lineWidth", 1);
                self.setKeyVal("surface", "marker", "none");
            end

            %%%%
            %%%% Set the coloring settings.
            %%%%

            self.cenabled = false;
            self.cenabled = self.cenabled || (self.type.has.scatter && self.colormap.enabled);
            self.cenabled = self.cenabled || (self.type.has.line && self.surface.enabled && self.colormap.enabled);
            self.cenabled = self.cenabled || (self.type.is.histogram2 && self.colormap.enabled);
            self.cenabled = self.cenabled || (self.type.is.diffusion && self.colormap.enabled);
            self.cenabled = self.cenabled || (self.type.is.heatmap && self.colormap.enabled);

            % if ~self.cenabled && self.type.has.line
            %     if self.type.is.d2; self.plot.enabled = true; end
            %     if self.type.is.d3; self.plot3.enabled = true; end
            % end

            if isprop(self, "axes")
                if isfield(self.axes, "box") && isempty(self.axes.box); self.axes.box = "on"; end
                if isfield(self.axes, "xgrid") && isempty(self.axes.xgrid); self.axes.xgrid = "on"; end
                if isfield(self.axes, "ygrid") && isempty(self.axes.ygrid); self.axes.ygrid = "on"; end
                if isfield(self.axes, "zgrid") && isempty(self.axes.zgrid); self.axes.zgrid = "on"; end
            end

            %if  (self.type.is.diffusion || self.type.is.histogram2) && isempty(self.zlabel.txt)
            %    self.zlabel.txt = "Density";
            %end

        end % function

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Convert the components of the input component ``comp``
        %>  of the parent object into a cell array of key-val pairs.<br>
        %>
        %>  \details
        %>  This is a dynamic method of the class [pm.vis.axes.Axes](@ref Axes).<br>
        %>  This method is used internally by the subclasses to convert the parent object
        %>  attributes to input arguments of MATLAB intrinsic visualization functions.<br>
        %>
        %>  \param[in]  comp    :   The input scalar MATLAB string representing the name of a ``struct``
        %>                          component of the parent object, whose fields names and values are to
        %>                          be returned as subsequent pairs in the output ``hash`` cell array.<br>
        %>
        %>  \return
        %>  ``hash``            :   The output cell array containing the pairs of ``field-name, field-value``
        %>                          of the input MATLAB struct ``comp``.<br>
        %>
        %>  \interface{comp2hash}
        %>  \code{.m}
        %>
        %>      a = pm.vis.axes.Axes(ptype);
        %>
        %>      hash = a.comp2hash(comp);
        %>
        %>  \endcode
        %>
        %>  \warning
        %>  This method has side-effects by manipulating
        %>  the existing attributes of the parent object.<br>
        %>
        %>  \example{comp2hash}
        %>  \code{.m}
        %>
        %>      a = pm.vis.axes.Axes("histogram", "histogram", {"nbins", 5});
        %>      hash = a.comp2hash(a.histogram)
        %>
        %>  \endcode
        %>
        %>  \final{comp2hash}
        %>
        %>  \author
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function hash = comp2hash(self, comp)

            excludes = {"enabled"};
            if strcmp(comp, "axes")
                excludes = [excludes(:); "labels"; "parent"; "ncol"; "nrow"];
            elseif strcmp(comp, "colorbar")
                excludes = [excludes(:); "width"; "height"];
            elseif strcmp(comp, "colormap")
                excludes = [excludes(:); "map"];
            elseif strcmp(comp, "contour3") || strcmp(comp, "contourf") || strcmp(comp, "contour")
                % "color" is a keyword, but we set it as an explicit argument.
                excludes = [excludes(:); "levels"; "lineSpec"; "color"];
           %elseif strcmp(comp, "heatmap")
           %    excludes = [excludes(:)];
            elseif strcmp(comp, "histfit")
                excludes = [excludes(:); "dist"; "nbins"];
            elseif strcmp(comp, "histogram")
                excludes = [excludes(:); "edges"; "nbins"];
            elseif strcmp(comp, "histogram2")
                excludes = [excludes(:); "nbins"; "xedges"; "yedges"];
            elseif strcmp(comp, "legend")
                excludes = [excludes(:); "labels"];
            elseif strcmp(comp, "plot") || strcmp(comp, "plot3")
                excludes = [excludes(:); "size"]; % color;
            elseif strcmp(comp, "scatter") || strcmp(comp, "scatter3")
                excludes = [excludes(:); "marker"; "filled"; "color"; "size"];
            elseif strcmp(comp, "surface")
                excludes = [excludes(:); "color"; "size"];
            elseif strcmp(comp, "title")
                excludes = [excludes(:); "titletext"; "subtitletext"];
            elseif strcmp(comp, "xlabel") || strcmp(comp, "ylabel") || strcmp(comp, "zlabel")
                excludes = [excludes(:); "txt"];
            elseif ~strcmp(comp, "heatmap")
                disp("comp");
                disp( comp );
                error   ( newline ...
                        + "Internal library error: Unrecognized MATLAB function name" + newline ...
                        + "as ``comp`` argument of object of class ``pm.vis.axes.Axes``." + newline ...
                        + newline ...
                        );
            end
            unique = true;
            onlyfull = true;
            hash = pm.matlab.hashmap.struct2hash(self.(comp), excludes, unique, onlyfull);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end