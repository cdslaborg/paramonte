%>  \brief
%>  This is the abstract class for generating instances of axes
%>  with various types of plots from one or more columns of data.<br>
%>
%>  \details
%>  While it can be directly used for visualization tasks,
%>  this class primarily serves as the superclass for the
%>  visualization subclasses accessible to the end users.<br>
%>
%>  \note
%>  In case of [pm.vis.SubplotContour](@ref SubplotContour), [pm.vis.SubplotContour3](@ref SubplotContour3), and
%>  [pm.vis.SubplotContourf](@ref SubplotContourf) visualization objects, the density of the input data
%>  is first computed and smoothed by a Gaussian kernel density estimator and then passed to the
%>  MATLAB intrinsic plotting functions ``contour``, ``contourf``, ``contour3``.<br>
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
%>      <li>    ``colx`` (standing for x-columns; available for all axes types)<br>
%>
%>              Optional property that determines the columns of
%>              the specified dataframe to serve as the x-values.
%>              It can have multiple forms:<br>
%>              <ol>
%>                  <li>    a numeric or cell array of column indices in the input ``dfref``.
%>                  <li>    a string or cell array of column names in ``dfref.Properties.VariableNames``.
%>                  <li>    a cell array of a mix of the above two.
%>                  <li>    a numeric range.
%>              </ol>
%>              If ``colx`` is empty,<br>
%>              <ol>
%>                  <li>  it will be set to the row indices of ``dfref`` for line/scatter axes types.
%>                  <li>  it will be set to all columns of ``dfref`` for density axes types.
%>              </ol>
%>
%>              \warning
%>              In all cases, ``colx`` must have a length that is either
%>              0 (empty), 1, or equal to the length of ``coly`` or ``colz``.<br>
%>              If the length is 1, then ``colx`` will be plotted against
%>              data corresponding to each element of ``coly`` and ``colz``.<br>
%>
%>              \example{colx}
%>              \code{.m}
%>
%>                  Subplot.colx = [7, 8, 9]
%>                  Subplot.colx = ["sampleLogFunc", "sampleVariable1"]
%>                  Subplot.colx = {"sampleLogFunc", 9, "sampleVariable1"}
%>                  Subplot.colx = 7:9      # every column in the data frame starting from column #7 to #9
%>                  Subplot.colx = 7:2:20   # every other column in the data frame starting from column #7 to #20
%>
%>              \endcode
%>              <br>
%>
%>      <li>    ``coly`` (standing for y-columns; available for all axes types except
%>              [pm.vis.SubplotHeatmap](@ref SubplotHeatmap), [pm.vis.SubplotHistogram](@ref SubplotHistogram), [pm.vis.SubplotHistfit](@ref SubplotHistfit))<br>
%>
%>              Optional property that determines the columns of
%>              the specified dataframe to serve as the z-values.<br>
%>              It can have multiple forms:<br>
%>              <ol>
%>                  <li>    a numeric or cell array of column indices in the input ``dfref``.
%>                  <li>    a string or cell array of column names in ``dfref.Properties.VariableNames``.
%>                  <li>    a cell array of a mix of the above two.
%>                  <li>    a numeric range.
%>              </ol>
%>              If ``coly`` is empty,<br>
%>              <ol>
%>                  <li>    it will be set to the row indices of ``dfref`` for line/scatter axes types.
%>                  <li>    it will be set to all columns of ``dfref`` for density axes types.
%>              </ol>
%>
%>              \warning
%>              In all cases, ``coly`` must have a length that is either
%>              0 (empty), 1, or equal to the length of ``colx`` or ``colz``.<br>
%>              If the length is 1, then ``coly`` will be plotted against
%>              data corresponding to each element of ``colx`` and ``colz``.<br>
%>
%>              \example{coly}
%>              \code{.m}
%>
%>                  Subplot.coly = [7, 8, 9]
%>                  Subplot.coly = ["sampleLogFunc", "sampleVariable1"]
%>                  Subplot.coly = {"sampleLogFunc", 9, "sampleVariable1"}
%>                  Subplot.coly = 7:9      # every column in the data frame starting from column #7 to #9
%>                  Subplot.coly = 7:2:20   # every other column in the data frame starting from column #7 to #20
%>
%>              \endcode
%>              <br>
%>
%>      <li>    ``colz`` (standing for z-columns; available only for 3D plots, e.g.,
%>              [pm.vis.SubplotLine3](@ref SubplotLine3), [pm.vis.SubplotScatter3](@ref SubplotScatter3),
%>              [pm.vis.SubplotLineScatter3](@ref SubplotLineScatter3))<br>
%>
%>              Optional property that determines the columns of
%>              the specified dataframe to serve as the z-values.<br>
%>              It can have multiple forms:<br>
%>              <ol>
%>                  <li>    a numeric or cell array of column indices in the input ``dfref``.
%>                  <li>    a string or cell array of column names in ``dfref.Properties.VariableNames``.
%>                  <li>    a cell array of a mix of the above two.
%>                  <li>    a numeric range.
%>              </ol>
%>              If ``colz`` is empty,<br>
%>              <ol>
%>                  <li>    it will be set to the row indices of ``dfref`` for line/scatter axes types.
%>                  <li>    it will be set to all columns of ``dfref`` for density axes types.
%>              </ol>
%>
%>              \warning
%>              In all cases, ``colz`` must have a length that is either
%>              0 (empty), 1, or equal to the length of ``colx`` or ``coly``.<br>
%>              If the length is 1, then ``colz`` will be plotted against
%>              data corresponding to each element of ``colx`` and ``coly``.<br>
%>
%>              \example{colz}
%>              \code{.m}
%>
%>                  Subplot.colz = [7, 8, 9]
%>                  Subplot.colz = ["sampleLogFunc", "sampleVariable1"]
%>                  Subplot.colz = {"sampleLogFunc", 9, "sampleVariable1"}
%>                  Subplot.colz = 7:9      % every column in the data frame starting from column #7 to #9
%>                  Subplot.colz = 7:2:20   % every other column in the data frame starting from column #7 to #20
%>
%>              \endcode
%>              <br>
%>
%>      <li>    ``colc`` (standing for color-columns; available only for 2D/3D line/scatter axes types, e.g.,
%>              [pm.vis.SubplotLine3](@ref SubplotLine3), [pm.vis.SubplotScatter3](@ref SubplotScatter3),
%>              [pm.vis.SubplotLineScatter3](@ref SubplotLineScatter3), [pm.vis.SubplotLine](@ref SubplotLine),
%>              [pm.vis.SubplotScatter](@ref SubplotScatter), [pm.vis.SubplotLineScatter](@ref SubplotLineScatter))<br>
%>
%>              Optional property that determines the columns of the input ``dfref`` to
%>              use as the color mapping values for each line/point element in the plot.<br>
%>              It can have multiple forms:<br>
%>              <ol>
%>                  <li>    a numeric or cell array of column indices in the input ``dfref``.
%>                  <li>    a string or cell array of column names in ``dfref.Properties.VariableNames``.
%>                  <li>    a cell array of a mix of the above two.
%>                  <li>    a numeric range.
%>              </ol>
%>              The default value is the indices of the rows of the input ``dfref``.<br>
%>              \example{colc}
%>              \code{.m}
%>
%>                  Subplot.colc = [7,8,9]
%>                  Subplot.colc = ["sampleLogFunc", "sampleVariable1"]
%>                  Subplot.colc = {"sampleLogFunc", 9, "sampleVariable1"}
%>                  Subplot.colc = 7:9      # every column in the data frame starting from column #7 to #9
%>                  Subplot.colc = 7:2:20   # every other column in the data frame starting from column #7 to #20
%>
%>              \endcode
%>              <br>
%>  </ol>
%>
%>  \note
%>  For example usage, see the documentation of the constructor of
%>  the class [pm.vis.Subplot::Subplot](@ref Subplot::Subplot).<br>
%>
%>  \devnote
%>  While dynamic addition of class attributes is not ideal, the current
%>  design was deemed unavoidable and best, given the constraints of the
%>  MATLAB language and visualization tools.<br>
%>
%>  \see
%>  [pm.vis](@ref \psldir/main/+pm/+vis)<br>
%>  [pm.vis.axes.Axes](@ref Axes)<br>
%>  [pm.matlab.Handle](@ref Handle)<br>
%>  [pm.vis.Cascade](@ref Cascade)<br>
%>  [pm.vis.Subplot](@ref Subplot)<br>
%>  [pm.vis.Triplex](@ref Triplex)<br>
%>  [pm.vis.Figure](@ref Figure)<br>
%>  [pm.vis.Plot](@ref Plot)<br>
%>  [pm.vis.Tile](@ref Tile)<br>
%>
%>  \final
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef Subplot < pm.vis.axes.Axes

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = public)
        %>
        %>  ``dfref``
        %>
        %>  A scalar object of class [pm.container.DataFrame](@ref DataFrame)
        %>  containing (a reference to) the user-specified data to visualize.<br>
        %>
        dfref = [];
        %>
        %>  ``fout``
        %>
        %>  A MATLAB ``struct`` whose fields are the outputs of
        %>  various plotting tools used to make the current axis.<br>
        %>
        fout = struct();
        %>
        %>  ``rows``
        %>
        %>  A numeric vector that serves as a storage for an arbitrary subset of indices
        %>  of the rows of the input dataframe reference ``dfref`` to the class constructor .<br>
        %>  It can be either:<br>
        %>  <ol>
        %>      <li>    a numeric range, or,
        %>      <li>    a list of row indices of the ``dfref``.
        %>  </ol>
        %>
        %>  \example{rows}
        %>  \code{fout}
        %>
        %>      s = pm.vis.Subplot();
        %>
        %>      s.rows = 15 : -2 : 8;
        %>      s.rows = [12, 46, 7, 8, 9, 4, 7, 163];
        %>
        %>  \endcode
        %>
        %>  \warning
        %>  If ``rows`` is empty, the default will be all rows of the ``dfref``.<br>
        %>
        %>  \note
        %>  The [pm.container.DataFrame.rowslog()](@ref DataFrame::rowslog) method of
        %>  this class can be used to generate logarithmically-spaced
        %>  row indices of the target dataframe.<br>
        %>
        rows = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Construct and return an object of class [pm.vis.Subplot](@ref Subplot).<br>
        %>
        %>  \details
        %>  This is the custom constructor of the class [pm.vis.Subplot](@ref Subplot).<br>
        %>
        %>  \param[in]  ptype       :   See the documentation of the corresponding input argument of the
        %>                              superclass constructor [pm.vis.axes.Axes::Axes](@ref Axes::Axes).<br>
        %>  \param[in]  dfref       :   The input MATLAB matrix or table containing the data to plot or
        %>                              a function handle that returns such a MATLAB matrix or table.<br>
        %>                              Specifying a function handle is superior to specifying the
        %>                              data directly, because the function handle will always use
        %>                              the most updated version of the user table or matrix.<br>
        %>                              (**optional**. The default is an empty table.)
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.<br>
        %>                              The input ``varargin`` can also contain the components
        %>                              of the ``subplot`` component of the parent object.<br>
        %>
        %>  \return
        %>  ``self``            :   The output scalar object of class [pm.vis.Subplot](@ref Subplot).<br>
        %>
        %>  \interface{Subplot}
        %>  \code{.m}
        %>
        %>      s = pm.vis.Subplot(ptype, dfref);
        %>      s = pm.vis.Subplot(ptype, dfref, varargin);
        %>
        %>  \endcode
        %>
        %>  \warning
        %>  In all cases, ``colc`` must have a length that is either
        %>  0, or 1, or equal to the maximum lengths of ``(colx, coly, colz)``.<br>
        %>  If the length is 1, then all data will be plotted with the same color
        %>  mapping determined by values specified by the elements of ``colc``.<br>
        %>  If it is an empty object having length 0, then the default value will be used.<br>
        %>
        %>  \note
        %>  In case of [pm.vis.SubplotContour](@ref SubplotContour)/[pm.vis.SubplotContourf](@ref SubplotContourf)/[pm.vis.SubplotContour3](@ref SubplotContour3)
        %>  visualization objects, the density of the input data is first computed and smoothed
        %>  by a Gaussian kernel density estimator and then passed to the MATLAB intrinsic
        %>  plotting functions , ``contour``, ``contourf``, ``contour3``.<br>
        %>
        %>  \example{Subplot}
        %>  \include{lineno} example/vis/Subplot/main.m
        %>  \vis{Subplot}
        %>  <br>\image html example/vis/Subplot/SubplotLine.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotLine3.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotScatter.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotScatter3.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotLineScatter.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotLineScatter3.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotHistfit.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotHistogram.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotHistogram2.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotContour.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotContourf.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotContour3.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotHeatmap.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotEllipse3.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotEllipse.1.png width=700
        %>
        %>  \final{Subplot}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 22 2024, 6:45 PM, University of Texas at Arlington<br>
        %>  \AmirShahmoradi, July 5 2024, 1:07 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        function self = Subplot(ptype, dfref, varargin)
            if  nargin < 2
                dfref = [];
            end
            if  nargin < 1
                ptype = [];
            end
            self@pm.vis.axes.Axes(ptype, varargin{:});
            self.dfref = pm.container.DataFrame(dfref);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Generate a plot from the selected columns
        %>  of the parent object component ``dfref``.<br>
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
        %>      s = pm.vis.Subplot(ptype, dfref)
        %>
        %>      s.make(varargin);
        %>      s.make();
        %>
        %>  \endcode
        %>
        %>  \warning
        %>  This method has side-effects by manipulating
        %>  the existing attributes of the parent object.<br>
        %>
        %>  \example{make}
        %>  \code{.m}
        %>
        %>      dfref = rand(1000, 3);
        %>      s = pm.vis.Subplot("scatter", dfref);
        %>      s.make("colx", 1, "coly", 2, "colc", 3)
        %>
        %>  \endcode
        %>
        %>  \example{Subplot}
        %>  \include{lineno} example/vis/Subplot/main.m
        %>  \vis{Subplot}
        %>  <br>\image html example/vis/Subplot/SubplotLine.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotLine3.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotScatter.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotScatter3.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotLineScatter.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotLineScatter3.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotHistfit.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotHistogram.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotHistogram2.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotContour.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotContourf.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotContour3.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotHeatmap.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotEllipse3.1.png width=700
        %>  <br>\image html example/vis/Subplot/SubplotEllipse.1.png width=700
        %>
        %>  \final{make}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 22 2024, 6:36 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, July 5 2024, 1:07 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        function make(self, varargin)

            self.premake(varargin{:});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: No component of ``self`` is allowed to appear to the left of assignment operator, except ``fout``.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%% Set what to plot.

            noPlotEnabled = false;
            noPlotEnabled = noPlotEnabled || (self.type.is.scatter && ~self.scatter.enabled);
            noPlotEnabled = noPlotEnabled || (self.type.is.scatter3 && ~self.scatter3.enabled);
            noPlotEnabled = noPlotEnabled || (self.type.is.line && ~self.plot.enabled && ~self.surface.enabled);
            noPlotEnabled = noPlotEnabled || (self.type.is.line3 && ~self.plot3.enabled && ~self.surface.enabled);
            noPlotEnabled = noPlotEnabled || (self.type.is.lineScatter && ~self.plot.enabled && ~self.surface.enabled && ~self.scatter.enabled);
            noPlotEnabled = noPlotEnabled || (self.type.is.lineScatter3 && ~self.plot3.enabled && ~self.surface.enabled && ~self.scatter3.enabled);
            noPlotEnabled = noPlotEnabled || (self.type.is.histogram2 && ~self.histogram2.enabled);
            noPlotEnabled = noPlotEnabled || (self.type.is.histogram && ~self.histogram.enabled);
            noPlotEnabled = noPlotEnabled || (self.type.is.contourf && ~self.contourf.enabled);
            noPlotEnabled = noPlotEnabled || (self.type.is.contour3 && ~self.contour3.enabled);
            noPlotEnabled = noPlotEnabled || (self.type.is.histfit && ~self.histfit.enabled);
            noPlotEnabled = noPlotEnabled || (self.type.is.contour && ~self.contour.enabled);
            noPlotEnabled = noPlotEnabled || (self.type.is.heatmap && ~self.heatmap.enabled);
            if  noPlotEnabled
                error   ( newline ...
                        + "All plot types are deactivated by the user. There is nothing to display." + newline ...
                        + "Enable at least one the following plotting components of the parent object." + newline  ...
                        + newline ...
                        );
            end

            dfcopy = self.dfref.copy();
            if  isempty(dfcopy)
                help("pm.vis.Subplot");
                error   ( newline ...
                        + "There is no data to plot." + newline ...
                        + "Without any input data, visualizations are impossible." + newline ...
                        + "The input argument ``dfref`` to the object constructor was likely empty." + newline ...
                        + "If you merely want to create a template of visualization specifications," + newline ...
                        + "then use the class [pm.vis.axes.Axes](@ref Axes) which is the superclass of this class." + newline ...
                        + "For more information, see the class documentation displayed above." + newline ...
                        + newline ...
                        );
            end

            %%%%
            %%%% check rows presence
            %%%%

            if pm.array.len(self.rows)
                rowsindex = self.rows(:);
            else
                rowsindex = [1 : self.dfref.nrow()]';
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% check columns presence
            %%%%%%%%%%%%%%%%%%%%%%%%%%%

            coldatx = []; colindx = []; colnamx = [];
            if  isprop(self, "colx") && pm.array.len(self.colx)
                [colindx, colnamx] = pm.str.locname(dfcopy.Properties.VariableNames, self.colx);
            elseif self.type.is.heatmap || self.type.is.d1
                colnamx = dfcopy.Properties.VariableNames;
                colindx = 1 : length(colnamx);
            else
                coldatx = rowsindex;
                colnamx = "Row Index";
            end

            coldaty = []; colindy = []; colnamy = [];
            if  isprop(self, "coly") && pm.array.len(self.coly)
                [colindy, colnamy] = pm.str.locname(dfcopy.Properties.VariableNames, self.coly);
            elseif self.type.is.density && ~self.type.is.histogram2
                if  isfield(self.(self.type.name), "normalization") && ~isempty(self.(self.type.name).normalization)
                    colnamy = string(self.(self.type.name).normalization);
                elseif self.type.is.d1
                    colnamy = "Count";
                else
                    colnamy = "Density";
                end
            elseif  self.type.is.heatmap
                if  isempty(dfcopy.Properties.RowNames)
                    dfcopy.Properties.RowNames = string(rowsindex);
                end
                colnamy = dfcopy.Properties.RowNames;
            else
                coldaty = rowsindex;
                colnamy = "Row Index";
            end

            coldatz = []; colindz = []; colnamz = [];
            if  isprop(self, "colz") && pm.array.len(self.colz)
                [colindz, colnamz] = pm.str.locname(dfcopy.Properties.VariableNames, self.colz);
            elseif  self.type.has.line || self.type.has.scatter
                if ~self.type.is.d3
                    % This is necessary for 2D surface plot.
                    coldatz = zeros(length(rowsindex), 1);
                else
                    coldatz = rowsindex;
                    colnamz = "Row Index";
                end
            elseif  self.type.is.density && self.type.is.triaxes
                colnamz = "Density";
            end

            coldatc = []; colindc = []; colnamc = [];
            if  isprop(self, "colc") && pm.array.len(self.colc)
                [colindc, colnamc] = pm.str.locname(dfcopy.Properties.VariableNames, self.colc);
            elseif self.type.has.line || self.type.has.scatter
                if ~isempty(colnamz) && ~(isempty(coldatz) && isempty(colindz))
                    coldatc = coldatz;
                    colindc = colindz;
                    colnamc = colnamz;
                else
                    coldatc = rowsindex;
                    colnamc = "Row Index";
                end
            elseif self.type.is.density && ~self.type.is.d1
                colnamc = "Density";
            end

            %%%% check the lengths are consistent.

            colindlen = max([ length(colindc) ...
                            , length(colindz) ...
                            , length(colindy) ...
                            , length(colindx) ...
                            ]);

            if length(colindc) ~= colindlen && 1 < length(colindc); error("length of colx must be either 1 or equal to the maximum of the lengths of coly, colz, colc."); end
            if length(colindz) ~= colindlen && 1 < length(colindz); error("length of coly must be either 1 or equal to the maximum of the lengths of colx, colz, colc."); end
            if length(colindy) ~= colindlen && 1 < length(colindy); error("length of colz must be either 1 or equal to the maximum of the lengths of colx, coly, colc."); end
            if length(colindx) ~= colindlen && 1 < length(colindx); error("length of colc must be either 1 or equal to the maximum of the lengths of colx, coly, colz."); end

            %%%% assign data in case of single column assignments.

            if length(colindx) == 1; coldatx = dfcopy{rowsindex, colindx}; end
            if length(colindy) == 1; coldaty = dfcopy{rowsindex, colindy}; end
            if length(colindz) == 1; coldatz = dfcopy{rowsindex, colindz}; end
            if length(colindc) == 1; coldatc = dfcopy{rowsindex, colindc}; end

            %%%% get keyword arguments.

            kws = struct();
            for prop =  [ "axes" ...
                        , "colorbar" ...
                        , "colormap" ...
                        , "contour" ...
                        , "contour3" ...
                        , "contourf" ...
                        , "legend" ...
                        , "heatmap" ...
                        , "histfit" ...
                        , "histogram" ...
                        , "histogram2" ...
                        , "plot" ...
                        , "plot3" ...
                        , "scatter" ...
                        , "scatter3" ...
                        , "surface" ...
                        , "title" ...
                        , "xlabel" ...
                        , "ylabel" ...
                        , "zlabel" ...
                        ]
                if  isprop(self, prop)
                    kws.(prop) = self.comp2hash(prop);
                end
            end

            for prop = ["scatter", "scatter3"]
                if  isprop(self, prop) && self.(prop).filled
                    kws.(prop) = [{"filled"}; kws.(prop)(:)];
                end
            end

            %%%% generate the kde2d density estimates for the contour plots.

            kde2dUpdateNeeded = self.type.is.diffusion && (1 < length(colindx) || 1 < length(colindy));
            if  self.type.is.diffusion
                kde2d = struct();
                if ~kde2dUpdateNeeded
                    [kde2d.bandwidth, kde2d.density, kde2d.crdx, kde2d.crdy] = pm.stats.hist.kde2d([coldatx(:), coldaty(:)], self.resolution);
                    kde2d.density(kde2d.density < self.maxnoise) = NaN;
                end
            end

            %%%% determine the multicolor nature of the plots.

            isMultiColorPlot = false;
            if ~self.cenabled && (self.type.has.scatter || self.type.is.diffusion)
                if  self.type.is.diffusion
                    prop = self.type.name;
                elseif self.type.is.d2
                    prop = "scatter";
                elseif self.type.is.d3
                    prop = "scatter3";
                else
                    error("Internal error occurred. A line/scatter plot must be either 2D or 3D.");
                end
                if  pm.array.len(self.(prop).color)
                    if  all(size(self.(prop).color) == [colindlen, 3])
                        isMultiColorPlot = true;
                        multiColorList = self.(prop).color;
                    elseif all(size(self.(prop).color) == [1, 3]) || self.type.is.diffusion
                        coldatc = self.(prop).color;
                    else
                        error   ( newline ...
                                + "The value specified for the 'color' property of the " + prop + " component" + newline ...
                                + "must be either an RGB triplet vector of size [1, 3], or a matrix of size [" + string(colindlen) + ", 3]" + newline ...
                                + "for the current set of variables that are selected to plot. " + newline ...
                                + "It can also be an empty object, in which case, the colors" + newline ...
                                + "of the objects in the plot will be chosen automatically." + newline ...
                                + newline ...
                                );
                    end
                else
                    multiColorList = lines(colindlen);
                    isMultiColorPlot = true;
                end
            elseif self.type.is.diffusion && ~pm.array.len(self.(self.type.name).color)
                coldatc = self.(self.type.name).color;
            end

            %%%% initiate the legend.

            legendary = isprop(self, "legend") && self.legend.enabled;
            if  legendary
                legempty = pm.array.len(self.legend.labels) == 0;
                lglabels = strings(colindlen, 1);
                icollgx = 1;
                icollgy = 1;
                icollgz = 1;
            end

            %%%% generate axes and set axes properties.

            if  isprop(self, "axes")
                self.fout.axes = gca;
                if  self.axes.enabled
                    % This setting must happen here, before any
                    % individual axes component settings below.
                    set(gca, kws.axes{:});
                end
                if  self.type.is.triaxes
                    view(3);
                end
                hold on;
            end

            %%%%
            %%%% Ensure points and lines are scaled properly.
            %%%% It is not clear whether it's worth doing so because
            %%%% most figure resizing happen after making it,
            %%%% which makes the following solution useless.
            %%%%

            % temp = [];
            % if ~strcmpi(self.fout.axes.Units, "normalized")
            %     temp = self.fout.axes.Units;
            %     self.fout.axes.Units = "normalized";
            % end
            % defaultSizeFig = [560, 420] % width, height in pixels.
            % defaultSizeAxes = [0.775, 0.815] % width, height in pixels (assuming no colorbar).
            % self.fout.axes.Position
            % if ~isempty(temp)
            %     self.fout.axes.Units = temp;
            % end

            %%%% Make plots.

            if  self.type.is.heatmap && self.heatmap.enabled

                if ~isempty(self.precision)
                    self.fout.heatmap = heatmap(colnamx, colnamy, round(dfcopy{colnamy, colnamx}, self.precision), kws.heatmap{:});
                else
                    self.fout.heatmap = heatmap(colnamx, colnamy, dfcopy{colnamy, colnamx}, kws.heatmap{:});
                end

            else

                if ~self.type.has.line && ~self.type.has.scatter
                    % line/scatter plots must be taken care of separately.
                    self.fout.(self.type.name) = cell(colindlen, 1);
                end

                for icol = 1 : colindlen

                    %%%% set the data to plot.

                    if 1 < length(colindx); if legendary; icollgx = icol; end; coldatx = dfcopy{rowsindex, colindx(icol)}; end
                    if 1 < length(colindy); if legendary; icollgy = icol; end; coldaty = dfcopy{rowsindex, colindy(icol)}; end
                    if 1 < length(colindz); if legendary; icollgz = icol; end; coldatz = dfcopy{rowsindex, colindz(icol)}; end

                    if  isMultiColorPlot
                        coldatc = multiColorList(icol, :);
                    elseif 1 < length(colindc)
                        coldatc = dfcopy{rowsindex, colindc(icol)};
                    end

                    if  self.type.is.diffusion
                        if ~isempty(coldatc)
                            kws.(self.type.name) = pm.matlab.hashmap.repKeyVal("color", coldatc, kws.(self.type.name));
                        end
                        if  kde2dUpdateNeeded
                            [kde2d.bandwidth, kde2d.density, kde2d.crdx, kde2d.crdy] = pm.stats.hist.kde2d([coldatx(:), coldaty(:)], self.resolution);
                            kde2d.density(kde2d.density < self.maxnoise) = NaN;
                        end
                    end

                    %%%% set the default legend.

                    if  legendary && legempty
                        if  1 < length(colindz)
                            suffix = "-" + colnamz(icollgz);
                        else
                            suffix = "";
                        end
                        if  length(colindx) < 2 && 1 < length(colindy)
                            lglabels(icol) = colnamy(icollgy) + suffix;
                        elseif 1 < length(colindx) && length(colindy) < 2
                            lglabels(icol) = colnamx(icollgx) + suffix;
                        else
                            lglabels(icol) = colnamx(icollgx) + "-" + colnamy(icollgy) + suffix;
                        end
                    end

                    %%%% add density plot.

                    if  self.type.is.contour
                        self.fout.contour{icol} = struct();
                        [self.fout.contour{icol}.matrix, self.fout.contour{icol}.handle] = contour  ( kde2d.crdx ...
                                                                                                    , kde2d.crdy ...
                                                                                                    , kde2d.density ...
                                                                                                    , self.contour.levels ...
                                                                                                    , kws.contour{:} ...
                                                                                                    );
                    elseif self.type.is.contour3
                        self.fout.contour3{icol} = struct();
                        [self.fout.contour3{icol}.matrix, self.fout.contour3{icol}.handle]= contour3( kde2d.crdx...
                                                                                                    , kde2d.crdy...
                                                                                                    , kde2d.density...
                                                                                                    , self.contour3.levels ...
                                                                                                    , kws.contour3{:} ...
                                                                                                    );
                    elseif self.type.is.contourf
                        self.fout.contourf{icol} = struct();
                        [self.fout.contourf{icol}.matrix, self.fout.contourf{icol}.handle]= contourf( kde2d.crdx ...
                                                                                                    , kde2d.crdy ...
                                                                                                    , kde2d.density ...
                                                                                                    , self.contourf.levels ...
                                                                                                    , kws.contourf{:} ...
                                                                                                    );
                    elseif self.type.is.histfit
                        self.fout.histfit{icol} = histfit(coldatx, self.histfit.nbins, self.histfit.dist, kws.histfit{:});
                        self.fout.histfit{icol}(1).FaceAlpha = 0.6; % make histogram face transparent 60%
                        self.fout.histfit{icol}(1).EdgeAlpha = 0; % make histogram lines invisible
                        self.fout.histfit{icol}(2).Color = "black"; % make the fitted line black
                    elseif self.type.is.histogram
                        if ~isempty(self.histogram.nbins)
                            self.fout.histogram{icol} = histogram(coldatx, self.histogram.nbins, kws.histogram{:});
                        elseif ~isempty(self.histogram.edges)
                            self.fout.histogram{icol} = histogram(coldatx, self.histogram.edges, kws.histogram{:});
                        else
                            self.fout.histogram{icol} = histogram(coldatx, kws.histogram{:});
                        end
                    elseif self.type.is.histogram2
                        self.fout.histogram2{icol} = histogram2(coldatx, coldaty, kws.histogram2{:});
                    end

                    %%%% add line plot.

                    if  self.type.has.line
                        if  self.type.is.d2 && self.plot.enabled
                            self.fout.plot{icol} = plot ( coldatx ...
                                                        , coldaty ...
                                                        , kws.plot{:} ...
                                                        );
                        elseif self.type.is.d3 && self.plot3.enabled
                            self.fout.plot3{icol} = plot3   ( coldatx ...
                                                            , coldaty ...
                                                            , coldatz ...
                                                            , kws.plot3{:} ...
                                                            );
                        end
                        if  self.surface.enabled && self.colormap.enabled
                            self.fout.surface{icol} = surface   ( "XData", [coldatx(:), coldatx(:)] ...
                                                                , "YData", [coldaty(:), coldaty(:)] ...
                                                                , "ZData", [coldatz(:), coldatz(:)] ...
                                                                , "CData", [coldatc(:), coldatc(:)] ...
                                                                , kws.surface{:} ...
                                                                );
                        end
                    end

                    %%%% add scatter plot.

                    if  self.type.has.scatter
                        if  self.type.is.d3 && self.scatter3.enabled
                            self.fout.scatter3{icol} = scatter3 ( coldatx ...
                                                                , coldaty ...
                                                                , coldatz ...
                                                                , self.scatter3.size ...
                                                                , coldatc ...
                                                                , self.scatter3.marker ...
                                                                , kws.scatter3{:} ...
                                                                );
                        elseif self.type.is.d2 && self.scatter.enabled
                            self.fout.scatter{icol} = scatter   ( coldatx ...
                                                                , coldaty ...
                                                                , self.scatter.size ...
                                                                , coldatc ...
                                                                , self.scatter.marker ...
                                                                , kws.scatter{:} ...
                                                                );
                        end
                    end

                end % loop over plots

            end % if heatmap plot.

            % add axis labels.
            if ~self.type.is.heatmap
                if isprop(self, "xlabel") && self.xlabel.enabled; if ~pm.array.len(self.xlabel.txt); txt = colnamx; else; txt = self.xlabel.txt; end; self.fout.xlabel = xlabel(txt, kws.xlabel{:}); end
                if isprop(self, "ylabel") && self.ylabel.enabled; if ~pm.array.len(self.ylabel.txt); txt = colnamy; else; txt = self.ylabel.txt; end; self.fout.ylabel = ylabel(txt, kws.ylabel{:}); end
                if isprop(self, "zlabel") && self.zlabel.enabled; if ~pm.array.len(self.zlabel.txt); txt = colnamz; else; txt = self.zlabel.txt; end; self.fout.zlabel = zlabel(txt, kws.zlabel{:}); end
            else
                if isprop(self, "xlabel") && self.xlabel.enabled; xlabel(self.xlabel.txt, kws.xlabel{:}); end
                if isprop(self, "ylabel") && self.ylabel.enabled; ylabel(self.ylabel.txt, kws.ylabel{:}); end
            end

            % add the colormap and colorbar.

            if ~self.cenabled && self.type.is.heatmap
                colormap(flipud(gray)); % self.fout.colormap =
            elseif ~self.cenabled && self.type.is.contourf
                colormap(self.fout.axes, flipud(gray)); % self.fout.colormap =
            elseif  self.cenabled
                if ~isempty(self.colormap.map) || (isstring(self.colormap.map) && self.colormap.map ~= "")
                    cmap = self.colormap.map;
                else
                    if  self.type.has.line || self.type.has.scatter
                        cmap = "winter";
                        %if  self.type.is.d3
                        %    self.colormap.map = "winter";
                        %else
                        %    self.colormap.map = "autumn";
                        %end
                    elseif self.type.is.heatmap
                        cmap = pm.vis.cmap.redblue();
                    elseif self.type.is.density
                        cmap = "parula";
                    end
                end
                if ~self.type.is.heatmap
                    colormap(self.fout.axes, cmap); % self.fout.colormap =
                else
                    colormap(cmap); % self.fout.colormap =
                end
            end

            if  isprop(self, "colorbar")
                if  self.colorbar.enabled && (self.cenabled || self.type.is.contourf)
                    if  isempty(self.colorbar.fontSize) %|| ~isa(self.colorbar.fontSize, "numeric")
                        kws.colorbar = [kws.colorbar(:); {"fontSize"}; {self.fout.axes.FontSize + 1}];
                    end
                    self.fout.colorbar = colorbar(kws.colorbar{:});
                    ylabel(self.fout.colorbar, colnamc(1));%, "fontSize", self.colorbar.fontSize, "interpreter", "none"
                else
                    colorbar('off');
                    self.fout.colorbar = [];
                end
            end

            % add legend.

            if  legendary && 0 < colindlen
                self.fout.legend = [];
                if ~legempty
                    lglabels = self.legend.labels;
                end
                if ~self.legend.enabled
                    legend(self.fout.axes, "off");
                elseif isa(self.legend.labels, "cell") || isa(self.legend.labels, "string")
                    if  isempty(self.legend.fontSize)
                        kws.legend = [kws.legend(:); {"fontSize"}; {self.fout.axes.FontSize + 1}];
                    end
                    subset = [];
                    if  self.type.has.scatter && self.type.is.d2 && self.scatter.enabled
                        subset = [self.fout.scatter{:}];
                    elseif self.type.has.scatter && self.type.is.d3 && self.scatter3.enabled
                        subset = [self.fout.scatter3{:}];

                    elseif self.type.is.histfit && self.histfit.enabled
                        nhandle = length(self.fout.histfit);
                        subset = cell(nhandle, 1);
                        for i = 1 : nhandle
                            subset{i} = self.fout.histfit{i}(1); % get the handle to the histogram component of histfit.
                        end
                        subset = [subset{:}];
                    end
                    if ~isempty(subset)
                        self.fout.legend = legend(subset, lglabels{:}, kws.legend{:});
                    else
                        self.fout.legend = legend(lglabels{:}, kws.legend{:});
                    end
                else
                    error   ( newline ...
                            + "The specified ""legend.labels"" must be a cell array of string values." ...
                            + newline ...
                            );
                end
            end

            %if  self.type.is.targetable
            %    self.target.fout.axes = self.fout.axes;
            %end

            %%%% add title.

            if  self.title.enabled
                titletext = {};
                for field = ["titletext", "subtitletext"] % do not change the order of elements here.
                    if  isfield(self.title, field) && pm.array.len(self.title.(field))
                        titletext = [titletext(:); self.title.(field)(:)];
                    end
                end
                if ~isempty(titletext)
                    if ~self.type.is.heatmap
                        %if ~isempty(kws)
                        %    kws.title = {titletext{:}, kws.title(:)};
                        %end
                        self.fout.title = title(self.fout.axes, titletext{:}, kws.title{:});
                    else
                        %%%% Neither input nor output arguments are not supported when using title with heatmap.
                        title(titletext);
                    end
                end
            end

            %%%% set limits.

            for prop = ["xlim", "ylim", "zlim"]
                if  isprop(self, prop) && ~isempty(self.(prop))
                    limit = get(self.fout.axes, prop);
                    for i = 1 : 2
                        if ~isnan(self.(prop)(i))
                            limit(i) = self.(prop)(i);
                        end
                    end
                    set(self.fout.axes, prop, limit);
                end
            end

            %%%% set axes scales.

            for prop = ["xscale", "yscale", "zscale"]
                if  isprop(self, prop) && ~isempty(self.(prop))
                    set(self.fout.axes, prop, self.(prop));
                end
            end

            if ~self.type.is.heatmap
                hold off;
                %grid on;
                %box on;
            end

        end % function

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Prepare the subplot for visualization.<br>
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
        %>      s = pm.vis.Subplot(ptype, dfref);
        %>
        %>      s.premake(varargin);
        %>      s.premake();
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
        %>      dfref = rand(1000, 3);
        %>      s = pm.vis.Subplot("scatter", dfref);
        %>      s.premake("colx", 1, "coly", 2, "colc", 3)
        %>
        %>  \endcode
        %>
        %>  \final{premake}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 22 2024, 6:41 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, July 5 2024, 1:07 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        function premake(self, varargin)
            premake@pm.vis.axes.Axes(self, varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Reset the properties of the plot to the original default settings.<br>
        %>
        %>  \brief
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
        %>      pm.vis.Subplot.reset() % reset the plot to the default settings.
        %>
        %>  \endcode
        %>
        %>  \final{reset}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 22 2024, 6:43 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, July 5 2024, 1:07 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        function reset(self, varargin)
            self.rows = [];
            self.fout = struct();
            self.newprop("colx", {});
            if  self.type.is.d2 || self.type.is.d3
                self.newprop("coly", {});
                if  self.type.is.d3
                    self.newprop("colz", {});
                end
            end
            reset@pm.vis.axes.Axes(self, varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

end