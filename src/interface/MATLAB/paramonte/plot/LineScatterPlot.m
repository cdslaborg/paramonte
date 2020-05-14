%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ParaMonte: plain powerful parallel Monte Carlo library.
%
%   Copyright (C) 2012-present, The Computational Data Science Lab
%
%   This file is part of the ParaMonte library.
%
%   ParaMonte is free software: you can redistribute it and/or modify it
%   under the terms of the GNU Lesser General Public License as published
%   by the Free Software Foundation, version 3 of the License.
%
%   ParaMonte is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public License
%   along with the ParaMonte library. If not, see,
%
%       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
%
%   ACKNOWLEDGMENT
%
%   As per the ParaMonte library license agreement terms,
%   if you use any parts of this library for any purposes,
%   we ask you to acknowledge the ParaMonte library's usage
%   in your work (education/research/industry/development/...)
%   by citing the ParaMonte library as described on this page:
%
%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   LineScatterPlot(dataFrame, plotType)
%
%   This is the LineScatterPlot class for generating instances of
%   line/scatter/lineScatter/line3/scatter3/lineScatter3/
%   based on a wide variety of MATLAB's builtin functions.
%
%   NOTE: This is a low-level ParaMonte class and is not meant
%   NOTE: to be directly instantiated by the user.
%
%   Parameters
%   ----------
%
%       dataFrame
%
%           a MATLAB data Table from which the selected data will be plotted.
%           This is a low-level internal argument and is not meant
%           to be accessed or be provided by the user.
%
%       plotType
%
%           a string representing the type of the plot to be made.
%           This is a low-level internal argument and is not meant
%           to be accessed or be provided by the user.
%
%   Attributes
%   ----------
%
%       ycolumns
%
%           optional property that determines the columns of dataFrame to serve as
%           the y-values. It can have multiple forms:
%
%               1.  a numeric or cell array of column indices in the input dataFrame.
%               2.  a string or cell array of column names in dataFrame.Properties.VariableNames.
%               3.  a cell array of a mix of the above two.
%               4.  a numeric range.
%
%           Example usage:
%
%               1.  ycolumns = [7,8,9]
%               2.  ycolumns = ["SampleLogFunc","SampleVariable1"]
%               3.  ycolumns = {"SampleLogFunc",9,"SampleVariable1"}
%               4.  ycolumns = 7:9      # every column in the data frame starting from column #7 to #9
%               5.  ycolumns = 7:2:20   # every other column in the data frame starting from column #7 to #20
%
%           WARNING: In all cases, ycolumns must have a length that is either 0, or 1, or equal
%           WARNING: to the length of xcolumns. If the length is 1, then ycolumns will be
%           WARNING: plotted against data corresponding to each element of ycolumns.
%           WARNING: If it is an empty object having length 0, then the default value will be used.
%
%           The default value is the names of all columns of the input dataFrame.
%
%       zcolumns (available only in line3/scatter3/lineScatter3 objects)
%
%           optional property that determines the columns of dataFrame to serve as
%           the z-values. It can have multiple forms:
%
%               1.  a numeric or cell array of column indices in the input dataFrame.
%               2.  a string or cell array of column names in dataFrame.Properties.VariableNames.
%               3.  a cell array of a mix of the above two.
%               4.  a numeric range.
%
%           Example usage:
%
%               1.  zcolumns = [7,8,9]
%               2.  zcolumns = ["SampleLogFunc","SampleVariable1"]
%               3.  zcolumns = {"SampleLogFunc",9,"SampleVariable1"}
%               4.  zcolumns = 7:9      # every column in the data frame starting from column #7 to #9
%               5.  zcolumns = 7:2:20   # every other column in the data frame starting from column #7 to #20
%
%           WARNING: In all cases, zcolumns must have a length that is either 0, or 1, or equal
%           WARNING: to the length of xcolumns. If the length is 1, then zcolumns will be
%           WARNING: plotted against data corresponding to each element of zcolumns.
%           WARNING: If it is an empty object having length 0, then the default value will be used.
%
%           The default value is the names of all columns of the input dataFrame.
%
%       ccolumns (standing for color-columns)
%
%           optional property that determines the columns of dataFrame to serve
%           as the color-mapping-values for each line/point element in the plot.
%           It can have multiple forms:
%
%               1.  a numeric or cell array of column indices in the input dataFrame.
%               2.  a string or cell array of column names in dataFrame.Properties.VariableNames.
%               3.  a cell array of a mix of the above two.
%               4.  a numeric range.
%
%           Example usage:
%
%               1.  ccolumns = [7,8,9]
%               2.  ccolumns = ["SampleLogFunc","SampleVariable1"]
%               3.  ccolumns = {"SampleLogFunc",9,"SampleVariable1"}
%               4.  ccolumns = 7:9      # every column in the data frame starting from column #7 to #9
%               5.  ccolumns = 7:2:20   # every other column in the data frame starting from column #7 to #20
%
%           WARNING: In all cases, ccolumns must have a length that is either 0, or 1, or equal
%           WARNING: to the maximum lengths of (xcolumns,ycolumns,zcolumns). If the length is 1, then all data
%           WARNING: will be plotted with the same color mapping determined by values specified by the elements
%           WARNING: of ccolumns. If it is an empty object having length 0, then the default value will be used.
%
%           The default value is the indices of the rows of the input dataFrame.
%
%       colormap
%
%           A MATLAB struct() property with two components:
%
%               1. enabled: logical value. If `true`, the colormap will be applied to the plot
%               1. values: a string or any other value that the colormap function of MATLAB accepts as input.
%
%           Example usage:
%
%               1.  colormap = "autumn"
%               1.  colormap = "winter"
%
%           If colormap is not provided or is empty, the default will be "autumn".
%
%       colorbar_kws
%
%           A MATLAB struct() whose components' values are passed to MATLAB's colorbar function.
%           If your desired attribute is missing from the fieldnames of colorbar_kws, simply add
%           a new field named as the attribute and assign the desired value to it.
%
%           Example usage:
%
%               colorbar_kws.enabled = true % add colorbar
%               colorbar_kws.location = "west"
%
%           If a desired property is missing among the struct fields, simply add the field
%           and its value to colorbar_kws.
%
%           WARNING: keep in mind that MATLAB keyword arguments are case-INsensitive.
%           WARNING: therefore make sure you do not add the keyword as multiple different fields.
%           WARNING: For example, colorbar_kws.color and colorbar_kws.Color are the same,
%           WARNING: and only one of the two will be processed.
%
%       plot_kws (available only in line/lineScatter/line3/lineScatter3 objects)
%
%           A MATLAB struct() whose fields (with the exception of few, e.g., enabled, singleOptions, ...)
%           are directly passed to the `plot()` function of MATLAB (or plot3() if plot is 3d).
%           This property exists only if the object is instantiated as a line/lineScatter/line3/lineScatter3 object.
%
%           Example usage:
%
%               plot_kws.enabled = true; % add plot()
%               plot_kws.linewidth = 2;
%
%           If a desired property is missing among the struct fields, simply add the field
%           and its value to plot_kws.
%
%           WARNING: keep in mind that MATLAB keyword arguments are case-INsensitive.
%           WARNING: therefore make sure you do not add the keyword as multiple different fields.
%           WARNING: For example, plot_kws.color and plot_kws.Color are the same,
%           WARNING: and only one of the two will be processed.
%
%       surface_kws (available only in line/lineScatter/line3/lineScatter3 objects)
%
%           A MATLAB struct() whose fields (with the exception of few, e.g., enabled, singleOptions, ...)
%           are directly passed to the `surface()` function of MATLAB.
%           This property exists only if the object is instantiated as a line/lineScatter/line3/lineScatter3 object.
%           This property is used only if the line plot is colored via colormap.
%
%           Example usage:
%
%               surface_kws.enabled = true; % add surface()
%               surface_kws.linewidth = 2;
%
%           If a desired property is missing among the struct fields, simply add the field
%           and its value to surface_kws.
%
%           WARNING: keep in mind that MATLAB keyword arguments are case-INsensitive.
%           WARNING: therefore make sure you do not add the keyword as multiple different fields.
%           WARNING: For example, surface_kws.linewidth and surface_kws.Linewidth are the same,
%           WARNING: and only one of the two will be processed.
%
%       scatter_kws (available only in scatter/lineScatter/scatter3/lineScatter3 objects)
%
%           A MATLAB struct() whose fields (with the exception of few, e.g., enabled, singleOptions, ...)
%           are directly passed to the `scatter()` function of MATLAB.
%           This property exists only if the object is instantiated as a scatter/lineScatter/scatter3/lineScatter3 object.
%
%           Example usage:
%
%               scatter_kws.enabled = true; % add scatter()
%               scatter_kws.marker = ".";
%               scatter_kws.size = 10;
%
%           If a desired property is missing among the struct fields, simply add the field
%           and its value to scatter_kws.
%
%           WARNING: keep in mind that MATLAB keyword arguments are case-INsensitive.
%           WARNING: therefore make sure you do not add the keyword as multiple different fields.
%           WARNING: For example, scatter_kws.marker and scatter_kws.Marker are the same,
%           WARNING: and only one of the two will be processed.
%
%       target
%
%           an object of class Target_class for adding target values to the plots.
%           For more information, see the documentation for Target.
%
%   Superclass Attributes
%   ----------------------
%
%       See the documentation for the BasePlot class
%
%   Returns
%   -------
%
%       an object of LineScatterPlot class
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
classdef LineScatterPlot < BasePlot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Access = public)

        ccolumns
        colorbar_kws
        colormap
        target

    end

    properties (Access = protected, Hidden)
        %dfref = [];
        %isdryrun = [];
        is3d
        plotType
        isLinePlot = false;
        isScatterPlot = false;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function reset(self)

            reset@BasePlot(self);
            prop="xcolumns"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            prop="ycolumns"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            self.xcolumns = {};
            self.ycolumns = {};
            self.ccolumns = {};
            self.colormap = [];

            self.isLinePlot = false;
            if contains(self.plotType,"line")
                self.isLinePlot = true;
                prop="plot_kws"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                prop="surface_kws"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            end
            self.isScatterPlot = false;
            if contains(self.plotType,"scatter")
                self.isScatterPlot = true;
                prop="scatter_kws"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            end
            if contains(self.plotType,"3")
                prop="zcolumns"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                self.is3d = true;
            else
                self.is3d = false;
            end

            if self.isLinePlot

                self.plot_kws = struct();
                self.plot_kws.enabled = true;
                self.plot_kws.linewidth = {};
                self.plot_kws.singleOptions = {};
                self.plot_kws.color = {};

                self.surface_kws = struct();
                if self.isScatterPlot
                    self.surface_kws.enabled = false;
                    self.plot_kws.color = uint8([200 200 200 150]);;
                else
                    self.surface_kws.enabled = true;
                    self.surface_kws.facecolor = "none";
                    self.surface_kws.edgecolor = "flat";
                    self.surface_kws.edgealpha = 0.5;
                    self.surface_kws.linestyle = "-";
                    self.surface_kws.marker = "none";
                end

            end

            if self.isScatterPlot
                self.scatter_kws = struct();
                self.scatter_kws.marker = [];
                self.scatter_kws.singleOptions = {};
                self.scatter_kws.size = [];
                self.scatter_kws.enabled = true;
            end

            self.colorbar_kws = struct();
            self.colorbar_kws.enabled = true;
            self.colorbar_kws.fontsize = [];

            self.isdryrun = true;
            self.plot();
            self.isdryrun = false;

            self.target = Target_class();

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = LineScatterPlot(varargin) % expected input arguments: dataFrame, plotType

            self = self@BasePlot(varargin{1});
            self.plotType = lower(string(varargin{2}));
            self.reset()
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function helpme(self,varargin)
            %
            %   Open the documentation for the input object's name in string format, otherwise, 
            %   open the documentation page for the class of the object owning the helpme() method.
            %
            %   Parameters
            %   ----------
            %
            %       This function takes at most one string argument, 
            %       which is the name of the object for which help is needed.
            %
            %   Returns
            %   -------
            %
            %       None. 
            %
            %   Example
            %   -------
            %
            %       helpme("plot")
            %
            methodNotFound = true;
            if nargin==2
                if strcmpi(varargin{1},"target")
                    cmd = "doc Target_class";
                    methodNotFound = false;
                else
                    methodList = ["plot","helpme"];
                    for method = methodList
                        if strcmpi(varargin{1},method)
                            methodNotFound = false;
                            cmd = "doc self." + method;
                        end
                    end
                end
            elseif nargin~=1
                error("The helpme() method takes at most one argument that must be string.");
            end
            if methodNotFound
                cmd = "doc self";
            end
            eval(cmd);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function plot(self,varargin)
            %
            %   Generate a plot from the selected columns of the object's dataFrame.
            %
            %   Parameters
            %   ----------
            %
            %       Any property,value pair of the object.
            %       If the property is a struct(), then its value must be given as a cell array,
            %       with consecutive elements representing the struct's property-name,property-value pairs.
            %       Note that all of these property-value pairs can be also directly set directly via the 
            %       object's attributes, before calling the plot() method.
            %
            %   Returns
            %   -------
            %
            %       None. However, this method causes side-effects by manipulating
            %       the existing attributes of the object.
            %
            %   Example
            %   -------
            %
            %       plot("ycolumns",[8])
            %       plot("ycolumns",7:10)
            %       plot( "gcf_kws", {"enabled",true,"color","none"} )
            %

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% parse arguments
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.parseArgs(varargin{:});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% set what to plot
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if self.isLinePlot
                % activate at least one plot in the figure
                if ~(self.plot_kws.enabled || self.surface_kws.enabled)
                    warning ( newline ...
                            + "Both line surface() and plot() plot types have been disabled by the user. There is nothing to display. " ...
                            + "To add at least one plot, set at least one the following components of the line-plot, " + newline  + newline  ...
                            + "To add at least one plot, set at least one the following components of the line-plot, " + newline  + newline  ...
                            + "    self.plot_kws.enabled = true; % to generate single color monochromatic line plots" + newline ...
                            + "    self.surface_kws.enabled = true; % to generate colorful, color-mapped line plots" + newline  + newline  ...
                            + "You can also pass these arguments at the time of calling the plot() method:" + newline  + newline  ...
                            + "    self.plot(""plot_kws"",{""enabled"",true}); % to generate single color monochromatic line plots" + newline ...
                            + "    self.plot(""surface_kws"",{""enabled"",true}); % to generate single color monochromatic line plots"+ newline ...
                            + newline ...
                            );
                end
            end

            cEnabled =  ( self.isLinePlot && self.surface_kws.enabled ) || ( self.isScatterPlot && self.scatter_kws.enabled );
            lgEnabled = self.legend_kws.enabled && ~cEnabled;

            if self.isScatterPlot && self.scatter_kws.enabled
                key = "marker"; val = ".";
                if isfield(self.scatter_kws,key) && isempty(self.scatter_kws.(key))
                    self.scatter_kws.(key) = val;
                end
                key = "size"; val = 10;
                if isfield(self.scatter_kws,key) && isempty(self.scatter_kws.(key))
                    self.scatter_kws.(key) = val;
                end
            end

            if self.isLinePlot && self.plot_kws.enabled
                key = "linewidth"; val = 1;
                if isfield(self.plot_kws,key) && isempty(self.plot_kws.(key))
                    self.plot_kws.(key) = val;
                end
            end

            if cEnabled && ~getVecLen(self.colormap)
                if self.is3d
                    self.colormap = "winter";
                else
                    self.colormap = "autumn";
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if self.isdryrun; return; end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % generate figure and axes if needed

            if self.gcf_kws.enabled
                gcf_kws_cell = convertStruct2Cell(self.gcf_kws,{"enabled","singleOptions"});
                if isfield(self.gcf_kws,"singleOptions"); gcf_kws_cell = { gcf_kws_cell{:} , self.gcf_kws.singleOptions{:} }; end
                self.currentFig.gcf = figure( gcf_kws_cell{:} );
            else
                set(0, "CurrentFigure", gcf);
                self.currentFig.gcf = gcf;
            end

            % check rows presence

            if getVecLen(self.rows)
                rowindex = self.rows;
            else
                rowindex = 1:1:length(self.dfref{:,1});
            end
            %ndata = length(rowindex);

            % check columns presence

            if getVecLen(self.xcolumns)
                [xcolnames, xcolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.xcolumns);
            else
                xcolindex = [];
                xcolnames = "Count";
                xdata = rowindex;
            end

            [ycolnames, ycolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.ycolumns);

            if self.is3d && getVecLen(self.zcolumns)
                [zcolnames, zcolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.zcolumns);
            else
                zcolindex = [];
                zcolnames = "Count";
                zdata = rowindex;
            end

            % set color data

            if cEnabled
                if getVecLen(self.ccolumns)
                    [ccolnames, ccolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.ccolumns);
                else
                    ccolindex = [];
                    ccolnames = "Count";
                    cdata = 1:1:length(self.dfref{rowindex,1});
                end
            else
                ccolindex = [];
                ccolnames = [];
                cdata = [];
            end

            % check the lengths are consistent

            xcolindexlen = length(xcolindex);
            ycolindexlen = length(ycolindex);
            zcolindexlen = length(zcolindex);
            ccolindexlen = length(ccolindex);
            maxLenColumns = max (   [ xcolindexlen ...
                                    , ycolindexlen ...
                                    , zcolindexlen ...
                                    , ccolindexlen ...
                                    ] ...
                                );

            if xcolindexlen~=maxLenColumns && xcolindexlen>1; error("length of xcolumns must be either 1 or equal to the maximum of the lengths of ycolumns, zcolumns, ccolumns."); end
            if ycolindexlen~=maxLenColumns && ycolindexlen>1; error("length of ycolumns must be either 1 or equal to the maximum of the lengths of xcolumns, zcolumns, ccolumns."); end
            if zcolindexlen~=maxLenColumns && zcolindexlen>1; error("length of zcolumns must be either 1 or equal to the maximum of the lengths of xcolumns, ycolumns, ccolumns."); end
            if ccolindexlen~=maxLenColumns && ccolindexlen>1; error("length of ccolumns must be either 1 or equal to the maximum of the lengths of xcolumns, ycolumns, zcolumns."); end

            % assign data in case of single column assignments

            if xcolindexlen==1
                xdata = self.dfref{rowindex,xcolindex};
            end
            if ycolindexlen==1
                ydata = self.dfref{rowindex,ycolindex};
            end
            if zcolindexlen==1
                zdata = self.dfref{rowindex,zcolindex};
            end
            if ccolindexlen==1
                cdata = self.dfref{rowindex,ccolindex};
            end

            %%%%%%%%%%%%%%%%%%%%%%%
            % get keyword arguments
            %%%%%%%%%%%%%%%%%%%%%%%

            if self.isLinePlot
                plot_kws_cell = convertStruct2Cell(self.plot_kws,{"enabled","singleOptions"});
            end
            if self.isScatterPlot
                scatter_kws_cell = convertStruct2Cell(self.scatter_kws,{"enabled","singleOptions","color","size"});
            end
            if self.isLinePlot
                surface_kws_cell = convertStruct2Cell(self.surface_kws,{"enabled","singleOptions","color","size"});
            end

            %%%%%%%%%%%%%%%%%%%%%%%
            % add line/scatter plot
            %%%%%%%%%%%%%%%%%%%%%%%

            box on; grid on;

            lglabels = [];
            if cEnabled
                colormap(self.colormap);
            else
                if self.isScatterPlot && self.scatter_kws.enabled
                    scatter_colors = lines(maxLenColumns);
                end
            end

            if (self.isScatterPlot && self.scatter_kws.enabled) || (self.isLinePlot && (self.surface_kws.enabled || self.plot_kws.enabled))

                if ~self.is3d
                    zdata = zeros(length(rowindex(:)),1);
                end

                for icol = 1:maxLenColumns

                    if xcolindexlen>1
                        xdata = self.dfref{rowindex,xcolindex(icol)};
                    end
                    if ycolindexlen>1
                        ydata = self.dfref{rowindex,ycolindex(icol)};
                    end
                    if zcolindexlen>1
                        zdata = self.dfref{rowindex,zcolindex(icol)};
                    end
                    if ccolindexlen>1
                        cdata = self.dfref{rowindex,ccolindex(icol)};
                    end
                    if lgEnabled && ~self.is3d
                        if xcolindexlen<2 && ycolindexlen>=1
                            lglabels = [ lglabels , ycolnames(icol) ];
                        elseif xcolindexlen>1 && ycolindexlen<2
                            lglabels = [ lglabels , xcolnames(icol) ];
                        else
                            lglabels = [ lglabels , xcolnames(icol)+"-"+ycolnames(icol) ];
                        end
                    end

                    % add plot

                    if self.isLinePlot
                        if self.plot_kws.enabled
                            if self.is3d
                                self.currentFig.plot3 = plot3   ( xdata ...
                                                                , ydata ...
                                                                , zdata ...
                                                                , plot_kws_cell{:} ...
                                                                );
                            else
                                self.currentFig.plot = plot ( xdata ...
                                                            , ydata ...
                                                            , plot_kws_cell{:} ...
                                                            );
                            end
                            hold on;
                        end
                        if self.surface_kws.enabled && getVecLen(self.colormap)
                            self.currentFig.surface = surface   ( "XData",[xdata(:) xdata(:)] ...
                                                                , "YData",[ydata(:) ydata(:)] ...
                                                                , "ZData",[zdata(:) zdata(:)] ...
                                                                , "CData",[cdata(:) cdata(:)] ...
                                                                , surface_kws_cell{:} ...
                                                                );
                            if self.is3d; view(3); end
                            hold on;
                        end
                    end

                    if self.isScatterPlot && self.scatter_kws.enabled
                        if cEnabled
                            if self.is3d
                                self.currentFig.scatter3 = scatter3 ( xdata ...
                                                                    , ydata ...
                                                                    , zdata ...
                                                                    , self.scatter_kws.size ...
                                                                    , cdata ...
                                                                    , scatter_kws_cell{:} ...
                                                                    );
                            else
                                self.currentFig.scatter = scatter   ( xdata ...
                                                                    , ydata ...
                                                                    , self.scatter_kws.size ...
                                                                    , cdata ...
                                                                    , scatter_kws_cell{:} ...
                                                                    );
                            end
                        else
                            if self.is3d
                                %plot_kws = {};
                                %if ~isa(self.plot_kws,"cell"); plot_kws = self.plot_kws;
                                self.currentFig.scatter3 = scatter3 ( xdata ...
                                                                    , ydata ...
                                                                    , zdata ...
                                                                    , self.scatter_kws.size ...
                                                                    , scatter_colors(icol,:) ...
                                                                    , scatter_kws_cell{:} ...
                                                                    );
                               %self.currentFig.plot = plot3( xdata ...
                               %                            , ydata ...
                               %                            , zdata ...
                               %                            , scatter_kws_cell{:} ...
                               %                            );
                            else
                                %plot_kws = {};
                                %if ~isa(self.plot_kws,"cell"); plot_kws = self.plot_kws;
                                self.currentFig.scatter = scatter   ( xdata ...
                                                                    , ydata ...
                                                                    , self.scatter_kws.size ...
                                                                    , scatter_colors(icol,:) ...
                                                                    , scatter_kws_cell{:} ...
                                                                    );
                               %self.currentFig.plot = plot ( xdata ...
                               %                            , ydata ...
                               %                            , scatter_kws_cell{:} ...
                               %                            );
                            end
                        end
                        hold on;
                    end

                end % loop plot

                self.currentFig.gca = gca;

            end

            % add axis labels

            if xcolindexlen>1
                self.currentFig.xlabel = xlabel("Variable Values");
            else
                self.currentFig.xlabel = xlabel(xcolnames(1));
            end

            if ycolindexlen>1
                self.currentFig.ylabel = ylabel("Variable Values");
            else
                self.currentFig.ylabel = ylabel(ycolnames(1));
            end

            if self.is3d
            if zcolindexlen>1
                self.currentFig.zlabel = zlabel("Variable Values");
            else
                self.currentFig.zlabel = zlabel(zcolnames(1));
            end
            end

            % add line colorbar

            if self.colorbar_kws.enabled && cEnabled
                if isempty(self.colorbar_kws.fontsize) || ~isa(self.colorbar_kws.fontsize,"numeric")
                    self.colorbar_kws.fontsize = self.currentFig.ylabel.FontSize;
                end
                colorbar_kws_cell = convertStruct2Cell(self.colorbar_kws,{"enabled","singleOptions"});
                self.currentFig.colorbar = colorbar(colorbar_kws_cell{:});
                ylabel(self.currentFig.colorbar,ccolnames(1),"fontsize",self.colorbar_kws.fontsize);
            else
                colorbar('off');
                self.currentFig.colorbar = [];
            end

            if ~self.is3d || (self.is3d && self.legend_kws.enabled)
                self.doBasePlotStuff(lgEnabled,lglabels)
            end

            grid on; hold off;

        end % function plot

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % classdef LinePlot