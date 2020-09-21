%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   MIT License
%%%%
%%%%   ParaMonte: plain powerful parallel Monte Carlo library.
%%%%
%%%%   Copyright (C) 2012-present, The Computational Data Science Lab
%%%%
%%%%   This file is part of the ParaMonte library.
%%%%
%%%%   Permission is hereby granted, free of charge, to any person obtaining a
%%%%   copy of this software and associated documentation files (the "Software"),
%%%%   to deal in the Software without restriction, including without limitation
%%%%   the rights to use, copy, modify, merge, publish, distribute, sublicense,
%%%%   and/or sell copies of the Software, and to permit persons to whom the
%%%%   Software is furnished to do so, subject to the following conditions:
%%%%
%%%%   The above copyright notice and this permission notice shall be
%%%%   included in all copies or substantial portions of the Software.
%%%%
%%%%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
%%%%   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
%%%%   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
%%%%   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
%%%%   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
%%%%   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
%%%%   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%
%%%%   ACKNOWLEDGMENT
%%%%
%%%%   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
%%%%   As per the ParaMonte library license agreement terms, if you use any parts of
%%%%   this library for any purposes, kindly acknowledge the use of ParaMonte in your
%%%%   work (education/research/industry/development/...) by citing the ParaMonte
%%%%   library as described on this page:
%%%%
%%%%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   LineScatterPlot(dataFrame, plotType, resetExternal)
%
%   This is the LineScatterPlot class for generating instances of
%   line/scatter/lineScatter/line3/scatter3/lineScatter3/
%   based on a wide variety of MATLAB's builtin functions.
%
%       NOTE
%
%           This is a low-level ParaMonte class and is not meant
%           to be directly instantiated by the user.
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
%       xcolumns
%
%           Optional property that determines the columns of dataFrame 
%           to serve as the x-values. It can have multiple forms:
%
%               1.  a numeric or cell array of column indices in the input dataFrame.
%               2.  a string or cell array of column names in dataFrame.Properties.VariableNames.
%               3.  a cell array of a mix of the above two.
%               4.  a numeric range.
%
%           Example usage:
%
%               1.  xcolumns = [7,8,9]
%               2.  xcolumns = ["SampleLogFunc","SampleVariable1"]
%               3.  xcolumns = {"SampleLogFunc",9,"SampleVariable1"}
%               4.  xcolumns = 7:9      # every column in the data frame starting from column #7
%               5.  xcolumns = 7:2:20   # every other column in the data frame starting from column #7
%
%           WARNING
%
%               In all cases, xcolumns must have a length that is either 0, or 1, 
%               or equal to the length of ycolumns. If the length is 1, then xcolumns 
%               will be plotted against data corresponding to each element of ycolumns.
%               If it is an empty object having length 0, then a default value will be used.
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
%           WARNING
%
%               In all cases, ycolumns must have a length that is either 0, or 1, 
%               or equal to the length of xcolumns. If the length is 1, then ycolumns 
%               will be plotted against data corresponding to each element of xcolumns.
%               If it is an empty object having length 0, then a default value will be used.
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
%           WARNING: to the length of xcolumns.
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
%       colormap (available only in histogram2/contourf/contour3/contour objects) 
%
%           A MATLAB struct() with the following components:
%
%               enabled
%
%                   A logical value. If `true`, the colormap will be applied 
%                   to the plot. 
%
%               values
%
%                   A string or any other value that the colormap function 
%                   of MATLAB accepts as input.
%
%           Example usage:
%
%               1.  colormap.enabled = true;
%               1.  colormap.values = "winter";
%               1.  colormap.values = "autumn";
%
%       colorbar (available only in histogram2/contourf/contour3/contour objects) 
%
%           A MATLAB struct() with the following components:
%
%               enabled
%
%                   A logical value. If `true`, the colorbar will be applied to the plot.
%
%               kws
%
%                   A MATLAB struct() whose components' values are passed to 
%                   MATLAB's colorbar() function. If your desired attribute is 
%                   missing from the fieldnames of colorbar.kws, simply add 
%                   a new field named as the attribute and assign the desired 
%                   value to it.
%
%           Example usage:
%
%               colorbar.enabled = true; % add colorbar
%               colorbar.kws.location = "west";
%
%           WARNING
%
%               Keep in mind that MATLAB keyword arguments are case-INsensitive.
%               Hence, sure you do not add the keyword as multiple different fields.
%               For example, colorbar.kws.color and colorbar.kws.Color are the same,
%               and only one of the two will be processed.
%
%       plot (available only in line/lineScatter/line3/lineScatter3 objects)
%
%           A MATLAB struct() with the following components: 
%
%               enabled
%
%                   A boolean indicating whether a call to the 
%                   ``colorbar()`` function of MATLAB should be made or not.
%                   If a call is made, a new colorbar will be added to the figure.
%
%               kws
%
%                   A MATLAB struct() whose fields are directly passed to 
%                   the ``plot()`` or ``plot3()`` function of MATLAB.
%
%           This property exists only if the object is instantiated 
%           as a line/lineScatter/line3/lineScatter3 object. 
%
%           Example usage:
%
%               plot.enabled = true; % enable plot()
%               plot.kws.lineWidth = 2;
%
%           If a desired property is missing among the struct fields, 
%           simply add the field and its value to ``plot.kws``.
%
%           WARNING
%
%               Keep in mind that MATLAB keyword arguments are case-INsensitive.
%               therefore make sure you do not add the keyword as multiple different fields.
%               For example, ````plot.kws.color`` and ``plot.kws.Color`` are the same, 
%               and only one of the two will be processed.
%
%       surface (available only in line/lineScatter/line3/lineScatter3 objects)
%
%           A MATLAB struct() with the following components: 
%
%               enabled
%
%                   A boolean indicating whether a call to the 
%                   ``surface()`` function of MATLAB should be made or not.
%
%               kws
%
%                   A MATLAB ``struct()`` whose fields are directly passed to 
%                   the ``surface()`` function of MATLAB.
%
%           This property exists only if the object is instantiated 
%           as a line/lineScatter/line3/lineScatter3 object. 
%
%           Example usage:
%
%               surface.enabled = true; % add surface()
%               surface.kws.lineWidth = 2;
%
%           If a desired property is missing among the struct fields, 
%           simply add the field and its value to ``surface.kws``.
%
%           WARNING
%
%               Keep in mind that MATLAB keyword arguments are case-INsensitive.
%               therefore make sure you do not add the keyword as multiple different fields.
%               For example, ````surface.kws.lineWidth`` and ``surface.kws.Linewidth`` are the same, 
%               and only one of the two will be processed.
%
%       scatter (available only in scatter/lineScatter/scatter3/lineScatter3 objects)
%
%           A MATLAB struct() with the following components: 
%
%               enabled
%
%                   A boolean indicating whether a call to the 
%                   ``scatter()`` function of MATLAB should be made or not.
%
%               kws
%
%                   A MATLAB ``struct()`` whose fields are directly passed to 
%                   the ``scatter()`` or ``scatter3()`` function of MATLAB.
%
%           This property exists only if the object is instantiated 
%           as a line/lineScatter/line3/lineScatter3 object. 
%
%           Example usage:
%
%               scatter.enabled = true; % add scatter()
%               scatter.kws.marker = ".";
%               scatter.size = 10; % set the point size 
%               scatter.color = "red"; % set the points color 
%
%           If a desired property is missing among the struct fields, simply add the field
%           and its value to scatter.kws.
%
%           WARNING
%
%               Keep in mind that MATLAB keyword arguments are case-INsensitive.
%               therefore make sure you do not add the keyword as multiple different fields.
%               For example, ````surface.kws.lineWidth`` and ``surface.kws.Linewidth`` are the same, 
%               and only one of the two will be processed.
%
%       target
%
%           An object of class Target_class for adding target values to the plots.
%           For more information, see the documentation for the Target_class().
%
%   Superclass Attributes
%   ---------------------
%
%       See the documentation for the ``BasePlot`` class.
%
%   Returns
%   -------
%
%       An object of ``LineScatterPlot`` class.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
classdef LineScatterPlot < BasePlot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Access = public)
        ccolumns
        colorbar
        colormap
    end

    properties (Hidden)
        % number of data points to plot
        ndata
        xlim = [];
        ylim = [];
        zlim = [];
        ccolnames
        ccolindex
        xcolnames
        xcolindex
        ycolnames
        ycolindex
        zcolnames
        zcolindex
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function resetInternal(self)

            resetInternal@BasePlot(self);
            prop="xcolumns"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            prop="ycolumns"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            self.xcolumns = {};
            self.ycolumns = {};
            self.ccolumns = {};

            self.colormap = struct();
            self.colormap.enabled = true;
            self.colormap.values = [];

            if self.type.isLine
                prop="plot"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                prop="surface"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            end

            if self.type.isScatter
                prop="scatter"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            end

            if self.type.is3d
                prop="zcolumns"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            end

            if self.type.isLine

                self.plot = struct();
                self.plot.enabled = true;
                self.plot.kws = struct();
                self.plot.kws.lineWidth = {};
                self.plot.kws.color = {};

                self.surface.kws = struct();
                if self.type.isScatter
                    self.surface.enabled = false;
                    self.plot.kws.color = [];
                else
                    self.surface.enabled = true;
                end

            end

            if self.type.isScatter
                self.scatter = struct();
                self.scatter.marker = [];
                self.scatter.filled = [];
                self.scatter.color = [];
                self.scatter.size = [];
                self.scatter.enabled = true;
                self.scatter.kws = struct();
                %NO%self.scatter.kws.markerEdgeColor = struct();
                %NO%self.scatter.kws.markerFaceColor = struct();
            end

            self.colorbar = struct();
            self.colorbar.enabled = true;
            self.colorbar.kws = struct();
            self.colorbar.kws.fontSize = [];

            self.isdryrun = true;
            self.make();
            self.isdryrun = false;

            self.legend.enabled = false;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = LineScatterPlot(plotType, dataFrame, resetExternal)
            if nargin<3; resetExternal = []; end
            self = self@BasePlot(plotType, dataFrame, resetExternal);
            if nargin<3; self.resetExternal = @self.resetInternal; end
            self.resetInternal();
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

        function make(self,varargin)
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
            %       object's attributes, before calling the make() method.
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
            %       make("ycolumns",[8])
            %       make("ycolumns",7:10)
            %

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% parse arguments
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.parseArgs(varargin{:});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% set what to plot
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            noPlotEnabled   =   (~self.type.isLine && self.type.isScatter) && ~(self.scatter.enabled) ...
                            ||  (self.type.isLine && ~self.type.isScatter) && ~(self.plot.enabled || self.surface.enabled) ...
                            ||  (self.type.isLine && self.type.isScatter) && ~(self.plot.enabled || self.surface.enabled || self.scatter.enabled);
            if noPlotEnabled
                warning ( newline ...
                        + "All plots are deactivated by the user. There is nothing to display. " ...
                        + "To add at least one plot, set at least one the following components of the line/scatter-plot, " + newline  ...
                        + newline ...
                        + "    self.plot.enabled = true; % to generate single color monochromatic line plots" + newline ...
                        + "    self.surface.enabled = true; % to generate colorful, color-mapped line plots" + newline ...
                        + "    self.scatter.enabled = true; % to generate colorful, color-mapped scatter plots" + newline ...
                        + newline ...
                        );
            end

            if self.type.isScatter % && self.scatter.enabled
                fname = "scatter";
                key = "marker"; val = "o"; if ~isfield(self.(fname),key) || isempty(self.(fname).(key)); self.(fname).(key) = val; end
                key = "filled"; val = true; if ~isfield(self.(fname),key) || isempty(self.(fname).(key)); self.(fname).(key) = val; end
                key = "size"; val = 5; if ~isfield(self.(fname),key) || isempty(self.(fname).(key)); self.(fname).(key) = val; end
                %NO%key = "markerEdgeColor"; val = []; if ~isfield(self.(fname).kws,key) || isempty(self.(fname).kws.(key)); self.(fname).kws.(key) = val; end
                %NO%key = "markerFaceColor"; val = []; if ~isfield(self.(fname).kws,key) || isempty(self.(fname).kws.(key)); self.(fname).kws.(key) = val; end
                if self.type.isLine
                    fname = "plot";
                    key = "color"; val = uint8([200 200 200 150]); if ~isfield(self.(fname).kws,key) || isempty(self.(fname).kws.(key)); self.(fname).kws.(key) = val; end
                end
            end

            if self.type.isLine % && self.plot.enabled
                if self.surface.enabled; self.plot.enabled = false; end
                fname = "plot";
                key = "lineWidth"; val = 1; if ~isfield(self.(fname).kws,key) || isempty(self.(fname).kws.(key)); self.(fname).kws.(key) = val; end
            end

            if self.type.isLine % && self.surface.enabled
                fname = "surface";
                keyList = ["faceColor","edgeColor","edgeAlpha","lineStyle","marker"];
                valueList = {"none","flat",0.5,"-","none"};
                for i = 1:length(keyList)
                    key = keyList(i); val = 1; if ~isfield(self.(fname).kws,key) || isempty(self.(fname).kws.(key)); self.(fname).kws.(key) = valueList{i}; end
                end
            end

            if self.colormap.enabled && ~getVecLen(self.colormap.values)
                self.colormap.values = "winter";
                %if self.type.is3d
                %    self.colormap.values = "winter";
                %else
                %    self.colormap.values = "autumn";
                %end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if self.isdryrun; return; end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            cEnabled = (self.type.isLine && self.colormap.enabled && self.surface.enabled) || (self.type.isScatter && self.colormap.enabled);
            if self.type.isLine && ~cEnabled
                self.plot.enabled = true;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% generate figure and axes if needed
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if self.figure.enabled
                figure_kws_cell = convertStruct2Cell(self.figure.kws,{"enabled","singleOptions"});
                %if isfield(self.figure.kws,"singleOptions"); figure_kws_cell = { figure_kws_cell{:} , self.figure.kws.singleOptions{:} }; end
                self.currentFig.figure = figure( figure_kws_cell{:} );
            else
                set(0, "CurrentFigure", gcf);
                self.currentFig.figure = gcf;
                hold on;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% check rows presence
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if getVecLen(self.rows)
                self.rowsindex = self.rows;
            else
                self.rowsindex = 1:1:length(self.dfref{:,1});
            end
            self.ndata = length(self.rowsindex);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% check columns presence
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if getVecLen(self.xcolumns)
                [self.xcolnames, self.xcolindex] = getColNamesIndex(self.dfref.Properties.VariableNames, self.xcolumns);
            else
                self.xcolindex = [];
                self.xcolnames = "Count";
                xdata = self.rowsindex;
            end

            [self.ycolnames, self.ycolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.ycolumns);

            if self.type.is3d
                if getVecLen(self.zcolumns)
                    [self.zcolnames, self.zcolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.zcolumns);
                else
                    self.zcolindex = [];
                    self.zcolnames = "Count";
                    zdata = self.rowsindex;
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% set color data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if self.colormap.enabled
                if getVecLen(self.ccolumns)
                    [self.ccolnames, self.ccolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.ccolumns);
                else
                    self.ccolindex = [];
                    self.ccolnames = "Count";
                    cdata = 1:1:length(self.dfref{self.rowsindex,1});
                end
            else
                self.ccolindex = [];
                self.ccolnames = [];
                cdata = [];
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% check the lengths are consistent
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            xcolindexlen = length(self.xcolindex);
            ycolindexlen = length(self.ycolindex);
            zcolindexlen = length(self.zcolindex);
            ccolindexlen = length(self.ccolindex);
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

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% assign data in case of single column assignments
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if xcolindexlen==1; xdata = self.dfref{self.rowsindex,self.xcolindex}; end
            if ycolindexlen==1; ydata = self.dfref{self.rowsindex,self.ycolindex}; end
            if zcolindexlen==1; zdata = self.dfref{self.rowsindex,self.zcolindex}; end
            if ccolindexlen==1; cdata = self.dfref{self.rowsindex,self.ccolindex}; end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% get keyword arguments
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if self.type.isLine
                plot_kws_cell = convertStruct2Cell(self.plot.kws,{"enabled","singleOptions"});
            end
            if self.type.isScatter
                scatter_kws_cell = convertStruct2Cell(self.scatter.kws,{"enabled","singleOptions","color","cdata","size"});
                if self.scatter.filled; scatter_kws_cell = {"filled", scatter_kws_cell{:}}; end
            end
            if self.type.isLine
                surface_kws_cell = convertStruct2Cell(self.surface.kws,{"enabled","singleOptions","color","size"});
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% add line/scatter plot
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            lglabels = [];

            if ~self.type.is3d
                self.target.currentFig = struct();
                self.target.counter.xline = 0;
                self.target.counter.yline = 0;
                self.target.counter.scatter = 0;
            end

            if (self.type.isScatter && self.scatter.enabled) || (self.type.isLine && (self.surface.enabled || self.plot.enabled))

                isMultiColorScatterPlot = false;
                if ~self.colormap.enabled && self.type.isScatter
                    if getVecLen(self.scatter.color)
                            if all(size(self.scatter.color)==[maxLenColumns,3])
                                isMultiColorScatterPlot = true;
                                scatterMultiColorList = self.scatter.color;
                            elseif all(size(self.scatter.color)==[1,3])
                                cdata = self.scatter.color;
                            else
                                error   ( "The value specified for the 'scatter.color' property of the Scatter Plot object must be either " ...
                                        + "an RGB triplet vector of size [1,3], or a matrix of size [" + string(maxLenColumns) + ",3] for the " ...
                                        + "current set of variables that are selected to plot. It can also be an empty object, in which case, " ...
                                        + "the colors of the objects in the plot will be chosen automatically." ...
                                        );
                            end
                    else
                        scatterMultiColorList = lines(maxLenColumns);
                        isMultiColorScatterPlot = true;
                    end
                end

                if ~self.type.is3d
                    zdata = zeros(self.ndata,1);
                end

                lgxicol = 1;
                lgyicol = 1;
                lgzicol = 1;

                counter.plot = 0;
                counter.plot3 = 0;
                counter.surface = 0;
                counter.scatter = 0;
                counter.scatter3 = 0;

                for icol = 1:maxLenColumns

                    if xcolindexlen>1; lgxicol = icol; xdata = self.dfref{self.rowsindex,self.xcolindex(icol)}; end
                    if ycolindexlen>1; lgyicol = icol; ydata = self.dfref{self.rowsindex,self.ycolindex(icol)}; end
                    if zcolindexlen>1; lgzicol = icol; zdata = self.dfref{self.rowsindex,self.zcolindex(icol)}; end

                    if isMultiColorScatterPlot
                        cdata = scatterMultiColorList(icol,:);
                    else
                        if ccolindexlen>1; cdata = self.dfref{self.rowsindex,self.ccolindex(icol)}; end
                    end

                    if self.legend.enabled % && ~self.type.is3d
                        if xcolindexlen<2 && ycolindexlen>1
                            lglabels = [ lglabels , self.ycolnames(lgyicol) ];
                        elseif xcolindexlen>1 && ycolindexlen<2
                            lglabels = [ lglabels , self.xcolnames(lgxicol) ];
                        else
                            lglabels = [ lglabels , self.xcolnames(lgxicol)+"-"+self.ycolnames(lgyicol) ];
                        end
                        if zcolindexlen>1
                            lglabels(end) = lglabels(end) + "-" + self.zcolnames(lgzicol);
                        end
                    end

                    % add line plot

                    if self.type.isLine

                        if self.plot.enabled
                            if self.type.is3d
                                counter.plot3 = counter.plot3 + 1;
                                self.currentFig.plot3{counter.plot3} = plot3( xdata ...
                                                                            , ydata ...
                                                                            , zdata ...
                                                                            , plot_kws_cell{:} ...
                                                                            );
                            else
                                counter.plot = counter.plot + 1;
                                self.currentFig.plot{counter.plot} = plot   ( xdata ...
                                                                            , ydata ...
                                                                            , plot_kws_cell{:} ...
                                                                            );
                            end
                            hold on;
                        end

                        if self.surface.enabled && self.colormap.enabled
                            counter.surface = counter.surface + 1;
                            self.currentFig.surface{counter.surface} = surface  ( "XData",[xdata(:) xdata(:)] ...
                                                                                , "YData",[ydata(:) ydata(:)] ...
                                                                                , "ZData",[zdata(:) zdata(:)] ...
                                                                                , "CData",[cdata(:) cdata(:)] ...
                                                                                , surface_kws_cell{:} ...
                                                                                );
                            if self.type.is3d; view(3); end
                            hold on;
                        end

                    end

                    if self.type.isScatter && self.scatter.enabled
                        if self.type.is3d
                            counter.scatter3 = counter.scatter3 + 1;
                            self.currentFig.scatter3{counter.scatter3} = scatter3   ( xdata ...
                                                                                    , ydata ...
                                                                                    , zdata ...
                                                                                    , self.scatter.size ...
                                                                                    , cdata ...
                                                                                    , self.scatter.marker ...
                                                                                    , scatter_kws_cell{:} ...
                                                                                    );
                        else
                            counter.scatter = counter.scatter + 1;
                            self.currentFig.scatter{counter.scatter} = scatter  ( xdata ...
                                                                                , ydata ...
                                                                                , self.scatter.size ...
                                                                                , cdata ...
                                                                                , self.scatter.marker ...
                                                                                , scatter_kws_cell{:} ...
                                                                                );
                        end
                        hold on;
                    end

                end % loop plot

                self.currentFig.axes = gca;
                if cEnabled
                    colormap(self.currentFig.axes,self.colormap.values);
                end

            end

            % add axis labels

            if xcolindexlen>1
                self.currentFig.xlabel = xlabel("Variable Values", "Interpreter", "none");
            else
                self.currentFig.xlabel = xlabel(self.xcolnames(1), "Interpreter", "none");
            end

            if ycolindexlen>1
                self.currentFig.ylabel = ylabel("Variable Values", "Interpreter", "none");
            else
                self.currentFig.ylabel = ylabel(self.ycolnames(1), "Interpreter", "none");
            end

            if self.type.is3d
            if zcolindexlen>1
                self.currentFig.zlabel = zlabel("Variable Values", "Interpreter", "none");
            else
                self.currentFig.zlabel = zlabel(self.zcolnames(1), "Interpreter", "none");
            end
            end

            % add line colorbar

            if self.colorbar.enabled && cEnabled
                if isempty(self.colorbar.kws.fontSize) || ~isa(self.colorbar.kws.fontSize,"numeric")
                    self.colorbar.kws.fontSize = self.currentFig.ylabel.FontSize;
                end
                colorbar_kws_cell = convertStruct2Cell(self.colorbar.kws,{"enabled","singleOptions"});
                self.currentFig.colorbar = colorbar(colorbar_kws_cell{:});
                ylabel(self.currentFig.colorbar, self.ccolnames(1), "fontsize", self.colorbar.kws.fontSize, "Interpreter", "none");
            else
                colorbar('off');
                self.currentFig.colorbar = [];
            end

            if self.legend.enabled && (~isfield(self.legend,"labels") || isempty(self.legend.labels)); self.legend.labels = lglabels; end
            self.doBasePlotStuff();

            for fname = ["xlim","ylim","zlim"]
                if ~isempty(self.(fname))
                    limit = get(self.currentFig.axes, fname);
                    for i = 1:2
                        if ~isnan(self.(fname)(i))
                            limit(i) = self.(fname)(i);
                        end
                    end
                    set(self.currentFig.axes, fname, limit);
                end
            end

            %box on; grid on; hold off;
            hold off;

        end % function plot

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % classdef LineScatterPlot
