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
%   GridPlot(dataFrame, plotType)
%
%   This is the GridPlot class for generating instances of square grid plots
%   comprised of histogram/histfit/line/scatter/lineScatter subplots
%   based on a wide variety of ParaMonte's plotting functions.
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
%       columns
%
%           a string or array of strings or cell array of chars representing the
%           names of the columns from the dataFrame to add to the plot.
%
%   Attributes
%   ----------
%
%       ccolumn (standing for color-columns)
%
%           Optional property that determines the column of dataFrame to serve
%           as the color-mapping-values for each line/point element in the line/scatter
%           plot. It can be either a char vector or a string, or the index of the column
%           of interest from the input dataFrame.
%
%           Example usage:
%
%               1.  ccolumn = 7 % column #7 of the dataFrame
%               2.  ccolumn = "SampleLogFunc" % the column of the dataFrame with this name
%
%       colormap
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
%       colorbar
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
%       title
%
%           A MATLAB struct() with the following components:
%
%               enabled
%
%                   A logical value. If `true`, the title will be applied to the plot.
%
%               kws
%
%                   A MATLAB struct() whose components' values are passed to 
%                   MATLAB's title() function. If your desired attribute is 
%                   missing from the fieldnames of title.kws, simply add 
%                   a new field named as the attribute and assign the desired 
%                   value to it.
%
%           Example usage:
%
%               title.enabled = true; % add title
%               title.kws.location = "west";
%
%           WARNING
%
%               Keep in mind that MATLAB keyword arguments are case-INsensitive.
%               Hence, sure you do not add the keyword as multiple different fields.
%               For example, title.kws.fontSize and title.kws.FONTSIZE are the same,
%               and only one of the two will be processed.
%
%       template
%
%           A MATLAB struct() containing the templates of the ParaMonte plot objects 
%           and their properties will will be used to set the properties of similar 
%           types of plots in the subplots of the grid plot.
%
%       plotType
%
%           A MATLAB struct() containing the types of the plots to be plotted in each 
%           section (upper, lower triangles) and the diagonal elements of the GridPlot.
%
%      WARNING: Adding 3d subplots to the GridPlot is not supported anymore as of 
%      WARNING: ParaMonte 2.0.0. Although it is theoretically possible to add 3d subplots 
%      WARNING: to the gridplot (line3, scatter3, linescatter3), there is no practical 
%      WARNING: use for them. Therefore, you should use them at your own risk.
%      WARNING: Perhaps the most meaningful scenario would be when the third
%      WARNING: Z-axis variable is the same column as the ccolumn of the grid-plot.
%      WARNING: If you find the 3D plots useful, or find bugs with the 3d subplots,
%      WARNING: please report it at: https://github.com/cdslaborg/paramonte/issues
%
%   Superclass Attributes
%   ---------------------
%
%       See the documentation for the BasePlot class
%
%   Returns
%   -------
%
%       an object of GridPlot class
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
classdef GridPlot < BasePlot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Access = protected, Hidden)
        %dfref = [];
        %isdryrun = [];
        plotTypeListLen
        plotTypeList
        lowerEnabled
        upperEnabled
        diagEnabled
        diagYTick
        ccolnames
        ccolindex
        colnames
        colindex
    end

    properties (Access = public)

        title
        columns
        ccolumn
        plotType
        %line_kws
        %hist_kws
        %histogram2
        %histfit_kws
        %scatter_kws
        %scatter3_kws
        %lineScatter_kws
        %lineScatter3_kws
        colorbar
        colormap
        template

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function resetInternal(self)

            resetInternal@BasePlot(self);

            trianglePlotList = ["contour", "contourf", "lineScatter"]; %, "line", "scatter", "contour3", "line3", "scatter3", "linescatter3"];
            diagonalPlotList = ["histogram", "histfit"];

            self.plotTypeList = [trianglePlotList, diagonalPlotList];
            self.plotTypeListLen = length(self.plotTypeList);

            self.plotType.upper.enabled = true;
            self.plotType.upper.value = [];
            self.plotType.upper.names = trianglePlotList;

            self.plotType.lower.enabled = true;
            self.plotType.lower.value = [];
            self.plotType.lower.names = trianglePlotList;

            self.plotType.diag.enabled = true;
            self.plotType.diag.value = [];
            self.plotType.diag.names = diagonalPlotList;

            self.axes.main.margin.bottom = 0.07; % 0.14;
            self.axes.main.margin.right = 0.10;
            self.axes.main.margin.left = 0.07; %0.1;
            self.axes.main.margin.top = 0.0;
            self.axes.main.origin.x = 0.02;
            self.axes.main.origin.y = 0.00;
            %self.axes.main.height = 1 - self.axes.main.margin.top;
            %self.axes.main.width = 1 - self.axes.main.margin.right;
            self.axes.main.nrow = length(self.columns);
            self.axes.main.ncol = self.axes.main.nrow;
            self.axes.kws = struct();

            self.axes.subplot.interspace = 0.015; % space between subplots
            %self.axes.subplot.height = (1-self.axes.main.margin.top-self.axes.main.margin.bottom)/self.axes.main.ncol - self.axes.subplot.interspace;
            %self.axes.subplot.width  = (1-self.axes.main.margin.left-self.axes.main.margin.right)/self.axes.main.nrow - self.axes.subplot.interspace;

            self.title = struct();
            self.title.text = "";
            %self.title.subtext = "";
            self.title.enabled = false;
            self.title.kws = struct();
            self.title.kws.fontSize = 13;

            self.colormap = struct();
            self.colormap.enabled = true;
            self.colormap.values = [];

            self.colorbar = struct();
            self.colorbar.enabled = true;
            self.colorbar.kws = struct();
            self.colorbar.kws.fontSize = 12;
            %self.colorbar.origin.x = 1 - self.axes.main.margin.right;
            %self.colorbar.origin.y = self.axes.main.margin.bottom;
            %self.colorbar.height = self.axes.main.nrow * ( self.axes.subplot.height + self.axes.subplot.interspace ) - self.axes.subplot.interspace;
            self.colorbar.width = 0.03;

            self.figure.height = 900;
            %self.figure.width = self.figure.height * ( 1 ...
            %                                         - self.axes.main.margin.bottom ...
            %                                         + self.axes.main.margin.right ...
            %                                         + self.axes.main.margin.left ...
            %                                         - self.axes.main.margin.top ...
            %                                         + self.colorbar.width ...
            %                                         );

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% setup subplot templates
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            templateLocal = struct();

            orangered = getRGB("deep carrot orange");

            % target

            templateLocal.target = struct();
            templateLocal.target.vline = struct();
            templateLocal.target.vline.kws = struct();
            templateLocal.target.hline = struct();
            templateLocal.target.hline.kws = struct();
            templateLocal.target.scatter = struct();
            templateLocal.target.scatter.kws = struct();

            templateLocal.target.vline.kws.lineWidth = 1;
            templateLocal.target.vline.kws.lineStyle = "-";
            templateLocal.target.vline.kws.color = orangered;

            templateLocal.target.hline.kws.lineWidth = 1;
            templateLocal.target.hline.kws.lineStyle = "-";
            templateLocal.target.hline.kws.color = orangered;

            templateLocal.target.scatter.size = 40;
            templateLocal.target.scatter.color = orangered;

            % contour

            templateLocal.contour = struct();
            templateLocal.contour.enabled = [];
            templateLocal.contour.kws = struct();

            % contourf

            templateLocal.contourf = struct();
            templateLocal.contourf.enabled = [];
            templateLocal.contourf.kws = struct();

            % histogram

            templateLocal.histogram = struct();
            templateLocal.histogram.enabled = [];
            templateLocal.histogram.kws = struct();

            % histfit

            templateLocal.histfit = struct();
            templateLocal.histfit.enabled = [];
            templateLocal.histfit.kws = struct();

            % line / scatter / lineScatter

            templateLocal.scatter = struct();
            templateLocal.scatter.enabled = [];
            templateLocal.scatter.kws = struct();

            templateLocal.plot = struct();
            templateLocal.plot.enabled = [];
            templateLocal.plot.kws = struct();

            templateLocal.surface = struct();
            templateLocal.surface.enabled = [];
            templateLocal.surface.kws = struct();

            % figure

            templateLocal.figure = struct();
            templateLocal.figure.kws = struct();
            templateLocal.figure.enabled = false;

            % axes template

            templateLocal.axes = struct();
            templateLocal.axes.kws = struct();

            % legend template

            templateLocal.legend = struct();
            templateLocal.legend.enabled = false;
            templateLocal.legend.kws = struct();

            % colorbar template

            templateLocal.colorbar = struct();
            templateLocal.colorbar.enabled = false;
            templateLocal.colorbar.kws = struct();

            % colormap template

            templateLocal.colormap = struct();
            templateLocal.colormap.enabled = false;
            templateLocal.colormap.values = [];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% setup subplots template
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.template = struct();

            self.template.histfit = struct();
            self.template.histfit.axes = templateLocal.axes;
            self.template.histfit.figure = templateLocal.figure;
            self.template.histfit.histfit = templateLocal.histfit;
            self.template.histfit.legend = templateLocal.legend;

            self.template.histogram = struct();
            self.template.histogram.axes = templateLocal.axes;
            self.template.histogram.figure = templateLocal.figure;
            self.template.histogram.histogram = templateLocal.histogram;
            self.template.histogram.legend = templateLocal.legend;

            self.template.contour = struct();
            self.template.contour.gridSize = 2^9;
            self.template.contour.noiseLevel = 0.001;
            self.template.contour.axes = templateLocal.axes;
            self.template.contour.figure = templateLocal.figure;
            self.template.contour.contour = templateLocal.contour;
            self.template.contour.colorbar = templateLocal.colorbar;
            self.template.contour.colormap = templateLocal.colormap;
            self.template.contour.legend = templateLocal.legend;

            self.template.contourf = struct();
            self.template.contourf.gridSize = 2^9;
            self.template.contourf.noiseLevel = 0.001;
            self.template.contourf.axes = templateLocal.axes;
            self.template.contourf.figure = templateLocal.figure;
            self.template.contourf.contourf = templateLocal.contourf;
            self.template.contourf.colorbar = templateLocal.colorbar;
            self.template.contourf.colormap = templateLocal.colormap;
            self.template.contourf.legend = templateLocal.legend;

            %self.template.line = struct();
            %self.template.line.figure = templateLocal.figure;
            %self.template.line.plot = templateLocal.plot;
            %self.template.line.colorbar = templateLocal.colorbar;
            %self.template.scatter = struct();
            %self.template.scatter.figure = templateLocal.figure;
            %self.template.scatter.scatter = templateLocal.scatter;
            %self.template.scatter.colorbar = templateLocal.colorbar;

            self.template.lineScatter = struct();
            self.template.lineScatter.axes = templateLocal.axes;
            self.template.lineScatter.figure = templateLocal.figure;
            self.template.lineScatter.plot = templateLocal.plot;
            self.template.lineScatter.surface = templateLocal.surface;
            self.template.lineScatter.scatter = templateLocal.scatter;
            self.template.lineScatter.colorbar = templateLocal.colorbar;
            self.template.lineScatter.colormap = templateLocal.colormap;
            self.template.lineScatter.legend = templateLocal.legend;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %self.columns = {};
            self.ccolumn = {};

            %for i = 1:self.plotTypeListLen
            %    self.axes.subplot.(self.plotTypeList{i}) = struct();
            %end

            self.updateLayout();

            self.isdryrun = true;
            self.make();
            self.isdryrun = false;

            % self.target = GridTarget(self.dfref, self.columns, self.currentFig);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = GridPlot(dataFrame, column, resetExternal)
            if nargin<3; resetExternal = []; end
            self = self@BasePlot("grid", dataFrame, resetExternal);
            if nargin<3; self.resetExternal = @self.resetInternal; end
            if nargin<2
                error   ( "GridPlot requires two input arguments: " + newline + newline ...
                        + "    the data Table " + newline ...
                        + "    and a cell array of column names to plot " + newline + newline ...
                        + "The column names is required in order to set the layout of the Grid plot." ...
                        );
            else
                self.columns = column;
            end
            self.resetInternal()
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
                if strcmpi(varargin{1},"update")
                    cmd = "doc self.updateLayout";
                    methodNotFound = false;
                else
                    methodList = ["plot","helpme","addTarget","show","hide","rotateAxesLabels","exportFig"];
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
            %       If the property is a struct();, then its value must be given as a cell array,
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
            %       make("ccolumn",8)
            %       make("colormap","autumn")
            %

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% parse arguments
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.parseArgs(varargin{:});

            % ensure the number of subplots is consistent with the existing layout

            if length(self.columns)~=self.axes.main.nrow || length(self.columns)~=self.axes.main.ncol
                self.currentFig.subplotList = [];
                self.updateLayout();
            end

            self.diagYTick = cell(self.axes.main.nrow,1);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% set plot types to make
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if self.plotType.lower.enabled
                if isempty(self.plotType.lower.value)
                    self.plotType.lower.value = "contour"; %"histogram2";
                else
                    self.plotType.lower.value = string(self.plotType.lower.value);
                end
            end

            if self.plotType.upper.enabled
                if isempty(self.plotType.upper.value)
                    self.plotType.upper.value = "lineScatter";
                else
                    self.plotType.upper.value = string(self.plotType.upper.value);
                end
            end

            if self.plotType.diag.enabled
                if isempty(self.plotType.diag.value)
                    self.plotType.diag.value = "histogram";
                else
                    self.plotType.diag.value = string(self.plotType.diag.value);
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% adjust figure position
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.figure.kws.position =  [ 50  ... origin x
                                        , 100 ... origin y
                                        , self.figure.width  ... width
                                        , self.figure.height ... height
                                        ];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% set what to plot
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if getVecLen(self.ccolumn)>1
                error( "the ccolumn attribute should be set to a single column name from the input data Table.");
            end

            %if self.colormap.enabled && ~getVecLen(self.colormap.values)
            %    self.colormap.values = "winter";
            %end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% check subplot layout
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            columnsIsPresent = getVecLen(self.columns);
            for i = 1:self.plotTypeListLen

                % set template legend properties

                if ~isfield(self.template.(self.plotTypeList(i)).legend,"enabled") || isempty(self.template.(self.plotTypeList(i)).legend.enabled)
                    self.template.(self.plotTypeList(i)).legend.enabled = false;
                end

                % set template rows to plot

                if ~getVecLen(self.rows)
                    self.rows = 1:length(self.dfref{:,1});
                end
                self.template.(self.plotTypeList(i)).rows = self.rows;

                %%%% set figure properties

                fname = "figure";
                self.template.(self.plotTypeList(i)).(fname).enabled = false;
                if ~isfield(self.template.(self.plotTypeList(i)).(fname).kws,"kws") || isempty(self.template.(self.plotTypeList(i)).(fname).kws)
                    self.template.(self.plotTypeList(i)).(fname).kws = struct();
                end

                %%%% set axes properties

                fname = "axes";
                if ~isfield(self.template.(self.plotTypeList(i)).(fname).kws,"kws") || isempty(self.template.(self.plotTypeList(i)).(fname).kws)
                    self.template.(self.plotTypeList(i)).(fname).kws = struct();
                end
                propList = ["box", "xgrid", "ygrid", "fontSize"];
                valueList = {"on", "on", "on", 11};
                for j = 1:length(propList)
                    prop = propList(j);
                    if ~isfield(self.template.(self.plotTypeList(i)).(fname).kws,prop) || isempty(self.template.(self.plotTypeList(i)).(fname).kws.(prop))
                        self.template.(self.plotTypeList(i)).(fname).kws.(prop) = valueList{j};
                    end
                end
                if contains(self.plotTypeList(i),"3")
                    self.template.(self.plotTypeList(i)).(fname).kws.zgrid = "on";
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%% set template plot properties
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if strcmpi(self.plotTypeList(i),"histogram")

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%% histogram
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    fname = "histogram";
                    %if ~isfield(self.template.(self.plotTypeList(i)).(fname),"enabled") || isempty(self.template.(self.plotTypeList(i)).(fname).enabled) % condition creates vicious sticky memory effect
                        self.template.(self.plotTypeList(i)).(fname).enabled = strcmp(self.plotType.diag.value, fname);
                    %end
                    if ~isfield(self.template.(self.plotTypeList(i)).(fname),"kws") || isempty(self.template.(self.plotTypeList(i)).(fname).kws)
                        self.template.(self.plotTypeList(i)).(fname).kws = struct();
                    end
                    propList = ["edgeColor", "facealpha", "faceColor", "edgeColor"];
                    valueList = {"none", 0.6, "auto" "none"};
                    for j = 1:length(propList)
                        prop = propList(j);
                        if ~isfield(self.template.(self.plotTypeList(i)).(fname).kws,prop) || isempty(self.template.(self.plotTypeList(i)).(fname).kws.(prop))
                            self.template.(self.plotTypeList(i)).(fname).kws.(prop) = valueList{j};
                        end
                    end

                elseif strcmpi(self.plotTypeList(i),"histfit")

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%% histfit
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    fname = "histfit";
                    %if ~isfield(self.template.(self.plotTypeList(i)).(fname),"enabled") || isempty(self.template.(self.plotTypeList(i)).(fname).enabled) % condition creates vicious sticky memory effect
                        self.template.(self.plotTypeList(i)).(fname).enabled = strcmp(self.plotType.diag.value, fname);
                    %end
                    if ~isfield(self.template.(self.plotTypeList(i)).(fname),"kws") || isempty(self.template.(self.plotTypeList(i)).(fname).kws)
                        self.template.(self.plotTypeList(i)).(fname).kws = struct();
                    end
                    propList = ["dist", "nbins"];
                    valueList = {"Normal", []};
                    for j = 1:length(propList)
                        prop = propList(j);
                        if ~isfield(self.template.(self.plotTypeList(i)).(fname),prop) || isempty(self.template.(self.plotTypeList(i)).(fname).(prop))
                            self.template.(self.plotTypeList(i)).(fname).(prop) = valueList{j};
                        end
                    end

                else

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%% histogram2
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    if strcmpi(self.plotTypeList(i),"histogram2")

                        fname = "histogram2";
                        %if ~isfield(self.template.(self.plotTypeList(i)).(fname),"enabled") || isempty(self.template.(self.plotTypeList(i)).(fname).enabled) % condition creates vicious sticky memory effect
                            self.template.(self.plotTypeList(i)).(fname).enabled = strcmp(self.plotType.upper.value,fname) || strcmp(self.plotType.lower.value,fname);
                        %end
                        if ~isfield(self.template.(self.plotTypeList(i)).(fname),"kws") || isempty(self.template.(self.plotTypeList(i)).(fname).kws)
                            self.template.(self.plotTypeList(i)).(fname).kws = struct();
                        end
                        propList = ["numbins", "showemptybins", "edgeColor", "faceColor", "displayStyle"];
                        valueList = {[100 100], "off", "none", "flat", "bar3"};
                        for j = 1:length(propList)
                            prop = propList(j);
                            if ~isfield(self.template.(self.plotTypeList(i)).(fname).kws,prop) || isempty(self.template.(self.plotTypeList(i)).(fname).kws.(prop))
                                self.template.(self.plotTypeList(i)).(fname).kws.(prop) = valueList{j};
                            end
                        end

                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%% contour / contourf
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    if contains(lower(self.plotTypeList(i)),"contour")
                        fname = self.plotTypeList(i);
                        %if ~isfield(self.template.(self.plotTypeList(i)).(fname),"enabled") || isempty(self.template.(self.plotTypeList(i)).(fname).enabled) % condition creates vicious sticky memory effect
                            self.template.(self.plotTypeList(i)).(fname).enabled = strcmp(self.plotType.upper.value,fname) || strcmp(self.plotType.lower.value,fname);
                        %end
                        if ~isfield(self.template.(self.plotTypeList(i)).(fname),"kws") || isempty(self.template.(self.plotTypeList(i)).(fname).kws)
                            self.template.(self.plotTypeList(i)).(fname).kws = struct();
                        end
                        propList = ["levels"];
                        valueList = {8};
                        for j = 1:length(propList)
                            prop = propList(j);
                            if ~isfield(self.template.(self.plotTypeList(i)).(fname),prop) || isempty(self.template.(self.plotTypeList(i)).(fname).(prop))
                                self.template.(self.plotTypeList(i)).(fname).(prop) = valueList{j};
                            end
                        end
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%% line / surface
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    if contains(lower(self.plotTypeList(i)),"line")

                        fname = "plot";
                        if ~isfield(self.template.(self.plotTypeList(i)).(fname),"enabled") || isempty(self.template.(self.plotTypeList(i)).(fname).enabled)
                            self.template.(self.plotTypeList(i)).(fname).enabled = true;
                        end
                        if ~isfield(self.template.(self.plotTypeList(i)).(fname),"kws") || isempty(self.template.(self.plotTypeList(i)).(fname).kws)
                            self.template.(self.plotTypeList(i)).(fname).kws = struct();
                        end
                        propList = ["lineWidth", "color"];
                        valueList = {0.75, uint8([200 200 200 100])};
                        for j = 1:length(propList)
                            prop = propList(j);
                            if ~isfield(self.template.(self.plotTypeList(i)).(fname).kws,prop) || isempty(self.template.(self.plotTypeList(i)).(fname).kws.(prop))
                                self.template.(self.plotTypeList(i)).(fname).kws.(prop) = valueList{j};
                            end
                        end

                        fname = "surface";
                        if ~isfield(self.template.(self.plotTypeList(i)).(fname),"enabled") || isempty(self.template.(self.plotTypeList(i)).(fname).enabled)
                            self.template.(self.plotTypeList(i)).(fname).enabled = false;
                        end
                        if ~isfield(self.template.(self.plotTypeList(i)).(fname),"kws") || isempty(self.template.(self.plotTypeList(i)).(fname).kws)
                            self.template.(self.plotTypeList(i)).(fname).kws = struct();
                        end
                        keyList = ["faceColor","edgeColor","edgeAlpha","lineStyle","marker"];
                        valueList = {"none","flat",0.5,"-","none"};
                        for j = 1:length(keyList)
                            key = keyList(j); val = valueList{j};
                            if ~isfield(self.template.(self.plotTypeList(i)).(fname),key) || isempty(self.template.(self.plotTypeList(i)).(fname).(key))
                                self.template.(self.plotTypeList(i)).(fname).(key) = val;
                            end
                        end

                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%% scatter
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    if contains(lower(self.plotTypeList(i)),"scatter")
                        fname = "scatter";
                        if ~isfield(self.template.(self.plotTypeList(i)).(fname),"enabled") || isempty(self.template.(self.plotTypeList(i)).(fname).enabled)
                            self.template.(self.plotTypeList(i)).(fname).enabled = true;
                        end
                        if ~isfield(self.template.(self.plotTypeList(i)).(fname),"kws") || isempty(self.template.(self.plotTypeList(i)).(fname).kws)
                            self.template.(self.plotTypeList(i)).(fname).kws = struct();
                        end
                        propList = ["marker", "size", "color", "filled"];
                        valueList = {"o", 6, [], true};
                        for j = 1:length(propList)
                            prop = propList(j);
                            if ~isfield(self.template.(self.plotTypeList(i)).(fname),prop) || isempty(self.template.(self.plotTypeList(i)).(fname).(prop))
                                self.template.(self.plotTypeList(i)).(fname).(prop) = valueList{j};
                            end
                        end
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%% colormap / colorbar
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    if ~( strcmpi(self.plotTypeList(i),"histogram") || strcmpi(self.plotTypeList(i),"histfit") )

                        %%%% set template properties colormap

                        fname = "colormap";
                        if ~isfield(self.template.(self.plotTypeList(i)),fname) || isempty(self.template.(self.plotTypeList(i)).(fname))
                            self.template.(self.plotTypeList(i)).(fname) = struct();
                        end
                        if ~isfield(self.template.(self.plotTypeList(i)).(fname),"enabled") || isempty(self.template.(self.plotTypeList(i)).(fname).enabled)
                            self.template.(self.plotTypeList(i)).(fname).enabled = self.colormap.enabled;
                        end

                        % NOTE: The colormap of the subplots is set to the main colormap values only if it is non-empty.

                        if ~isfield(self.colormap,"values") || ~getVecLen(self.colormap.values)
                            prop = "values";
                            if ~isfield(self.template.(self.plotTypeList(i)).colormap,prop) || isempty(self.template.(self.plotTypeList(i)).colormap.(prop))
                                if strcmpi(self.plotTypeList(i),"contour")
                                    self.template.(self.plotTypeList(i)).colormap.(prop) = [];
                                elseif strcmpi(self.plotTypeList(i),"contourf")
                                    self.template.(self.plotTypeList(i)).colormap.(prop) = flipud(gray);
                                elseif contains(lower(self.plotTypeList(i)),"line") || contains(lower(self.plotTypeList(i)),"scatter")
                                    self.template.(self.plotTypeList(i)).colormap.(prop) = "winter";
                                end
                            end
                        else
                            self.template.(self.plotTypeList(i)).colormap.values = self.colormap.values;
                        end
                        self.template.(self.plotTypeList(i)).colormap.enabled = self.colormap.enabled;

                        if contains(self.plotTypeList(i),"3") % && ~strcmpi(self.plotTypeList(i),"histogram2")
                            self.template.(self.plotTypeList(i)).zcolumns = self.ccolumn;
                        end

                        self.template.(self.plotTypeList(i)).ccolumns = self.ccolumn;

                        %%%% set template properties colorbar: by default, no colorbar allowed for the subplots

                        prop = "colorbar";
                        if ~isfield(self.template.(self.plotTypeList(i)),prop) || isempty(self.template.(self.plotTypeList(i)).(prop))
                            self.template.(self.plotTypeList(i)).(prop) = struct();
                            self.template.(self.plotTypeList(i)).colorbar.enabled = false;
                        end

                        fname = "colorbar";
                        if ~isfield(self.template.(self.plotTypeList(i)).(fname),"enabled") || isempty(self.template.(self.plotTypeList(i)).(fname).enabled)
                            self.template.(self.plotTypeList(i)).(fname).enabled = false;
                        end
                        if ~isfield(self.template.(self.plotTypeList(i)).(fname),"kws") || isempty(self.template.(self.plotTypeList(i)).(fname).kws)
                            self.template.(self.plotTypeList(i)).(fname).kws = struct();
                        end

                    end
                    %prop = "kws";
                    %defaultEnabled = ~isfield(self.template.(self.plotTypeList(i)).colorbar,prop) || isempty(self.template.(self.plotTypeList(i)).colorbar.(prop));
                    %if defaultEnabled; self.template.(self.plotTypeList(i)).(prop) = struct(); end
                    %prop = "enabled";
                    %defaultEnabled = ~isfield(self.template.(self.plotTypeList(i)).colorbar.kws,prop) || isempty(self.template.(self.plotTypeList(i)).colorbar.kws.(prop));
                    %if defaultEnabled; self.template.(self.plotTypeList(i)).colorbar.kws.(prop) = false; end

                end % set plot properties

            end % loop

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% check columns presence
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if ~columnsIsPresent
                try
                    self.columns = self.dfref.Properties.VariableNames;
                catch % this happens when the dataFrame is empty, for example, when there is no sample in the output files.
                    self.columns = [];
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% set main axes properties
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            key = "position"; val = self.getMainAxisPosition(); if isfield(self.axes.kws,key) && isempty(self.axes.kws.(key)); self.axes.kws.(key) = val; end
            key = "Color"; val = "none"; if isfield(self.axes.kws,key) && isempty(self.axes.kws.(key)); self.axes.kws.(key) = val; end
            key = "Xlim"; val = [0 1]; if isfield(self.axes.kws,key) && isempty(self.axes.kws.(key)); self.axes.kws.(key) = val; end
            key = "Ylim"; val = [0 1]; if isfield(self.axes.kws,key) && isempty(self.axes.kws.(key)); self.axes.kws.(key) = val; end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if self.isdryrun; return; end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% generate figure and axes if needed
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if self.figure.enabled
                figure_kws_cell = convertStruct2Cell(self.figure.kws,{"enabled","singleOptions"});
                %if isfield(self.figure.kws,"singleOptions"); figure_kws_cell = { figure_kws_cell{:} }; end
                self.currentFig.figure = figure( figure_kws_cell{:} );
            else
                self.currentFig.figure = gcf;
            end
            set(0, "CurrentFigure", self.currentFig.figure);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% set columns to plot
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if getVecLen(self.columns)
                [self.colnames, self.colindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.columns);
            else
                error("the component ""columns"" of the GridPlot object appears to be empty.");
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% construct plot handles
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            cornerList = fieldnames(self.plotType);
            cornerListLen = length(cornerList);

            for i = 1:cornerListLen
                corner = string(cornerList{i});
                if ~any( contains(self.plotType.(corner).names, self.plotType.(corner).value) )
                    error   ( newline ...
                            + "Unrecognized plot type requested for the " + corner + " component of ``plotType`` " + newline ...
                            + "component of the GridPlot object. Possible values include: " + newline ...
                            + join(self.plotType.(corner).names,", ") + newline ...
                            );
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% generate axes and subplots
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            axes_kws_cell = convertStruct2Cell(self.axes.kws,{"Parent","parent"});
            self.currentFig.axes = axes ( "Parent", self.currentFig.figure ...
                                        , axes_kws_cell{:} ...
                                        );
            axis(self.currentFig.axes,"off");

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% set up subplot types and classes
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            axisCount = self.axes.main.nrow * self.axes.main.ncol;
            msgSuffix = " out of " + string(axisCount);

            for i = 1:cornerListLen
                corner = string(cornerList{i});
                if self.plotType.(corner).enabled
                    % setup lowerPlotClass, upperPlotClass, diagPlotClass
                    if strcmpi(self.plotType.(corner).value,"line") ||  strcmpi(self.plotType.(corner).value,"scatter") || strcmpi(self.plotType.(corner).value,"lineScatter")
                        thisPlotClass = "LineScatterPlot";
                    else
                        thisPlotClass = "DensityPlot";
                    end
                    if strcmp(corner,"upper")
                        upperPlotClass = thisPlotClass;
                    elseif strcmp(corner,"lower")
                        lowerPlotClass = thisPlotClass;
                    elseif strcmp(corner,"diag")
                        diagPlotClass = thisPlotClass;
                    else
                        error   ( newline ...
                                + "Unrecognized field name for the ``plotType`` component of the GridPlot object." ...
                                + newline ...
                                );
                    end
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% adjust axes limits
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            columnsRange.min = zeros(self.axes.main.nrow,1);
            columnsRange.max = zeros(self.axes.main.nrow,1);
            for i = 1:self.axes.main.nrow
                minvalue = min(self.dfref.(self.colnames(i)));
                maxvalue = max(self.dfref.(self.colnames(i)));
                margin = 0.1 * (maxvalue - minvalue);
                columnsRange.min(i) = minvalue - margin;
                columnsRange.max(i) = maxvalue + margin;
            end

            targetRGB = getRGB("deep carrot orange");

            iaxis = 0;
            self.currentFig.subplotList = cell(self.axes.main.ncol,self.axes.main.nrow);
            for irow = 1:self.axes.main.nrow
                for icol = 1:self.axes.main.ncol

                    jrow = self.axes.main.nrow-irow+1;
                    %jcol = self.axes.main.ncol-icol+1;
                    isDiag = irow==icol && self.plotType.diag.enabled;
                    isLower = irow>icol && self.plotType.lower.enabled;
                    isUpper = irow<icol && self.plotType.upper.enabled;

                    iaxis = iaxis + 1;

                    if isLower || isUpper || isDiag

                        disp( "generating subplot #" + string(iaxis+": ("+irow+","+icol+")") + msgSuffix )

                        currentSubplotAxis = axes   ( "position" ...
                                                    , self.getSubplotAxisPosition(jrow,icol) ...
                                                    );

                        hold on;

                        %%%% generate plot object

                        if isLower
                            currentPlotType = self.plotType.lower.value;
                            currentPlotTypeClass = lowerPlotClass;
                        elseif isUpper
                            currentPlotType = self.plotType.upper.value;
                            currentPlotTypeClass = upperPlotClass;
                        else
                            currentPlotType = self.plotType.diag.value;
                            currentPlotTypeClass = diagPlotClass;
                        end
                        self.currentFig.subplotList{irow,icol} = eval(currentPlotTypeClass+"(currentPlotType,self.dfref)");
                        self.currentFig.subplotList{irow,icol}.xcolumns = self.colnames(icol);
                        if ~self.currentFig.subplotList{irow,icol}.type.is1d
                            self.currentFig.subplotList{irow,icol}.ycolumns = self.colnames(irow);
                        end

                        %%%% set subplot properties

                        self.currentFig.subplotList{irow,icol} = copyComponents ( self.template.(currentPlotType) ...
                                                                                , self.currentFig.subplotList{irow,icol} ...
                                                                                );

                        %fieldList = fieldnames(self.template.(currentPlotType));
                        %%valueList = struct2cell(self.template.(currentPlotType));
                        %for i = 1:length(fieldList)
                        %    if isa(self.template.(currentPlotType).(fieldList{i}),"struct")
                        %        subFieldList = fieldnames( self.template.(currentPlotType).(fieldList{i}) );
                        %        for j = 1:length(subFieldList)
                        %            if isa(self.template.(currentPlotType).(fieldList{i}).(subFieldList{j}),"struct")
                        %                subSubFieldList = fieldnames( self.template.(currentPlotType).(fieldList{i}).(subFieldList{j}) );
                        %                for k = 1:length(subSubFieldList)
                        %                    self.currentFig.subplotList{irow,icol}.(fieldList{i}).(subFieldList{j}).(subSubFieldList{k}) = self.template.(currentPlotType).(fieldList{i}).(subFieldList{j}).(subSubFieldList{k});
                        %                end
                        %            else
                        %                self.currentFig.subplotList{irow,icol}.(fieldList{i}).(subFieldList{j}) = self.template.(currentPlotType).(fieldList{i}).(subFieldList{j});
                        %            end
                        %        end
                        %    else
                        %        self.currentFig.subplotList{irow,icol}.(fieldList{i}) = self.template.(currentPlotType).(fieldList{i});
                        %    end
                        %end
                        %if strcmp(currentPlotType,"histogram2")
                        %    %colormap();
                        %    %self.template.(currentPlotType).colormap.values = flipud(cold());
                        %    %self.template.(currentPlotType).colormap.values = "autumn";
                        %    %colormap(self.template.(currentPlotType).colormap);
                        %end

                        self.currentFig.subplotList{irow,icol}.make();
                        self.currentFig.subplotList{irow,icol}.currentFig.axes = currentSubplotAxis;

                        %%%% set up subplot target

                        if ~self.currentFig.subplotList{irow,icol}.type.is3d
                            self.currentFig.subplotList{irow,icol}.target.scatter.color = targetRGB;
                            self.currentFig.subplotList{irow,icol}.target.vline.color = targetRGB;
                            self.currentFig.subplotList{irow,icol}.target.hline.color = targetRGB;
                            if strcmp(currentPlotType,"histogram") || strcmp(currentPlotType,"histfit")
                                self.currentFig.subplotList{irow,icol}.target.hline.enabled = false;
                                self.currentFig.subplotList{irow,icol}.target.scatter.enabled = false;
                            elseif strcmp(currentPlotType,"contourf")
                                self.currentFig.subplotList{irow,icol}.target.hline.enabled = false;
                                self.currentFig.subplotList{irow,icol}.target.vline.enabled = false;
                                self.currentFig.subplotList{irow,icol}.target.scatter.enabled = false;
                            end
                        end

                        %if irow~=self.axes.main.nrow
                        %    set(gca,'XLabel',[]);
                        %else
                        %    %axisLabelX = get(gca,'XLabel');
                        %    %set(axisLabelX,'rotation',45,'VerticalAlignment','top','HorizontalAlignment','right');
                        %end

                        % adjust the limits
%{
                        if icol==1
                            %axisLabelY = get(gca,'YLabel');
                            %set(axisLabelY,'rotation',45,'VerticalAlignment','middle','HorizontalAlignment','right');
                            if irow==1
                                xlower = self.currentFig.subplotList{1,1}.currentFig.axes.XLim(1);
                                xupper = self.currentFig.subplotList{1,1}.currentFig.axes.XLim(2);
                                xrange = xupper - xlower;
                                ylower = self.currentFig.subplotList{1,1}.currentFig.axes.YLim(1);
                                yupper = self.currentFig.subplotList{1,1}.currentFig.axes.YLim(2);
                                yrange = yupper - ylower;
                                newYTicks = (self.currentFig.subplotList{1,1}.currentFig.axes.XTick - xlower) * yrange / xrange + ylower;
                                yticks(newYTicks);
                                yticklabels(self.currentFig.subplotList{1,1}.currentFig.axes.XTickLabel);
                                ylabel(string(self.currentFig.subplotList{1,1}.xcolumns));
                                %self.currentFig.subplotList{1,1}.currentFig.axes.YTick = self.currentFig.subplotList{1,1}.currentFig.axes.XTick;
                                %self.currentFig.subplotList{1,1}.currentFig.axes.YTickLabel = self.currentFig.subplotList{1,1}.currentFig.axes.XTickLabel;
                            end
                        else
                            set(gca,'YLabel',[]);
                        end

                        %set(self.currentFig.subplotList{irow,icol}.currentFig.axes,"fontSize",self.axes.subplot.fontSize);
                        if icol ~= 1
                            set(self.currentFig.subplotList{irow,icol}.currentFig.axes,"YTickLabels",[]);
                        end
                        if irow ~= self.axes.main.nrow
                            set(self.currentFig.subplotList{irow,icol}.currentFig.axes,"XTickLabels",[]);
                        end

                        self.currentFig.subplotList{irow,icol}.currentFig.axes.XLim = [ columnsRange.min(icol) columnsRange.max(icol) ];
                        if ~( strcmpi(currentPlotType,"histogram") || strcmpi(currentPlotType,"histfit") )
                            self.currentFig.subplotList{irow,icol}.currentFig.axes.YLim = [ columnsRange.min(irow) columnsRange.max(irow) ];
                        end

                        %set(self.currentFig.axisList(iaxis)(iaxis),'box',self.axes.subplot.box);
                        grid(self.currentFig.subplotList{irow,icol}.currentFig.axes,self.axes.subplot.grid);
%}
                        set(gca,'ZLabel',[]);

                        hold off;

                        %%%% set the diagonal YTicks

                        if icol==irow
                            xlower = self.currentFig.subplotList{irow,icol}.currentFig.axes.XLim(1);
                            xupper = self.currentFig.subplotList{irow,icol}.currentFig.axes.XLim(2);
                            xrange = xupper - xlower;
                            ylower = self.currentFig.subplotList{irow,icol}.currentFig.axes.YLim(1);
                            yupper = self.currentFig.subplotList{irow,icol}.currentFig.axes.YLim(2);
                            yrange = yupper - ylower;
                            self.diagYTick{icol}.original.values = self.currentFig.subplotList{irow,icol}.currentFig.axes.YTick;
                            self.diagYTick{icol}.original.labels = self.currentFig.subplotList{irow,icol}.currentFig.axes.YTickLabel;
                            self.diagYTick{icol}.modified.values = (self.currentFig.subplotList{irow,icol}.currentFig.axes.XTick - xlower) * yrange / xrange + ylower;
                            self.diagYTick{icol}.modified.labels = self.currentFig.subplotList{irow,icol}.currentFig.axes.XTickLabel;
                        end

                    else

                        % ATTN: do not remove the extra space after "skipping".
                        disp( "skipping   subplot #" + string(iaxis+": ("+irow+","+icol+")") + msgSuffix );
                        self.currentFig.subplotList{irow,icol} = [];

                    end

                end

            end

            self.hideShowAxesLabels();

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% set up colorbar
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if self.colorbar.enabled && self.colormap.enabled
                self.show("colorbar");
            else
                colorbar('off');
                self.currentFig.colorbar = [];
            end

        end % function plot

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function hide(self,varargin)
            %
            %   Hide the requested section(s) of the gridplot
            %
            %   Parameters
            %   ----------
            %
            %       varargin (optional)
            %
            %           A comma-separated sequence of strings to char-vectors each of which can be one of the following:
            %
            %               "diag"      : corresponding to the diagonal subplots in the grid plot
            %               "upper"     : corresponding to the upper traingle of the grid plot
            %               "lower"     : corresponding to the lower traingle of the grid plot
            %               "colorbar"  : corresponding to the colorbar object of the grid plot
            %
            %           if no input argument is provided, the action is performed for all plot sections listed in the above.
            %
            %   Returns
            %   -------
            %
            %       None.
            %
            %   Example
            %   -------
            %
            %       hide() % hide all parts of the grid plot
            %       hide("diag","colorbar") % hide the diagonal subplots and the colorbar of the grid plot
            %       hide("lower") % hide the lower traingle subplots of the grid plot
            %
            self.hideShow("hide",varargin{:})
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function show(self,varargin)
            %
            %   Show the requested section(s) of the gridplot
            %
            %   Parameters
            %   ----------
            %
            %       varargin (optional)
            %
            %           A comma-separated sequence of strings to char-vectors each of which can be one of the following:
            %
            %               "diag"      : corresponding to the diagonal subplots in the grid plot
            %               "upper"     : corresponding to the upper traingle of the grid plot
            %               "lower"     : corresponding to the lower traingle of the grid plot
            %               "colorbar"  : corresponding to the colorbar object of the grid plot
            %
            %           if no input argument is provided, the action is performed for all plot sections listed in the above.
            %
            %   Returns
            %   -------
            %
            %       None.
            %
            %   Example
            %   -------
            %
            %       show() % show all parts of the grid plot
            %       show("diag","colorbar") % show the diagonal subplots and the colorbar of the grid plot
            %       show("lower") % show the lower traingle subplots of the grid plot
            %
            self.hideShow("show",varargin{:})
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function rotateAxesLabels(self, degreeX, degreeY)
            %
            %   Rotate the axes labels of the subplots of the grid plot
            %
            %   Parameters
            %   ----------
            %
            %       varargin (optional)
            %
            %           No input arguments:
            %
            %               Both X-axes and Y-axes labels are rotated by 45 degrees with respect to the horizontal line
            %
            %           One input arguments:
            %
            %               Both X-axes and Y-axes are rotated by the input numeric value (in degrees) with respect to the horizontal line
            %
            %           Two input arguments (spearated by comma):
            %
            %               The x-axes labels are rotated by the first input numeric value (in degrees) with respect to the horizontal line
            %               The Y-axes labels are rotated by the second input numeric value (in degrees) with respect to the horizontal line
            %
            %   Returns
            %   -------
            %
            %       None.
            %
            %   Example
            %   -------
            %
            %       rotateAxesLabels(20)    % rotate all axes labels by 20 degrees with respect to the horizontal line.
            %       rotateAxesLabels(20,70) % rotate all X-axes/Y-axes labels by 20/70 degrees with respect to the horizontal line.
            %       rotateAxesLabels()      % rotate all axes by 45 degrees with respect to the horizontal line.
            %
            if nargin==1
                degreeX = 45;
                degreeY = 45;
            elseif nargin==2
                degreeY = 45;
            elseif nargin~=3
                error   ( newline ...
                        + "The rotateAxesLabels() method gets at most two input numeric values representing the X-axis and Y-axis label orientations from the horizontal axis, in degrees. " ...
                        + "If only one argument is provided, the input degree value will be used for both X and Y axis label orientations. " + newline ...
                        + "Alternatively, if no argument is passed, then all subplot axis labels will be rotated by the default 45 degrees." ...
                        + newline ...
                        );
            end

            for icol = 1:self.axes.main.ncol
                for irow = 1:self.axes.main.nrow
                    axisLabelX = get(self.currentFig.subplotList{irow,icol}.currentFig.axes,"XLabel"); set(axisLabelX,"rotation",degreeX,"VerticalAlignment","top","HorizontalAlignment","right");
                    axisLabelY = get(self.currentFig.subplotList{irow,icol}.currentFig.axes,"YLabel"); set(axisLabelY,"rotation",degreeY,"VerticalAlignment","middle","HorizontalAlignment","right");
                end
            end

        end % rotateAxesLabels

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function addTarget(self,statistic,varargin)
            %
            %   Add target objects to the subplots of the grid plot.
            %
            %   Parameters
            %   ----------
            %
            %       statistic
            %
            %           A string or char vector with the following possible values:
            %
            %               "mode"      : use the location of the mode of the "SampleLogFunc" column of the dataFrame as the target value
            %               "mean"      : use the mean values of each pair of variables as the target value on each subplot.
            %               "median"    : use the median values of each pair of variables as the target value on each subplot.
            %
            %           If no or an empty input argument is provided, then by default, the target values to be
            %           added to the subplots will be the location of the mode of the data in the "SampleLogFunc"
            %           column of the input dataFrame (dfref).
            %
            %       varargin
            %
            %           pairs of key,value sequences are also possible as input. The possible keys include:
            %
            %               "hline", {}         : the values specified in the cell array, are directly passed to `hline`
            %                                   : component of Target_class object of each subplot.
            %               "vline", {}         : the values specified in the cell array, are directly passed to `vline`
            %                                   : component of Target_class object of each subplot.
            %               "scatter", {}       : the values specified in the cell array, are directly passed to `scatter`
            %                                   : component of Target_class object of each subplot.
            %               "values", values    : A numeric vector of same length as the number of variables used in the gridplot
            %                                   : (that is, the number of row or columns in the grid plot)
            %
            %   Returns
            %   -------
            %
            %       None. However, this method causes side-effects by manipulating the existing attributes
            %       of the `target` components of the `subplot` components of the `axes` component of `layout`.
            %
            %   Example
            %   -------
            %
            %       addTarget() % use the location of the mode of `SampleLogFunc` data-column as the target values on subplots
            %       addTarget("mode") % same as above
            %       addTarget("mean") % use the mean of the vaiables on each axis of each subplot as the target values
            %       addTarget("values",[1.3, 0.2, -0.5, 10]) % set the target values corresponding to each of the four variables in
            %                                                % the grid plot (assuming that there are 4 columns/rows of subplots)
            %
            if nargin==1
                statistic = "mode";
            else
                errorOccurred = false;
                if ~getVecLen(statistic) % isempty
                    key = "statistic";
                    [statistic, keyFound] = getKeyVal(key,varargin{:});
                    if ~keyFound
                        statistic = "mode";
                    end
                elseif isstring(statistic) || ischar(statistic)
                    errorOccurred = true;
                    statistic = string(statistic);
                    for possibleStatistic = ["mode","mean","median"]
                        if strcmpi(possibleStatistic,statistic)
                            errorOccurred = false;
                            break;
                        end
                    end
                else
                    errorOccurred = true;
                end
                if errorOccurred
                    error   ( newline ...
                            + "Invalid value specified for the input argument statistic, which must be either a string or char-vector. The input value is: "  ...
                            + string(statistic) ...
                            + newline ...
                            + "The possible values for the input statistic are: ""mode"", ""mean"", ""median""" ...
                            );
                end
            end

            % parse args

            props = struct();
            for i = 1:length(varargin)
                for kws = ["hline","vline","scatter","values"]
                    [val, keyFound] = getKeyVal(kws,varargin{:});
                    if keyFound
                        props.(kws) = val;
                    end
                end
            end

            % compute the statistics

            if getVecLen(self.columns)
                [self.colnames, self.colindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.columns);
            else
                error("the component ""columns"" of the GridPlot object appears to be empty.");
            end

            columnsLen = length(self.colnames);
            if isfield(props,"values")
                values = props.values;
            else
                values = zeros(columnsLen,1);
                statIsMode = strcmpi(statistic,"mode");
                statIsMean = strcmpi(statistic,"mean");
                statIsMedian = strcmpi(statistic,"median");
                if statIsMode
                    [~,irowMax] = max(self.dfref.SampleLogFunc);
                end
                for i = 1:columnsLen
                    if statIsMode
                        values(i) = self.dfref.(self.colnames(i))(irowMax);
                    elseif statIsMean
                        values(i) = mean(self.dfref.(self.colnames(i)));
                    elseif statIsMedian
                        values(i) = median(self.dfref.(self.colnames(i)));
                    end
                end
            end

            for icol = 1:columnsLen
                for irow = 1:columnsLen

                    set(self.currentFig.figure, "CurrentAxes", self.currentFig.subplotList{irow,icol}.currentFig.axes);

                    hold on;

                    if icol==irow; self.currentFig.subplotList{irow,icol}.target.hline.enabled = false; end
                    self.currentFig.subplotList{irow,icol}.target.values = [ values(icol) values(irow) ];
                    for kws = ["hline","vline","scatter"]
                        if isfield(props,kws)
                            for i = 1:2:length(props.(kws))
                                key = props.(kws){i};
                                val = props.(kws){i+1};
                                self.currentFig.subplotList{irow,icol}.target.(kws).(key) = val;
                            end
                        end
                    end
                    self.currentFig.subplotList{irow,icol}.target.make();

                    hold off;

                end
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function setAxesLabels(self, labels, varargin)
            %
            %   Set the axes labels of the subplots that are currently present in the GridPlot to the input user-provided values.
            %
            %   Parameters
            %   ----------
            %
            %       labels
            %
            %           A cell-array of strings or char-vectors, each of which represents the labels of
            %           the GridPlot's subplot axes, from the left of the GridPlot to the right.
            %           By default, if any element of the labels list is empty, the corresponding column name
            %           from the input dataFrame will be used as the corresponding axis label.
            %
            %           Example usage:
            %
            %               axesLabels = {"Sample Variable 1", [], "This Variable"};
            %
            %       varargin
            %
            %           A sequence of key, value pairs that are directly passed to MATLAB's xlabel() and ylabel() functions.
            %
            %           Example usage:
            %
            %               setAxesLabels( {"Variable_1", [], "Variable_3"} ... the axes labels (the second will be set to default).
            %                            , "fonrtSize", 13 ...
            %                            , "interpreter", "tex" ...
            %                            )
            %
            %   Returns
            %   -------
            %
            %       None. However, this method causes side-effects by manipulating the existing attributes
            %       of the `target` components of the `subplot` components of the `axes` component of `layout`.
            %
            %   Example
            %   -------
            %
            %       setAxesLabels( ["Variable_1", [], "Variable_3"] )% set the axes labels with the default properties.
            %
            %
            %       setAxesLabels( ["Variable_1", [], "Variable_3"] ... the axes labels (the second will be set to default).
            %                    , "fonrtSize", 13 ...
            %                    , "interpreter", "tex" ...
            %                    )
            %
            %       setAxesLabels([], "fonrtSize", 13, "interpreter", "tex") % use default axes labels, but set the properties.
            %
            labelsLen = length(labels);
            if labelsLen==0
                return; % nothing to do, return.
            else
                if ~iscell(labels) || labelsLen>self.axes.main.nrow
                    error   ( "The input labels must be a cell array of length " + string(self.axes.main.nrow) + " comprised of strings or char-vectors." ...
                            + "The input value is " + strjoin(string(labels),", ") ...
                            );
                end
            end
            for icol = 1:labelsLen % self.axes.main.ncol
                for irow = 1:labelsLen % self.axes.main.nrow
                    if ~isempty(self.currentFig.subplotList{irow,icol}) && strcmp(self.currentFig.subplotList{irow,icol}.currentFig.axes.Visible,"on")
                        if isstring(labels{icol}) || ischar(labels{icol})
                            if irow==labelsLen || (irow<labelsLen && ~isempty(self.currentFig.subplotList{irow+1,icol}) && strcmp(self.currentFig.subplotList{irow+1,icol}.currentFig.axes.Visible,"off"))
                                xlabel(self.currentFig.subplotList{irow,icol}.currentFig.axes, labels{icol}, varargin{:});
                            end
                        end
                        if isstring(labels{irow}) || ischar(labels{irow})
                            if icol==1 || (icol>1 && ~isempty(self.currentFig.subplotList{irow,icol-1}) && strcmp(self.currentFig.subplotList{irow,icol-1}.currentFig.axes.Visible,"off"))
                                ylabel(self.currentFig.subplotList{irow,icol}.currentFig.axes, labels{irow}, varargin{:});
                            end
                        end
                    end
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function setAxesLimits(self, limits)
            %
            %   Set the axes x-y limits of the subplots that are currently present in the GridPlot to the input user-provided values.
            %
            %   Parameters
            %   ----------
            %
            %       limits
            %
            %           A cell-array each element of which corresponds to one row of the GridPlot (from the top-left).
            %           Each element of the cell array is a real vector of length two whose values determine the lower
            %           and upper limits of the corresponding variable in the subplots of the GridPlot.
            %           If you wish to leave the limits for some of variables to their default values, set the
            %           corresponding elements of the input cell array to []. By default, if any element of the
            %           input cell is empty, the corresponding limit on the axis will remain unchanged.
            %
            %           Example usage:
            %
            %               limits =
            %               limits = {-10, [], 0};     % leave the limit on the second variable unchanged.
            %               limits = {};               % leave all limits unchanged.
            %
            %   Returns
            %   -------
            %
            %       None. However, this method causes side-effects by manipulating the existing attributes
            %       of the `target` components of the `subplot` components of the `axes` component of `layout`.
            %
            %   Example
            %   -------
            %
            %       setAxesLimits( {[-10,10], [-20, 0]} )
            %       setAxesLimits( {[-10,10], [], [-20, 0]} ) % the limits for the second variable will remain unchanged.
            %
            limitsLen = getVecLen(limits);
            if limitsLen==0
                return; % nothing to do, return.
            else
                if limitsLen>self.axes.main.nrow
                    error   ( "The input limits must be a string vector or a cell array of length " + string(self.axes.main.nrow) + "." ...
                            + "The input value is " + strjoin(string(limits),", ") ...
                            );
                end
            end
            limitsElementLen = zeros(limitsLen,1);
            for icol = 1:limitsLen % self.axes.main.ncol
                limitsElementLen(icol) = getVecLen(limits{icol});
                if limitsElementLen(icol)~=0 && limitsElementLen(icol)~=2
                    error   ( "The element " + string(icol) + " of the input cell array is invalid. " ...
                            + "The limits specified by the input cell must be either vectors of length two, or empty vectors. " ...
                            + "You have entered for the element " + string(icol) + ": " ...
                            + strjoin(string(limits{icol})) ...
                            );
                end
            end
            for icol = 1:limitsLen % self.axes.main.ncol
                for irow = 1:limitsLen % self.axes.main.nrow
                    if strcmp(self.currentFig.subplotList{irow,icol}.currentFig.axes.Visible,"on")
                        if limitsElementLen(icol)==2
                            xlim(self.currentFig.subplotList{irow,icol}.currentFig.axes, limits{icol}(:));
                        end
                        if limitsElementLen(irow)==2 && icol~=irow
                            ylim(self.currentFig.subplotList{irow,icol}.currentFig.axes, limits{irow}(:));
                        end
                    end
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function updateLayout(self)
            %
            %   Update the layout of the grid plot with the new changes.
            %
            %   Parameters
            %   ----------
            %
            %       None.
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
            %       For example, change the left/bottom margin of the main 
            %       axis of the figure to provide room for lengthy variable names.
            %       Then call the updateLayout() method of layout to reflect the changes.
            %
            self.axes.main.height = 1 - self.axes.main.margin.top;
            self.axes.main.width = 1 - self.axes.main.margin.right;
            self.axes.main.nrow = length(self.columns);
            self.axes.main.ncol = self.axes.main.nrow;
            self.axes.subplot.height = (1-self.axes.main.margin.top-self.axes.main.margin.bottom)/self.axes.main.ncol - self.axes.subplot.interspace;
            self.axes.subplot.width  = (1-self.axes.main.margin.left-self.axes.main.margin.right)/self.axes.main.nrow - self.axes.subplot.interspace;
            self.colorbar.origin.x = 1 - self.axes.main.margin.right;
            self.colorbar.origin.y = self.axes.main.margin.bottom;
            self.colorbar.height = self.axes.main.nrow * ( self.axes.subplot.height + self.axes.subplot.interspace ) - self.axes.subplot.interspace;
            self.figure.width    = self.figure.height *   ( 1 ...
                                                                - self.axes.main.margin.bottom ...
                                                                + self.axes.main.margin.right ...
                                                                + self.axes.main.margin.left ...
                                                                - self.axes.main.margin.top ...
                                                                + self.colorbar.width ...
                                                                );
            if isfield(self,"currentFig") && isfield(self.currentFig,"gca") && isgraphics(self.currentFig.axes)
                set(self.currentFig.axes,"position",self.getMainAxisPosition());
            end
            if any(strcmp(fieldnames(self.currentFig),"colorbar")) && isgraphics(self.currentFig.colorbar)
                set(self.currentFig.colorbar,"position",self.getColorbarAxisPosition());
            end
            if isfield(self.currentFig,"subplotList") && ~isempty(self.currentFig.subplotList)
                for irow = 1:self.axes.main.nrow
                    for icol = 1:self.axes.main.ncol
                        jrow = self.axes.main.nrow-irow+1;
                        if isgraphics(self.currentFig.subplotList{irow,icol}.currentFig.axes)
                            set ( self.currentFig.subplotList{irow,icol}.currentFig.axes, "position", self.getSubplotAxisPosition(jrow,icol) );
                        end
                    end
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access=public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function hideShow(self,varargin)

            diagRequested = false;
            lowerRequested = false;
            upperRequested = false;
            colorbarRequested = false;
            action = varargin{1};
            actionIsShow = strcmp(action,"show");
            actionIsHide = strcmp(action,"hide");
            if actionIsShow
                switchValue = "on";
            elseif actionIsHide
                switchValue = "off";
            end
            if nargin==2 % act on all parts of the figure
                diagRequested = true;
                lowerRequested = true;
                upperRequested = true;
                colorbarRequested = true;
            elseif nargin>2
                for i = 2:length(varargin)
                    part = string(varargin{i});
                    if strcmpi(part,"diag")
                        diagRequested  = true;
                    elseif strcmpi(part,"upper")
                        upperRequested = true;
                    elseif strcmpi(part,"lower")
                        lowerRequested = true;
                    elseif strcmpi(part,"colorbar") || strcmpi(part,"cbar")
                        colorbarRequested = true;
                    else
                        error   ( newline ...
                                + "Unrecognized input argument pass to the " + action + "() method: " + newline + newline ...
                                + join(part," ") + newline + newline ...
                                + "The " + action + "() method accepts only a list of comma-separated string values representing the ""part"" of the grid plot to " + action + ". " ...
                                + "Alternatively, if no argument is passed to " + action + "(), then all subplots will be affected." + newline ...
                                );
                    end
                end
            else
                error   ( newline ...
                        + "Internal error occurred. Not enough input arguments to hideShow()" ...
                        + newline ...
                        );
            end

            if colorbarRequested

                if actionIsShow

                    colorbar_kws_cell = convertStruct2Cell(self.colorbar.kws,{"enabled","singleOptions"});

                    if getVecLen(self.colormap.values)
                        colormap(self.currentFig.axes, self.colormap.values);
                    elseif strcmp(self.plotType.upper.value,"lineScatter") || strcmp(self.plotType.lower.value,"lineScatter")
                        if self.template.lineScatter.colormap.enabled && getVecLen(self.template.lineScatter.colormap.values)
                            colormap(self.currentFig.axes, self.template.lineScatter.colormap.values);
                        end
                        %if self.plotType.upper.enabled && strcmp(self.plotType.upper.value,"lineScatter"); indx = [1,2]; end
                        %if self.plotType.lower.enabled && strcmp(self.plotType.lower.value,"lineScatter"); indx = [2,1]; end
                        %colormap(self.currentFig.subplotList{indx(1),indx(2)}, self.colormap.values);
                    else
                        cornerList = ["upper","lower"];
                        for i = 1:length(cornerList)
                            cornerValue = self.plotType.(cornerList(i)).value;
                            if contains(cornerValue,"contour")
                                if self.template.(cornerValue).colormap.enabled && getVecLen(self.template.(cornerValue).colormap.values)
                                    colormap(self.currentFig.axes, self.template.(cornerValue).colormap.values);
                                    break;
                                end
                            end
                        end
                    end

                    % get column names

                    if getVecLen(self.ccolumn)
                        [self.ccolnames, self.ccolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.ccolumn);
                        ccolumnValues = self.dfref.(self.ccolnames);
                    else
                        ccolumnValues = self.rows;
                        self.ccolnames = "Count";
                        disp( newline ...
                            + "the component ""ccolumn"" of the GridPlot object appears to be empty. Using the data counts instead for colormapping..." ...
                            + newline ...
                            );
                    end

                    colorbarLimit = [ min(ccolumnValues) max(ccolumnValues) ]; % Colorbar range
                    caxis(self.currentFig.axes,colorbarLimit);
                    set(self.currentFig.axes,"CLim",colorbarLimit);
                    self.currentFig.colorbar = colorbar(self.currentFig.axes,colorbar_kws_cell{:});

                    colorbarLabel = self.ccolnames;
                    ylabel(self.currentFig.colorbar,colorbarLabel,"fontSize",self.colorbar.kws.fontSize);

                    set ( self.currentFig.colorbar ...
                        , "position" ...
                        , self.getColorbarAxisPosition() ...
                        , "fontSize", self.colorbar.kws.fontSize ...
                        ) ;

                elseif actionIsHide

                    %colorbar("off");
                    set(self.currentFig.colorbar,'visible','off');

                end

            end

            for icol = 1:self.axes.main.ncol
                for irow = 1:self.axes.main.nrow
                    if (icol<irow && lowerRequested) || (icol>irow && upperRequested) || (icol==irow && diagRequested)
                        set(self.currentFig.subplotList{irow,icol}.currentFig.axes,"visible",switchValue);
                        set(get(self.currentFig.subplotList{irow,icol}.currentFig.axes,"children"),"visible",switchValue);
                    end
                end
            end

            if lowerRequested || upperRequested || diagRequested
                self.hideShowAxesLabels();
            end

        end % function hideShow

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function hideShowAxesLabels(self)
            %
            % show or hide axis labels and ticks depending on the presence of the neighbor subplots.
            %
            for icol = 1:self.axes.main.ncol

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%% set the axes labels
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                for irow = 1:self.axes.main.nrow

                    subplotExists = ~isempty(self.currentFig.subplotList{irow,icol});

                    if subplotExists
                        if irow < self.axes.main.nrow
                            if ~isempty(self.currentFig.subplotList{irow+1,icol}) && strcmp(self.currentFig.subplotList{irow+1,icol}.currentFig.axes.Visible,"on")
                                set(self.currentFig.subplotList{irow,icol}.currentFig.axes, "XTickLabel", []);
                                self.currentFig.subplotList{irow,icol}.currentFig.axes.XLabel.String = "";
                            else
                                set(self.currentFig.subplotList{irow,icol}.currentFig.axes, "XTickLabelMode", "auto");
                                self.currentFig.subplotList{irow,icol}.currentFig.axes.XLabel.String = self.colnames(icol);
                            end
                        end

                        if icol > 1
                            if ~isempty(self.currentFig.subplotList{irow,icol-1}) && strcmp(self.currentFig.subplotList{irow,icol-1}.currentFig.axes.Visible,"on")
                                set(self.currentFig.subplotList{irow,icol}.currentFig.axes, "YTickLabel", [])
                                self.currentFig.subplotList{irow,icol}.currentFig.axes.YLabel.String = "";
                            elseif ~isempty(self.currentFig.subplotList{irow,icol})
                                set(self.currentFig.subplotList{irow,icol}.currentFig.axes, "YTickLabelMode", "auto")
                                self.currentFig.subplotList{irow,icol}.currentFig.axes.YLabel.String = self.colnames(irow);
                            end
                        end
                    end

                    if icol~=1 && subplotExists
                        axisLabelY = get(self.currentFig.subplotList{irow,icol}.currentFig.axes,"YLabel");
                        set(axisLabelY,"rotation",45,"VerticalAlignment","middle","HorizontalAlignment","right");
                    end
                    if irow~=self.axes.main.nrow && subplotExists
                        axisLabelX = get(self.currentFig.subplotList{irow,icol}.currentFig.axes,"XLabel");
                        set(axisLabelX,"rotation",45,"VerticalAlignment","top","HorizontalAlignment","right");
                    end

                end

            end % icol

            self.adjustAxesTicks();

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function adjustAxesTicks(self)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% reset the ytick labels of the diagonal subplots to x-axis ticks when there are subplots to their right.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for icol = 1:self.axes.main.ncol

                if ~isempty(self.currentFig.subplotList{icol,icol}) && strcmp(self.currentFig.subplotList{icol,icol}.currentFig.axes.Visible,"on")

                    % do NOT change the order of conditions in the if blocks. short-circuit is at work.
                    if ( icol==self.axes.main.ncol && (isempty(self.currentFig.subplotList{icol,icol-1}) || strcmp(self.currentFig.subplotList{icol,icol-1}.currentFig.axes.Visible,"off")) ) ...
                    || ( icol==1 && ( isempty(self.currentFig.subplotList{icol,icol+1}) || strcmp(self.currentFig.subplotList{icol,icol+1}.currentFig.axes.Visible,"off") ) ) ...
                    || ...
                    ( ...
                        icol~=1 && icol~=self.axes.main.ncol ...
                        && ...
                        ( isempty(self.currentFig.subplotList{icol,icol+1}) || strcmp(self.currentFig.subplotList{icol,icol+1}.currentFig.axes.Visible,"off") ) ...
                        && ...
                        ( isempty(self.currentFig.subplotList{icol,icol-1}) || strcmp(self.currentFig.subplotList{icol,icol-1}.currentFig.axes.Visible,"off") ) ...
                    )

                        % let the histograms have their own y-ranges

                        set(self.currentFig.subplotList{icol,icol}.currentFig.axes, "YTickMode", "auto", "YTickLabelMode", "auto");

                    elseif icol<self.axes.main.ncol

                        subplotCurrentFig = self.currentFig.subplotList{icol,icol}.currentFig;

                        % set the y-range and ticklabels according to the needs of the neighbor plots

                        if icol==1 || isempty(self.currentFig.subplotList{icol,icol-1}) || strcmp(self.currentFig.subplotList{icol,icol-1}.currentFig.axes.Visible,"off")

                            if ~isempty(self.currentFig.subplotList{end,icol})
                                neighborCurrentFig = self.currentFig.subplotList{end,icol}.currentFig;
                            else
                                neighborCurrentFig = subplotCurrentFig;
                            end
%{
                            xlower = subplotCurrentFig.axes.XLim(1);
                            xupper = subplotCurrentFig.axes.XLim(2);
                            xrange = xupper - xlower;

                            ylower = subplotCurrentFig.axes.YLim(1);
                            yupper = subplotCurrentFig.axes.YLim(2);
                            yrange = yupper - ylower;

                            newYTicks = (subplotCurrentFig.axes.XTick - xlower) * yrange / xrange + ylower;

                            yticks(subplotCurrentFig.axes, newYTicks);

                            %%%% set YTickLabels

                            retryEnabled = false;
                            if ~isfield(neighborCurrentFig,"axes") || isempty(neighborCurrentFig.axes.XTickLabel)
                                retryEnabled = true;
                            else
                                try
                                    yticklabels(subplotCurrentFig.axes, neighborCurrentFig.axes.XTickLabel);
                                catch
                                    retryEnabled = true;
                                end
                            end
                            if retryEnabled
                                yticklabels(subplotCurrentFig.axes, subplotCurrentFig.axes.XTickLabel);
                            end
%}
                            yticks(subplotCurrentFig.axes, self.diagYTick{icol}.modified.values);
                            yticklabels(subplotCurrentFig.axes, self.diagYTick{icol}.modified.labels);

                            %%%% set Y-axis label

                            retryEnabled = false;
                            if ~isfield(neighborCurrentFig,"xlabel") || isempty(neighborCurrentFig.xlabel.String)
                                retryEnabled = true;
                            else
                                try
                                    ylabel(subplotCurrentFig.axes, neighborCurrentFig.xlabel.String);
                                catch
                                    retryEnabled = true;
                                end
                            end
                            if retryEnabled
                                retryEnabled = false;
                                if isempty(subplotCurrentFig.xlabel.String)
                                    retryEnabled = true;
                                else
                                    try
                                        ylabel(subplotCurrentFig.axes, string(subplotCurrentFig.xlabel.String));
                                    catch
                                        retryEnabled = true;
                                    end
                                end
                                if retryEnabled
                                    ylabel(subplotCurrentFig.axes, string(self.currentFig.subplotList{icol,icol}.xcolumns));
                                end
                            end

                        end

                    end

                end

            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function position = getColorbarAxisPosition(self)
            position =  [ self.colorbar.origin.x ...
                        , self.colorbar.origin.y ...
                        , self.colorbar.width ...
                        , self.colorbar.height ...
                        ];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function position = getMainAxisPosition(self)
            position =  [ self.axes.main.origin.x ...
                        , self.axes.main.origin.y ...
                        , self.axes.main.width ...
                        , self.axes.main.height ...
                        ];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function position = getSubplotAxisPosition(self,irow,icol)
            position =  [ (icol-1)*(self.axes.subplot.interspace + self.axes.subplot.width)   + self.axes.main.margin.left ...
                        , (irow-1)*(self.axes.subplot.interspace + self.axes.subplot.height)  + self.axes.main.margin.bottom ...
                        , self.axes.subplot.width ...
                        , self.axes.subplot.height ...
                        ];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % classdef GridPlot
