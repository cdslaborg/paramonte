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
%   we ask you to acknowledge the use of the ParaMonte library
%   in your work (education/research/industry/development/...)
%   by citing the ParaMonte library as described on this page:
%
%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%
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
%           optional property that determines the column of dataFrame to serve
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
%           A string or any other value that the colormap function of MATLAB accepts as input.
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
%       layout
%
%           A MATLAB struct() containing an extensive amount of information about the layout
%           and the overall design of the grid plot, with the following components:
%
%               - layout.colorbar   :   the layout and design of the colorbar of the plot
%               - layout.subplot    :   the layout and design of the subplots of the plot
%               - layout.axes       :   the layout and design of the gridplot's axes and the subplots' axes
%               - layout.update()   :   a method that reflects the new layout changes into the grid-plot, when called.
%               - layout.plotType   :   a struct() with components: diag, lower, upper, each of which contains
%                                       the type of plot to be added to the corresponding section of the grid plot.
%                                       possible values for each section includes:
%
%                                           diag: histogram, histfit
%                                           lower: histogram2, contourf, contour, line, scatter, linescatter, line3, scatter3, linescatter3
%                                           upper: histogram2, contourf, contour, line, scatter, linescatter, line3, scatter3, linescatter3
%
%                                       WARNING: Although it is possible to add 3d subplots to the gridplot
%                                       WARNING: (line3, scatter3, linescatter3), there is no practical use
%                                       WARNING: for them. Therefore, you should use them at your own risk.
%                                       WARNING: Perhaps the most meaningful scenario would be when the third
%                                       WARNING: Z-axis variable is the same column as the ccolumn of the grid-plot.
%                                       WARNING: If you find the 3D plots useful, or find bugs with the 3d subplots,
%                                       WARNING: please report it at: https://github.com/cdslaborg/paramonte/issues
%
%   Superclass Attributes
%   ----------------------
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
        plotTypeList =  [ "line" ...
                        , "histfit" ...
                        , "histogram" ...
                        , "histogram2" ...
                        , "contourf" ...
                        , "contour" ...
                        , "scatter" ...
                        , "scatter3" ...
                        , "lineScatter" ...
                        , "lineScatter3" ...
                        ];
        plotTypeListLen
        lowerEnabled
        upperEnabled
        diagEnabled
        colnames
    end

    properties (Access = public)

        columns
        ccolumn
        %line_kws
        %hist_kws
        %histogram2_kws
        %histfit_kws
        %scatter_kws
        %scatter3_kws
        %lineScatter_kws
        %lineScatter3_kws
        colorbar_kws
        colormap
        layout

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function reset(self)

            reset@BasePlot(self);
            self.plotTypeListLen = length(self.plotTypeList);

            %self.columns = {};
            self.ccolumn = {};
            for i = 1:length(self.plotTypeList)
                self.layout.subplot.(self.plotTypeList{i}) = struct();
            end
            self.colormap = {};

            self.colorbar_kws = struct();
            self.colorbar_kws.enabled = true;
            self.colorbar_kws.fontsize = [];

            self.layout.plotType.upper = "lineScatter";
            self.layout.plotType.lower = "contour"; % "contourf"; % "histogram2";
            self.layout.plotType.diag = "histogram";
            self.layout.update = @self.updateLayout;

            self.layout.axis.main.margin.bottom = 0.07; % 0.14;
            self.layout.axis.main.margin.right = 0.10;
            self.layout.axis.main.margin.left = 0.07; %0.1;
            self.layout.axis.main.margin.top = 0.0;
            self.layout.axis.main.origin.x = 0.02;
            self.layout.axis.main.origin.y = 0.00;
            %self.layout.axis.main.height = 1 - self.layout.axis.main.margin.top;
            %self.layout.axis.main.width = 1 - self.layout.axis.main.margin.right;
            self.layout.axis.main.title.content = [];
            self.layout.axis.main.title.fontsize = 13;
            self.layout.axis.main.nrow = length(self.columns);
            self.layout.axis.main.ncol = self.layout.axis.main.nrow;

            self.layout.axis.subplot.interspace = 0.015; % space between subplots
            %self.layout.axis.subplot.height = (1-self.layout.axis.main.margin.top-self.layout.axis.main.margin.bottom)/self.layout.axis.main.ncol - self.layout.axis.subplot.interspace;
            %self.layout.axis.subplot.width  = (1-self.layout.axis.main.margin.left-self.layout.axis.main.margin.right)/self.layout.axis.main.nrow - self.layout.axis.subplot.interspace;
            self.layout.axis.subplot.box = "on";
            self.layout.axis.subplot.grid = "on";
            self.layout.axis.subplot.fontsize = 11;

            %self.layout.colorbar.origin.x = 1 - self.layout.axis.main.margin.right;
            %self.layout.colorbar.origin.y = self.layout.axis.main.margin.bottom;
            %self.layout.colorbar.height = self.layout.axis.main.nrow * ( self.layout.axis.subplot.height + self.layout.axis.subplot.interspace ) - self.layout.axis.subplot.interspace;
            self.layout.colorbar.width = 0.03;
            self.layout.colorbar.fontsize = 12;

            self.layout.figHeight   = 900;
            %self.layout.figWidth    = self.layout.figHeight *   ( 1 ...
            %                                                    - self.layout.axis.main.margin.bottom ...
            %                                                    + self.layout.axis.main.margin.right ...
            %                                                    + self.layout.axis.main.margin.left ...
            %                                                    - self.layout.axis.main.margin.top ...
            %                                                    + self.layout.colorbar.width ...
            %                                                    );

            for i = 1:self.plotTypeListLen
                self.layout.axis.subplot.(self.plotTypeList{i}) = struct();
            end

            self.updateLayout();

            self.isdryrun = true;
            self.plot();
            self.isdryrun = false;

            % self.target = GridTarget(self.dfref, self.columns, self.currentFig);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = GridPlot(varargin) % expected input arguments: dataFrame, dataFrame column names to plot

            self = self@BasePlot(varargin{1});
            try
                self.columns = varargin{2};
            catch
                error   ( "GridPlot requires two input arguments: " + newline + newline ...
                        + "    the data Table " + newline ...
                        + "    and a cell array of column names to plot " + newline + newline ...
                        + "The column names is required in order to set the layout of the Grid plot." ...
                        );
            end
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
            %       plot("ccolumn",8)
            %       plot("colormap","autumn")
            %       plot( "gcf_kws", {"enabled",true,"color","none"} )
            %

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% parse arguments
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %parseArgs(self,varargin{:});
            self.parseArgs(varargin{:});

            % ensure the number of subplots is consistent with the existing layout

            if length(self.columns)~=self.layout.axis.main.nrow || length(self.columns)~=self.layout.axis.main.ncol
                self.currentFig.subplotList = [];
                self.updateLayout();
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% adjust figure position
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.gcf_kws.position = [ 50                    ... origin x
                                    , 100                   ... origin y
                                    , self.layout.figWidth  ... width
                                    , self.layout.figHeight ... height
                                    ];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% set what to plot
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if getVecLen(self.ccolumn)>1
                error( "the ccolumn attribute should be set to a single column name from the input data Table.");
            end

            cEnabled = getVecLen(self.ccolumn) || isa(self.ccolumn,"cell");

            if cEnabled && ~getVecLen(self.colormap)
                self.colormap = "winter";
            end

            % check subplot data columns to plot

            rowsIsPresent = getVecLen(self.rows);
            columnsIsPresent = getVecLen(self.columns);
            ccolumnIsPresent = getVecLen(self.ccolumn);
            for i = 1:self.plotTypeListLen

                plotType_kws = self.plotTypeList(i) + "_kws";

                % check rows presence

                if rowsIsPresent
                    self.layout.subplot.(self.plotTypeList(i)).rows = self.rows;
                else
                    self.layout.subplot.(self.plotTypeList(i)).rows = {};
                end

                %self.layout.subplot.(self.plotTypeList(i)).xcolumns = {};

                % set color properties

                if ~( strcmpi(self.plotTypeList(i),"histogram") || strcmpi(self.plotTypeList(i),"histfit") )

                    %self.layout.subplot.(self.plotTypeList(i)).ycolumns = {};

                    if strcmpi(self.plotTypeList(i),"histogram2")
                        if ~isfield(self.layout.subplot.(self.plotTypeList(i)),"histogram2_kws") || isempty(self.layout.subplot.(self.plotTypeList(i)).histogram2_kws)
                            self.layout.subplot.(self.plotTypeList(i)).histogram2_kws = struct();
                        end
                        keyList = ["numbins", "showemptybins"];
                        valueList = {[100 100], "off"};
                        for j = 1:length(keyList)
                            prop = keyList(j);
                            defaultEnabled = ~isfield(self.layout.subplot.(self.plotTypeList(i)).histogram2_kws,prop) || isempty(self.layout.subplot.(self.plotTypeList(i)).histogram2_kws.(prop));
                            if defaultEnabled; self.layout.subplot.(self.plotTypeList(i)).histogram2_kws.(prop) = valueList{j}; end
                        end
                        prop = "colormap";
                        defaultEnabled = ~isfield(self.layout.subplot.(self.plotTypeList(i)),prop) || isempty(self.layout.subplot.(self.plotTypeList(i)).(prop));
                        if defaultEnabled; self.layout.subplot.(self.plotTypeList(i)).(prop) = {}; end
                    end

                    if strcmpi(self.plotTypeList(i),"contourf")
                        if ~isfield(self.layout.subplot.(self.plotTypeList(i)),"contourf_kws") || isempty(self.layout.subplot.(self.plotTypeList(i)).contourf_kws)
                            self.layout.subplot.(self.plotTypeList(i)).contourf_kws = struct();
                        end
                        prop = "colormap";
                        defaultEnabled = ~isfield(self.layout.subplot.(self.plotTypeList(i)),prop) || isempty(self.layout.subplot.(self.plotTypeList(i)).(prop));
                        if defaultEnabled; self.layout.subplot.(self.plotTypeList(i)).(prop) = gray; end
                    end

                    if strcmpi(self.plotTypeList(i),"contour")
                        if ~isfield(self.layout.subplot.(self.plotTypeList(i)),"contour_kws") || isempty(self.layout.subplot.(self.plotTypeList(i)).contourf_kws)
                            self.layout.subplot.(self.plotTypeList(i)).contourf_kws = struct();
                        end
                        prop = "colormap";
                        defaultEnabled = ~isfield(self.layout.subplot.(self.plotTypeList(i)),prop) || isempty(self.layout.subplot.(self.plotTypeList(i)).(prop));
                        if defaultEnabled; self.layout.subplot.(self.plotTypeList(i)).(prop) = {}; end
                    end

                    if ccolumnIsPresent
                        if contains(self.plotTypeList(i),"3")
                            self.layout.subplot.(self.plotTypeList(i)).zcolumns = self.ccolumn;
                        end
                        if ~(strcmpi(self.plotTypeList(i),"histogram2") || strcmpi(self.plotTypeList(i),"contourf") || strcmpi(self.plotTypeList(i),"contour"))
                            self.layout.subplot.(self.plotTypeList(i)).ccolumns = self.ccolumn;
                            if ~isfield(self.layout.subplot.(self.plotTypeList(i)),"plot_kws") || isempty(self.layout.subplot.(self.plotTypeList(i)).plot_kws)
                                self.layout.subplot.(self.plotTypeList(i)).plot_kws = struct();
                            end
                            keyList = ["linewidth", "color", "enabled"];
                            valueList = {0.75, uint8([200 200 200 100]), true};
                            for j = 1:length(keyList)
                                prop = keyList(j);
                                defaultEnabled = ~isfield(self.layout.subplot.(self.plotTypeList(i)).plot_kws,prop) || isempty(self.layout.subplot.(self.plotTypeList(i)).plot_kws.(prop));
                                if defaultEnabled; self.layout.subplot.(self.plotTypeList(i)).plot_kws.(prop) = valueList{j}; end
                            end
                            if ~isfield(self.layout.subplot.(self.plotTypeList(i)),"surface_kws") || isempty(self.layout.subplot.(self.plotTypeList(i)).surface_kws)
                                self.layout.subplot.(self.plotTypeList(i)).surface_kws = struct();
                            end
                            if ~isfield(self.layout.subplot.(self.plotTypeList(i)).surface_kws,"enabled") || isempty(self.layout.subplot.(self.plotTypeList(i)).surface_kws.enabled)
                                self.layout.subplot.(self.plotTypeList(i)).surface_kws.enabled = false;
                            end
                            prop = "colormap";
                            defaultEnabled = ~isfield(self.layout.subplot.(self.plotTypeList(i)),prop) || isempty(self.layout.subplot.(self.plotTypeList(i)).(prop));
                            if defaultEnabled; self.layout.subplot.(self.plotTypeList(i)).(prop) = self.colormap; end
                        end
                        prop = "colorbar_kws";
                        defaultEnabled = ~isfield(self.layout.subplot.(self.plotTypeList(i)),prop) || isempty(self.layout.subplot.(self.plotTypeList(i)).(prop));
                        if defaultEnabled; self.layout.subplot.(self.plotTypeList(i)).(prop) = struct(); end
                        prop = "enabled";
                        defaultEnabled = ~isfield(self.layout.subplot.(self.plotTypeList(i)).colorbar_kws,prop) || isempty(self.layout.subplot.(self.plotTypeList(i)).colorbar_kws.(prop));
                        if defaultEnabled; self.layout.subplot.(self.plotTypeList(i)).colorbar_kws.(prop) = false; end
                    end

                end

                % check gcf_kws/box/grid presence

                self.layout.subplot.(self.plotTypeList(i)).legend_kws.enabled = false;
                self.layout.subplot.(self.plotTypeList(i)).gcf_kws.enabled = false;
                self.layout.subplot.(self.plotTypeList(i)).gca_kws.box = self.layout.axis.subplot.box;
                self.layout.subplot.(self.plotTypeList(i)).gca_kws.fontsize = self.layout.axis.subplot.fontsize;

            end

            % check columns presence

            if ~columnsIsPresent
                try
                    self.columns = self.dfref.Properties.VariableNames;
                catch
                    self.columns = [];
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
                self.currentFig.gcf = gcf;
            end
            set(0, "CurrentFigure", self.currentFig.gcf);

            % set columns to plot

            if getVecLen(self.columns)
                [self.colnames, colindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.columns);
            else
                error("the component ""columns"" of the GridPlot object appears to be empty.");
            end

            % construct plot handles

            self.lowerEnabled = any(strcmp(self.plotTypeList,self.layout.plotType.lower));
            self.upperEnabled = any(strcmp(self.plotTypeList,self.layout.plotType.upper));
            self.diagEnabled = any(strcmp(self.plotTypeList,self.layout.plotType.diag));

            if self.lowerEnabled
                if isempty(self.layout.plotType.lower)
                    self.layout.plotType.lower = "contourf"; %"histogram2";
                else
                    self.layout.plotType.lower = string(self.layout.plotType.lower);
                end
            end

            if self.upperEnabled
                if isempty(self.layout.plotType.upper)
                    self.layout.plotType.upper = "lineScatter";
                else
                    self.layout.plotType.upper = string(self.layout.plotType.upper);
                end
            end

            if self.diagEnabled
                if isempty(self.layout.plotType.diag)
                    self.layout.plotType.diag = "histogram";
                else
                    self.layout.plotType.diag = string(self.layout.plotType.diag);
                end
            end

            % generate axes and subplots

            self.currentFig.gca = axes  ( "position" ...
                                        , self.getMainAxisPosition()  ...
                                        ..., "Position", mainPosition ...
                                        , "Xlim", [0 1], "Ylim", [0 1] ...
                                        , "Color", "none" ...
                                        , "Parent", self.currentFig.gcf ...
                                        );
            axis(self.currentFig.gca,"off");

            axisCount = self.layout.axis.main.nrow * self.layout.axis.main.ncol;
            msgSuffix = " out of " + string(axisCount);

            if self.lowerEnabled
                if contains(self.layout.plotType.lower,"line") || contains(self.layout.plotType.lower,"scatter")
                    lowerPlotClass = "LineScatterPlot";
                else
                    lowerPlotClass = "DensityPlot";
                end
            end
            if self.upperEnabled
                if contains(self.layout.plotType.upper,"line") || contains(self.layout.plotType.upper,"scatter")
                    upperPlotClass = "LineScatterPlot";
                else
                    upperPlotClass = "DensityPlot";
                end
            end
            if self.diagEnabled
                diagPlotClass = "DensityPlot";
            end

            % adjust limits

            columns.min = zeros(self.layout.axis.main.nrow,1);
            columns.max = zeros(self.layout.axis.main.nrow,1);
            for i = 1:self.layout.axis.main.nrow
                minvalue = min(self.dfref.(self.colnames(i)));
                maxvalue = max(self.dfref.(self.colnames(i)));
                margin = 0.1 * (maxvalue - minvalue);
                columns.min(i) = minvalue - margin;
                columns.max(i) = maxvalue + margin;
            end

            targetRGB = getRGB("deep carrot orange");

            iaxis = 0;
            self.currentFig.subplotList = cell(self.layout.axis.main.ncol,self.layout.axis.main.nrow);
            for irow = 1:self.layout.axis.main.nrow
                for icol = 1:self.layout.axis.main.ncol

                    jrow = self.layout.axis.main.nrow-irow+1;
                    jcol = self.layout.axis.main.ncol-icol+1;
                    isLower = irow>icol && self.lowerEnabled;
                    isUpper = irow<icol && self.upperEnabled;
                    isDiag = irow==icol && self.diagEnabled;

                    iaxis = iaxis + 1;

                    if isLower || isUpper || isDiag

                        disp( "generating subplot #" + string(iaxis+": ("+irow+","+icol+")") + msgSuffix )

                        currentSubplotAxis = axes   ( "position" ...
                                                    , self.getSubplotAxisPosition(jrow,icol) ...
                                                    );

                        hold on;

                        xcolumns = self.colnames(icol);
                        ycolumns = self.colnames(irow);

                        if isLower
                            currentPlotType = self.layout.plotType.lower;
                            currentPlotTypeClass = lowerPlotClass;
                            self.currentFig.subplotList{irow,icol} = eval(currentPlotTypeClass+"(self.dfref,currentPlotType)");
                            self.currentFig.subplotList{irow,icol}.xcolumns = xcolumns;
                            self.currentFig.subplotList{irow,icol}.ycolumns = ycolumns;
                        elseif isUpper
                            currentPlotType = self.layout.plotType.upper;
                            currentPlotTypeClass = upperPlotClass;
                            self.currentFig.subplotList{irow,icol} = eval(currentPlotTypeClass+"(self.dfref,currentPlotType)");
                            self.currentFig.subplotList{irow,icol}.xcolumns = xcolumns;
                            self.currentFig.subplotList{irow,icol}.ycolumns = ycolumns;
                        else
                            currentPlotType = self.layout.plotType.diag;
                            currentPlotTypeClass = diagPlotClass;
                            self.currentFig.subplotList{irow,icol} = eval(currentPlotTypeClass+"(self.dfref,currentPlotType)");
                            self.currentFig.subplotList{irow,icol}.xcolumns = xcolumns;
                        end

                        plotType_kws = currentPlotType + "_kws";
                        fieldList = fieldnames(self.layout.subplot.(currentPlotType));
                        %valueList = struct2cell(self.layout.subplot.(currentPlotType));
                        for i = 1:length(fieldList)
                            if isa(self.layout.subplot.(currentPlotType).(fieldList{i}),"struct")
                                fieldSubList = fieldnames(self.layout.subplot.(currentPlotType).(fieldList{i}));
                                for j = 1:length(fieldSubList)
                                    self.currentFig.subplotList{irow,icol}.(fieldList{i}).(fieldSubList{j}) = self.layout.subplot.(currentPlotType).(fieldList{i}).(fieldSubList{j});
                                end
                            else
                                self.currentFig.subplotList{irow,icol}.(fieldList{i}) = self.layout.subplot.(currentPlotType).(fieldList{i});
                            end
                        end

                        if strcmp(currentPlotType,"histogram2")
                            %colormap();
                            %self.layout.subplot.(currentPlotType).colormap = flipud(cold());
                            %self.layout.subplot.(currentPlotType).colormap = "autumn";
                            %colormap(self.layout.subplot.(currentPlotType).colormap);
                        end

                        self.currentFig.subplotList{irow,icol}.plot()
                        self.currentFig.subplotList{irow,icol}.currentFig.gca = currentSubplotAxis;

                        for fname = ["hline_kws","vline_kws","scatter_kws"]
                            self.currentFig.subplotList{irow,icol}.target.(fname).color = targetRGB;
                        end
                        if strcmp(currentPlotType,"histogram") || strcmp(currentPlotType,"histfit")
                            self.currentFig.subplotList{irow,icol}.target.hline_kws.enabled = false;
                            self.currentFig.subplotList{irow,icol}.target.scatter_kws.enabled = false;
                        elseif strcmp(currentPlotType,"histogram2") || strcmp(currentPlotType,"contourf")
                            self.currentFig.subplotList{irow,icol}.target.hline_kws.enabled = false;
                            self.currentFig.subplotList{irow,icol}.target.vline_kws.enabled = false;
                            self.currentFig.subplotList{irow,icol}.target.scatter_kws.enabled = false;
                        end

                        if irow==self.layout.axis.main.nrow
                            %axisLabelX = get(gca,'XLabel');
                            %set(axisLabelX,'rotation',45,'VerticalAlignment','top','HorizontalAlignment','right');
                        else
                            set(gca,'XLabel',[]);
                        end
                        if icol==1
                            %axisLabelY = get(gca,'YLabel');
                            %set(axisLabelY,'rotation',45,'VerticalAlignment','middle','HorizontalAlignment','right');
                            if irow==1
                                xlower = self.currentFig.subplotList{1,1}.currentFig.gca.XLim(1);
                                xupper = self.currentFig.subplotList{1,1}.currentFig.gca.XLim(2);
                                xrange = xupper - xlower;
                                ylower = self.currentFig.subplotList{1,1}.currentFig.gca.YLim(1);
                                yupper = self.currentFig.subplotList{1,1}.currentFig.gca.YLim(2);
                                yrange = yupper - ylower;
                                newYTicks = (self.currentFig.subplotList{1,1}.currentFig.gca.XTick - xlower) * yrange / xrange + ylower;
                                yticks(newYTicks);
                                yticklabels(self.currentFig.subplotList{1,1}.currentFig.gca.XTickLabel);
                                ylabel(string(self.currentFig.subplotList{1,1}.xcolumns));
                                %self.currentFig.subplotList{1,1}.currentFig.gca.YTick = self.currentFig.subplotList{1,1}.currentFig.gca.XTick;
                                %self.currentFig.subplotList{1,1}.currentFig.gca.YTickLabel = self.currentFig.subplotList{1,1}.currentFig.gca.XTickLabel;
                            end
                        else
                            set(gca,'YLabel',[]);
                        end
                        set(gca,'ZLabel',[]);

                        set(self.currentFig.subplotList{irow,icol}.currentFig.gca,"FontSize",self.layout.axis.subplot.fontsize);
                        if icol ~= 1
                            set(self.currentFig.subplotList{irow,icol}.currentFig.gca,"YTickLabels",[]);
                        end
                        if irow ~= self.layout.axis.main.nrow
                            set(self.currentFig.subplotList{irow,icol}.currentFig.gca,"XTickLabels",[]);
                        end

                        self.currentFig.subplotList{irow,icol}.currentFig.gca.XLim = [ columns.min(icol) columns.max(icol) ];
                        if ~( strcmpi(currentPlotType,"histogram") || strcmpi(currentPlotType,"histfit") )
                            self.currentFig.subplotList{irow,icol}.currentFig.gca.YLim = [ columns.min(irow) columns.max(irow) ];
                        end

                        %set(self.currentFig.axisList(iaxis)(iaxis),'box',self.layout.axis.subplot.box);
                        hold off;

                    else

                        disp( "skipping subplot #" + string(iaxis+": ("+irow+","+icol+")") + msgSuffix )
                        self.currentFig.subplotList{irow,icol} = [];

                    end

                    grid(self.currentFig.subplotList{irow,icol}.currentFig.gca,self.layout.axis.subplot.grid);

                end
            end

            % set up colorbar

            if self.colorbar_kws.enabled && cEnabled
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

        function rotateAxesLabels(self,varargin)
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
                degreeX = varargin{1};
                degreeY = varargin{1};
            elseif nargin==3
                degreeX = varargin{1};
                degreeY = varargin{2};
            else
                error   ( newline ...
                        + "The rotateAxesLabels() method gets at most two input numeric values representing the X-axis and Y-axis label orientations from the horizontal axis, in degrees. " ...
                        + "If only one argument is provided, the input degree value will be used for both X and Y axis label orientations. " + newline ...
                        + "Alternatively, if no argument is passed, then all subplot axis labels will be rotated by the default 45 degrees." ...
                        + newline ...
                        );
            end

            for icol = 1:self.layout.axis.main.ncol
                for irow = 1:self.layout.axis.main.nrow
                    axisLabelX = get(self.currentFig.subplotList{irow,icol}.currentFig.gca,"XLabel"); set(axisLabelX,"rotation",degreeX,"VerticalAlignment","top","HorizontalAlignment","right");
                    axisLabelY = get(self.currentFig.subplotList{irow,icol}.currentFig.gca,"YLabel"); set(axisLabelY,"rotation",degreeY,"VerticalAlignment","middle","HorizontalAlignment","right");
                end
            end

        end % function rotateAxesLabels

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
            %               "hline_kws", {}     : the values specified in the cell array, are directly passed to `hline_kws`
            %                                   : component of Target_class object of each subplot.
            %               "vline_kws", {}     : the values specified in the cell array, are directly passed to `vline_kws`
            %                                   : component of Target_class object of each subplot.
            %               "scatter_kws", {}   : the values specified in the cell array, are directly passed to `scatter_kws`
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
                for kws = ["hline_kws","vline_kws","scatter_kws","values"]
                    [val, keyFound] = getKeyVal(kws,varargin{:});
                    if keyFound
                        props.(kws) = val;
                    end
                end
            end

            % compute the statistics

            if getVecLen(self.columns)
                [self.colnames, colindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.columns);
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

                    set(self.currentFig.gcf, "CurrentAxes", self.currentFig.subplotList{irow,icol}.currentFig.gca);

                    hold on;

                    self.currentFig.subplotList{irow,icol}.target.value = [ values(icol) values(irow) ];
                    for kws = ["hline_kws","vline_kws","scatter_kws"]
                        if isfield(props,kws)
                            for i = 1:2:length(props.(kws))
                                key = props.(kws){i};
                                val = props.(kws){i+1};
                                self.currentFig.subplotList{irow,icol}.target.(kws).(key) = val;
                            end
                        end
                    end
                    self.currentFig.subplotList{irow,icol}.target.plot();

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
                if ~iscell(labels) || labelsLen>self.layout.axis.main.nrow
                    error   ( "The input labels must be a cell array of length " + string(self.layout.axis.main.nrow) + " comprised of strings or char-vectors." ...
                            + "The input value is " + strjoin(string(labels),", ") ...
                            );
                end
            end
            for icol = 1:labelsLen % self.layout.axis.main.ncol
                for irow = 1:labelsLen % self.layout.axis.main.nrow
                    if strcmp(self.currentFig.subplotList{irow,icol}.currentFig.gca.Visible,"on")
                        if ~isnumeric(labels{icol})
                            if irow==labelsLen || (irow<labelsLen && strcmp(self.currentFig.subplotList{irow+1,icol}.currentFig.gca.Visible,"off"))
                                xlabel(self.currentFig.subplotList{irow,icol}.currentFig.gca, labels{icol}, varargin{:});
                            end
                        end
                        if ~isnumeric(labels{irow})
                            if icol==1 || (icol>1 && strcmp(self.currentFig.subplotList{irow,icol-1}.currentFig.gca.Visible,"off"))
                                ylabel(self.currentFig.subplotList{irow,icol}.currentFig.gca, labels{irow}, varargin{:});
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
                if limitsLen>self.layout.axis.main.nrow
                    error   ( "The input limits must be a string vector or a cell array of length " + string(self.layout.axis.main.nrow) + "." ...
                            + "The input value is " + strjoin(string(limits),", ") ...
                            );
                end
            end
            limitsElementLen = zeros(limitsLen,1);
            for icol = 1:limitsLen % self.layout.axis.main.ncol
                limitsElementLen(icol) = getVecLen(limits{icol});
                if limitsElementLen(icol)~=0 && limitsElementLen(icol)~=2
                    error   ( "The element " + string(icol) + " of the input cell array is invalid. " ...
                            + "The limits specified by the input cell must be either vectors of length two, or empty vectors. " ...
                            + "You have entered for the element " + string(icol) + ": " ...
                            + strjoin(string(limits{icol})) ...
                            );
                end
            end
            for icol = 1:limitsLen % self.layout.axis.main.ncol
                for irow = 1:limitsLen % self.layout.axis.main.nrow
                    if strcmp(self.currentFig.subplotList{irow,icol}.currentFig.gca.Visible,"on")
                        if limitsElementLen(icol)==2
                            xlim(self.currentFig.subplotList{irow,icol}.currentFig.gca, limits{icol}(:));
                        end
                        if limitsElementLen(irow)==2 && icol~=irow
                            ylim(self.currentFig.subplotList{irow,icol}.currentFig.gca, limits{irow}(:));
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
            methodName = varargin{1};
            methodNameIsShow = strcmp(methodName,"show");
            methodNameIsHide = strcmp(methodName,"hide");
            if methodNameIsShow
                switchValue = "on";
            elseif methodNameIsHide
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
                                + "Unrecognized input argument pass to the " + methodName + "() method: " + newline + newline ...
                                + join(part," ") + newline + newline ...
                                + "The " + methodName + "() method accepts only a list of comma-separated string values representing the ""part"" of the grid plot to " + methodName + ". " ...
                                + "Alternatively, if no argument is passed to " + methodName + "(), then all subplots will be affected." + newline ...
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

                if methodNameIsShow

                    colorbar_kws_cell = convertStruct2Cell(self.colorbar_kws,{"enabled","singleOptions"});

                    colormap(self.currentFig.gca,self.colormap);

                    % get column names

                    if getVecLen(self.ccolumn)
                        [ccolname, ccolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.ccolumn);
                    else
                        error("the component ""ccolumn"" of the GridPlot object appears to be empty.");
                    end

                    ccolumnValues = self.dfref.(ccolname);
                    colorbarLimit = [ min(ccolumnValues) max(ccolumnValues) ]; % Colorbar range
                    caxis(self.currentFig.gca,colorbarLimit);
                    set(self.currentFig.gca,"CLim",colorbarLimit);
                    self.currentFig.colorbar = colorbar(self.currentFig.gca,colorbar_kws_cell{:});

                    colorbarLabel = ccolname;
                    ylabel(self.currentFig.colorbar,colorbarLabel,"fontsize",self.layout.colorbar.fontsize);

                    set ( self.currentFig.colorbar ...
                        , "position" ...
                        , self.getColorbarAxisPosition() ...
                        , "FontSize", self.layout.colorbar.fontsize ...
                        ) ;

                elseif methodNameIsHide

                    %colorbar("off");
                    set(self.currentFig.colorbar,'visible','off');

                end

            end

            for icol = 1:self.layout.axis.main.ncol
                for irow = 1:self.layout.axis.main.nrow
                    if (icol<irow && lowerRequested) || (icol>irow && upperRequested) || (icol==irow && diagRequested)
                        set(self.currentFig.subplotList{irow,icol}.currentFig.gca,"visible",switchValue);
                        set(get(self.currentFig.subplotList{irow,icol}.currentFig.gca,"children"),"visible",switchValue);
                    end
                end
            end

            if lowerRequested || upperRequested || diagRequested
                self.hideShowAxesLabels();
            end

        end % function hideShow

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function hideShowAxesLabels(self)
            % show or hide axis labels and ticks depending on the presence of the neighbor subplots
            for icol = 1:self.layout.axis.main.ncol
                for irow = 1:self.layout.axis.main.nrow
                    if irow < self.layout.axis.main.nrow
                        if strcmp(self.currentFig.subplotList{irow+1,icol}.currentFig.gca.Visible,"on")
                            set(self.currentFig.subplotList{irow,icol}.currentFig.gca, "XTickLabel", []);
                            self.currentFig.subplotList{irow,icol}.currentFig.gca.XLabel.String = "";
                        else
                            set(self.currentFig.subplotList{irow,icol}.currentFig.gca, "XTickLabelMode", "auto");
                            self.currentFig.subplotList{irow,icol}.currentFig.gca.XLabel.String = self.colnames(icol);
                        end
                    end
                    if icol > 1
                        if strcmp(self.currentFig.subplotList{irow,icol-1}.currentFig.gca.Visible,"on")
                            set(self.currentFig.subplotList{irow,icol}.currentFig.gca, "YTickLabel", [])
                            self.currentFig.subplotList{irow,icol}.currentFig.gca.YLabel.String = "";
                        else
                            set(self.currentFig.subplotList{irow,icol}.currentFig.gca, "YTickLabelMode", "auto")
                            self.currentFig.subplotList{irow,icol}.currentFig.gca.YLabel.String = self.colnames(irow);
                        end
                    end
                    if ~(icol==1)
                        axisLabelY = get(self.currentFig.subplotList{irow,icol}.currentFig.gca,"YLabel"); set(axisLabelY,"rotation",45,"VerticalAlignment","middle","HorizontalAlignment","right");
                    end
                    if ~(irow==self.layout.axis.main.nrow)
                        axisLabelX = get(self.currentFig.subplotList{irow,icol}.currentFig.gca,"XLabel"); set(axisLabelX,"rotation",45,"VerticalAlignment","top","HorizontalAlignment","right");
                    end
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function position = getColorbarAxisPosition(self)
            position =  [ self.layout.colorbar.origin.x ...
                        , self.layout.colorbar.origin.y ...
                        , self.layout.colorbar.width ...
                        , self.layout.colorbar.height ...
                        ];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function position = getMainAxisPosition(self)
            position =  [ self.layout.axis.main.origin.x ...
                        , self.layout.axis.main.origin.y ...
                        , self.layout.axis.main.width ...
                        , self.layout.axis.main.height ...
                        ];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function position = getSubplotAxisPosition(self,irow,icol)
            position =  [ (icol-1)*(self.layout.axis.subplot.interspace + self.layout.axis.subplot.width)   + self.layout.axis.main.margin.left ...
                        , (irow-1)*(self.layout.axis.subplot.interspace + self.layout.axis.subplot.height)  + self.layout.axis.main.margin.bottom ...
                        , self.layout.axis.subplot.width ...
                        , self.layout.axis.subplot.height ...
                        ];
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
            %       For xample, change the left/bottom margin of the main axis of the figure to provide room for lengthy variable names.
            %       Then call the update() method of layout to reflect the changes.
            %
            %           update()
            %
            self.layout.axis.main.height = 1 - self.layout.axis.main.margin.top;
            self.layout.axis.main.width = 1 - self.layout.axis.main.margin.right;
            self.layout.axis.main.nrow = length(self.columns);
            self.layout.axis.main.ncol = self.layout.axis.main.nrow;
            self.layout.axis.subplot.height = (1-self.layout.axis.main.margin.top-self.layout.axis.main.margin.bottom)/self.layout.axis.main.ncol - self.layout.axis.subplot.interspace;
            self.layout.axis.subplot.width  = (1-self.layout.axis.main.margin.left-self.layout.axis.main.margin.right)/self.layout.axis.main.nrow - self.layout.axis.subplot.interspace;
            self.layout.colorbar.origin.x = 1 - self.layout.axis.main.margin.right;
            self.layout.colorbar.origin.y = self.layout.axis.main.margin.bottom;
            self.layout.colorbar.height = self.layout.axis.main.nrow * ( self.layout.axis.subplot.height + self.layout.axis.subplot.interspace ) - self.layout.axis.subplot.interspace;
            self.layout.figWidth    = self.layout.figHeight *   ( 1 ...
                                                                - self.layout.axis.main.margin.bottom ...
                                                                + self.layout.axis.main.margin.right ...
                                                                + self.layout.axis.main.margin.left ...
                                                                - self.layout.axis.main.margin.top ...
                                                                + self.layout.colorbar.width ...
                                                                );
            if isfield(self,"currentFig") && isfield(self.currentFig,"gca") && isgraphics(self.currentFig.gca)
                set(self.currentFig.gca,"position",self.getMainAxisPosition());
            end
            if any(strcmp(fieldnames(self.currentFig),"colorbar")) && isgraphics(self.currentFig.colorbar)
                set(self.currentFig.colorbar,"position",self.getColorbarAxisPosition());
            end
            if isfield(self.currentFig,"subplotList") && ~isempty(self.currentFig.subplotList)
                for irow = 1:self.layout.axis.main.nrow
                    for icol = 1:self.layout.axis.main.ncol
                        jrow = self.layout.axis.main.nrow-irow+1;
                        if isgraphics(self.currentFig.subplotList{irow,icol}.currentFig.gca)
                            set ( self.currentFig.subplotList{irow,icol}.currentFig.gca, "position", self.getSubplotAxisPosition(jrow,icol) );
                        end
                    end
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % classdef GridPlot