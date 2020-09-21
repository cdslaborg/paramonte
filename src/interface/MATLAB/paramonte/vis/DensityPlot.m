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
%   DensityPlot(dataFrame, plotType, resetExternal)
%
%   This is the DensityPlot class for generating instances of
%   histfit/histogram/histogram2/contourf/contour3/contour plots, built upon MATLAB's 
%   builtin functions: histogram, histogram2, histfit, contourf, contour3, contour.
%
%   Note that all objects instantiated from this class generate density plots of the input data. 
%   The histfit/histogram objects generate 1-dimensional density maps.
%   The histogram2/contourf/contour3/contour objects generate 2-dimensional density maps.
%   In case of contourf/contour3/contour objects, the density of the input data is first 
%   computed and smoothed by a Gaussian kernel density estimator and then passed 
%   to the MATLAB intrinsic plotting functions contourf/contour3/contour.
%
%   NOTE
%
%       This is a low-level ParaMonte class and is not meant
%       to be directly instantiated by the user.
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
%       ycolumns (available only in 2D and 3D objects)
%
%           Optional property that determines the columns of dataFrame 
%           to serve as the y-values. It can have multiple forms:
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
%               4.  ycolumns = 7:9      # every column in the data frame starting from column #7
%               5.  ycolumns = 7:2:20   # every other column in the data frame starting from column #7
%
%           WARNING
%
%               In all cases, ycolumns must have a length that is either 0, or 1, 
%               or equal to the length of xcolumns. If the length is 1, then ycolumns 
%               will be plotted against data corresponding to each element of xcolumns.
%               If it is an empty object having length 0, then a default value will be used.
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
%       histfit (available only in histfit objects)
%
%           A MATLAB struct() with the following components:
%
%               enabled
%
%                   A logical value. If `true`, the histfit will be added to the plot.
%
%               kws
%
%                   A MATLAB struct() whose components' values are passed to 
%                   MATLAB's histfit() function. If your desired attribute is 
%                   missing from the fieldnames of histfit.kws, simply add 
%                   a new field named as the attribute and assign the desired 
%                   value to it.
%
%           Example usage:
%
%               histfit.enabled = true; % add histfit()
%
%           WARNING
%
%               Keep in mind that MATLAB keyword arguments are case-INsensitive.
%               Hence, sure you do not add the keyword as multiple different fields.
%               For example, histfit.kws.color and histfit.kws.Color are the same,
%               and only one of the two will be processed.
%
%       histogram (available only in histogram objects)
%
%           A MATLAB struct() with the following components:
%
%               enabled
%
%                   A logical value. If `true`, the histogram will be added to the plot.
%
%               kws
%
%                   A MATLAB struct() whose components' values are passed to 
%                   MATLAB's histogram() function. If your desired attribute is 
%                   missing from the fieldnames of histogram.kws, simply add 
%                   a new field named as the attribute and assign the desired 
%                   value to it.
%
%           Example usage:
%
%               histogram.enabled = true; % add histogram()
%               histogram.kws.edgeColor = "none";
%
%           WARNING
%
%               Keep in mind that MATLAB keyword arguments are case-INsensitive.
%               Hence, sure you do not add the keyword as multiple different fields.
%               For example, histogram.kws.edgeColor and histogram.kws.EdgeColor are the same,
%               and only one of the two will be processed.
%
%       histogram2 (available only in histogram2 objects)
%
%           A MATLAB struct() with the following components:
%
%               enabled
%
%                   A logical value. If `true`, the histogram2 will be added to the plot.
%
%               kws
%
%                   A MATLAB struct() whose components' values are passed to 
%                   MATLAB's histogram2() function. If your desired attribute is 
%                   missing from the fieldnames of histogram2.kws, simply add 
%                   a new field named as the attribute and assign the desired 
%                   value to it.
%
%           Example usage:
%
%               histogram2.enabled = true; % add histogram2()
%               histogram2.kws.edgeColor = "none";
%
%           WARNING
%
%               Keep in mind that MATLAB keyword arguments are case-INsensitive.
%               Hence, sure you do not add the keyword as multiple different fields.
%               For example, histogram2.kws.edgeColor and histogram2.kws.EdgeColor are the same,
%               and only one of the two will be processed.
%
%       contourf (available only in contourf objects)
%
%           A MATLAB struct() with the following components:
%
%               enabled
%
%                   A logical value. If `true`, the contourf will be added to the plot.
%
%               kws
%
%                   A MATLAB struct() whose components' values are passed to 
%                   MATLAB's contourf() function. If your desired attribute is 
%                   missing from the fieldnames of contourf.kws, simply add 
%                   a new field named as the attribute and assign the desired 
%                   value to it.
%
%           Example usage:
%
%               contourf.enabled = true; % add contourf()
%               contourf.kws.color = "none";
%
%           WARNING
%
%               Keep in mind that MATLAB keyword arguments are case-INsensitive.
%               Hence, sure you do not add the keyword as multiple different fields.
%               For example, contourf.kws.color and contourf.kws.Color are the same,
%               and only one of the two will be processed.
%
%       contour3 (available only in contour3 objects)
%
%           A MATLAB struct() with the following components:
%
%               enabled
%
%                   A logical value. If `true`, the contour3 will be added to the plot.
%
%               kws
%
%                   A MATLAB struct() whose components' values are passed to 
%                   MATLAB's contour3() function. If your desired attribute is 
%                   missing from the fieldnames of contour3.kws, simply add 
%                   a new field named as the attribute and assign the desired 
%                   value to it.
%
%           Example usage:
%
%               contour3.enabled = true; % add contour3()
%               contour3.kws.lineWidth = "none";
%
%           WARNING
%
%               Keep in mind that MATLAB keyword arguments are case-INsensitive.
%               Hence, sure you do not add the keyword as multiple different fields.
%               For example, contour3.kws.lineWidth and contour3.kws.WineWidth are the same,
%               and only one of the two will be processed.
%
%       contour (available only in contour objects)
%
%           A MATLAB struct() with the following components:
%
%               enabled
%
%                   A logical value. If `true`, the contour will be added to the plot.
%
%               kws
%
%                   A MATLAB struct() whose components' values are passed to 
%                   MATLAB's contour() function. If your desired attribute is 
%                   missing from the fieldnames of contour.kws, simply add 
%                   a new field named as the attribute and assign the desired 
%                   value to it.
%
%           Example usage:
%
%               contour.enabled = true; % add contour()
%               contour.kws.lineWidth = "none";
%
%           WARNING
%
%               Keep in mind that MATLAB keyword arguments are case-INsensitive.
%               Hence, sure you do not add the keyword as multiple different fields.
%               For example, contour.kws.lineWidth and contour.kws.WineWidth are the same,
%               and only one of the two will be processed.
%
%       gridSize (available in ``contour``, ``contourf``, ``contour3`` objects)
%
%           An integer indicating the grid resolution for discretization of the 
%           data during the kernel density estimation. It must be a power of 
%           two, otherwise it will be changed to the next power of two at the
%           time of using it. The default value is 512.
%
%           Example usage:
%
%                   gridSize = 512
%
%       noiseDensity (available in ``contour``, ``contourf``, ``contour3``)
%
%           A float indicating the threshold below which the kernel density 
%           estimate is considered to be noise and is rounded to zero.
%           The higher this value is, the less noise will be visible.
%
%           Example usage:
%
%                   noiseDensity = 1.e-5
%
%       target (available in 2D plots)
%
%           An object of class Target_class for adding target values to the plots.
%           For more information, see the documentation for the Target_class().
%
%   Superclass Attributes
%   ---------------------
%
%       See the documentation for the BasePlot class.
%
%   Returns
%   -------
%
%       An object of DensityPlot class.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
classdef DensityPlot < BasePlot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Access = public)
        xcolumns
    end

    properties (Access = protected, Hidden)
        noiseLevelDefault = 0.001;
        gridSizeDefault = 2^9;
        isHistogram = false;
        isHistogram2 = false;
        isContourf = false;
        isContour3 = false;
        isContour = false;
        isHistfit = false;
        xcolnames
        xcolindex
        ycolnames
        ycolindex
        ccolnames
        ccolindex
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function resetInternal(self)

            resetInternal@BasePlot(self);
            self.xcolumns = {};

            if self.type.isHistogram

                self.histogram.kws = struct();
                self.histogram.enabled = true;
                self.histogram.kws.edgeColor = {};

            elseif self.type.isHistogram2 || self.type.isContourf || self.type.isContour3 || self.type.isContour

                self.ycolumns = {};

                self.colormap = struct();
                self.colormap.enabled = true;
                self.colormap.values = [];

                if self.type.isHistogram2
                    self.histogram2.kws = struct();
                    self.histogram2.enabled = true;
                    self.histogram2.kws.edgeColor = {};
                    self.histogram2.kws.faceColor = {};
                    self.histogram2.kws.displayStyle = {};
                end

                %if self.type.isContourf || self.type.isContour3 || self.type.isContour
                %    self.gridSize = self.gridSizeDefault;
                %    self.gridSize = self.gridSizeDefault;
                %end

                if self.type.isContourf || self.type.isContour3 || self.type.isContour
                    self.(self.type.name).kws = struct();
                    self.(self.type.name).enabled = true;
                end

                self.colorbar.kws = struct();
                self.colorbar.enabled = true;
                self.colorbar.kws.fontSize = [];

                self.legend.enabled = false;

            elseif self.type.isHistfit

                %self.histfit.kws = struct();
                self.histfit.enabled = true;

            end

            self.isdryrun = true;
            self.make();
            self.isdryrun = false;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = DensityPlot(plotType, dataFrame, resetExternal)
            if nargin<3; resetExternal = []; end
            self = self@BasePlot(plotType, dataFrame, resetExternal);
            if nargin<3; self.resetExternal = @self.resetInternal; end
            if self.type.isHistogram
                prop="histogram"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            elseif self.type.isHistfit
                prop="histfit"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            elseif self.type.isHistogram2 || self.type.isContourf || self.type.isContour3 || self.type.isContour
                prop=self.type.name; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                prop="colorbar"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                prop="ycolumns"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                prop="colormap"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                if ~self.type.isHistogram2
                    prop="gridSize"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end; self.(prop) = [];
                    prop="noiseLevel"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end; self.(prop) = [];
                end
            else
                error("The requested plot type " + string(self.type.name) + " is not reqcognized.");
            end
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
            %       make("xcolumns",[8])
            %       make("xcolumns",7:10)
            %       make( "gcf_kws", {"enabled",true,"color","none"} )
            %
            %

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% parse arguments
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.parseArgs(varargin{:});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% set what to plot
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if self.type.isHistogram % && self.histogram.enabled
                fname = "histogram";
                key = "edgeColor"; val = "none"; if ~isfield(self.(fname).kws,key) || isempty(self.(fname).kws.(key)); self.(fname).kws.(key) = val; end
            end

            if self.type.isHistfit % && self.histogram.enabled
                fname = "histfit";
                key = "dist"; val = "Normal"; if ~isfield(self.(fname),key) || isempty(self.(fname).(key)); self.(fname).(key) = val; end
                key = "nbins"; val = []; if ~isfield(self.(fname),key) || isempty(self.(fname).(key)); self.(fname).(key) = val; end
            end

            if self.type.isHistogram2 % && self.histogram2.enabled
                fname = "histogram2";
                keyList = ["showemptybins", "edgeColor", "displayStyle"]; % "numbins"
                valueList = {"off", "none", "bar3"}; % [100 100]
                for j = 1:length(keyList)
                    key = keyList(j);
                    if ~isfield(self.(fname).kws,key) || isempty(self.(fname).kws.(key))
                        self.(fname).kws.(key) = valueList{j};
                    end
                end
                if self.colormap.enabled
                    key = "faceColor"; val = "flat"; if ~isfield(self.(fname).kws,key) || isempty(self.(fname).kws.(key)); self.(fname).kws.(key) = val; end
                else
                    key = "faceColor"; val = "auto"; self.(fname).kws.(key) = val;
                end
            end

            if self.type.isContour || self.type.isContourf || self.type.isContour3
                fname = self.type.name;
                if self.type.isContour3
                    self.(self.type.name).levels = 50;
                    key = "color"; val = []; if ~isfield(self.(fname).kws,key) || isempty(self.(fname).kws.(key)); self.(fname).kws.(key) = val; end
                else
                    self.(self.type.name).levels = 8;
                    if self.type.isContourf
                        key = "color"; val = "none"; if ~isfield(self.(fname).kws,key) || isempty(self.(fname).kws.(key)); self.(fname).kws.(key) = val; end
                        if ~self.colormap.enabled && self.contourf.enabled
                            self.colormap.values = flipud(gray);
                        end
                    else
                        key = "color"; val = []; if ~isfield(self.(fname).kws,key) || isempty(self.(fname).kws.(key)); self.(fname).kws.(key) = val; end
                    end
                end
                key = "gridSize"; val = self.gridSizeDefault; if isempty(self.(key)); self.(key) = val; end
                key = "noiseLevel"; val = self.noiseLevelDefault; if isempty(self.(key)); self.(key) = val; end
                key = "lineStyle"; val = "-"; if ~isfield(self.(fname).kws,key) || isempty(self.(fname).kws.(key)); self.(fname).kws.(key) = val; end
                key = "lineWidth"; val = 0.5; if ~isfield(self.(fname).kws,key) || isempty(self.(fname).kws.(key)); self.(fname).kws.(key) = val; end
                key = "showText"; val = "off"; if ~isfield(self.(fname).kws,key) || isempty(self.(fname).kws.(key)); self.(fname).kws.(key) = val; end
                key = "labelSpacing"; val = 144; if ~isfield(self.(fname).kws,key) || isempty(self.(fname).kws.(key)); self.(fname).kws.(key) = val; end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if self.isdryrun; return; end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% check columns presence
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if getVecLen(self.xcolumns)
                [self.xcolnames, self.xcolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.xcolumns);
            else
                self.self.xcolindex = [];
                self.xcolnames = "Count";
                xdata = self.rowsindex;
            end
            xcolindexlen = length(self.xcolindex);
            maxLenColumns = xcolindexlen;
            if xcolindexlen==1
                xdata = self.dfref{self.rowsindex,self.xcolindex};
            end

            if self.type.isHistogram2 || self.type.isContourf || self.type.isContour3 || self.type.isContour
                [self.ycolnames, self.ycolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.ycolumns);
                ycolindexlen = length(self.ycolindex);
                maxLenColumns = max(maxLenColumns,ycolindexlen);
                if ycolindexlen~=maxLenColumns && ycolindexlen>1
                    error("length of ycolumns must be either 1 or equal to the maximum of the lengths of xcolumns.");
                end
                if ycolindexlen==1
                    ydata = self.dfref{self.rowsindex,self.ycolindex};
                end
            end

            if xcolindexlen~=maxLenColumns && xcolindexlen>1
                error("length of xcolumns must be either 1 or equal to the maximum of the lengths of ycolumns.");
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% get keyword arguments
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if self.type.isHistogram && isfield(self.histogram,"kws")
                histogram_kws_cell = convertStruct2Cell(self.histogram.kws,{"enabled","singleOptions"});
            else
                histogram_kws_cell = {};
            end
            if self.type.isHistogram2 && isfield(self.histogram2,"kws")
                histogram2_kws_cell = convertStruct2Cell(self.histogram2.kws,{"enabled","singleOptions"});
            else
                histogram2_kws_cell = {};
            end
            if self.type.isHistfit && isfield(self.histfit,"kws")
                histfit_kws_cell = convertStruct2Cell(self.histfit.kws,{"enabled","singleOptions"});
            else
                histfit_kws_cell = {};
            end
            if (self.type.isContour3 || self.type.isContourf || self.type.isContour) && isfield(self.(self.type.name),"kws")
                kws_cell.(self.type.name) = convertStruct2Cell(self.(self.type.name).kws,{"enabled","levels","singleOptions"});
            else
                kws_cell.(self.type.name) = {};
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% add plots
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %box on; grid on;

            lglabels = [];

            if ((self.type.isContourf && self.contourf.enabled) ...
            ||  (self.type.isContour3 && self.contour3.enabled) ...
            ||  (self.type.isContour && self.contour.enabled) ) ...
            &&  (xcolindexlen==1 && ycolindexlen==1) ...
                self.currentFig.kde2d = struct();
                [ self.currentFig.kde2d.bandwidth ...
                , self.currentFig.kde2d.density ...
                , self.currentFig.kde2d.crdx ...
                , self.currentFig.kde2d.crdy ...
                ] = kde2d([xdata(:),ydata(:)],self.gridSize);
                self.currentFig.kde2d.density(self.currentFig.kde2d.density<self.noiseLevel) = NaN;
            end

            if (self.type.isHistogram && self.histogram.enabled) ...
            || (self.type.isHistogram2 && self.histogram2.enabled) ...
            || (self.type.isHistfit && self.histfit.enabled) ...
            || (self.type.isContourf && self.contourf.enabled) ...
            || (self.type.isContour3 && self.contour3.enabled) ...
            || (self.type.isContour && self.contour.enabled)

                isMultiColorContourPlot = false;
                if (self.type.isContour || self.type.isContourf || self.type.isContour3) && ~self.colormap.enabled
                    if getVecLen(self.(self.type.name).kws.color)
                            if all(size(self.(self.type.name).kws.color)==[maxLenColumns,3])
                                isMultiColorContourPlot = true;
                                contourMultiColorList = self.(self.type.name).kws.color;
                            elseif ~( all(size(self.(self.type.name).kws.color)==[1,3]) || self.type.isContourf )
                                error   ( "The value specified for the '" + self.type.name + ".kws.color' property of the plot object must be either " ...
                                        + "an RGB triplet vector of size [1,3], or a matrix of size [" + string(maxLenColumns) + ",3] for the " ...
                                        + "current set of variables that are selected to plot. It can also be an empty object, in which case, " ...
                                        + "the colors of the objects in the plot will be chosen automatically." ...
                                        );
                            end
                    else
                        contourMultiColorList = lines(maxLenColumns);
                        isMultiColorContourPlot = true;
                    end
                end

                counter.histfit = 0;
                counter.histogram = 0;
                counter.histogram2 = 0;
                counter.contour3 = 0;
                counter.contourf = 0;
                counter.contour = 0;

                for icol = 1:maxLenColumns

                    if isMultiColorContourPlot % true only for contour plots
                        self.(self.type.name).kws.color = contourMultiColorList(icol,:);
                        if isfield(self.(self.type.name),"kws")
                            kws_cell.(self.type.name) = convertStruct2Cell(self.(self.type.name).kws,{"enabled","levels","singleOptions"});
                        end
                    end

                    kde2dComputationNeeded = false;
                    if xcolindexlen>1
                        xdata = self.dfref{self.rowsindex,self.xcolindex(icol)};
                        if self.type.isContourf && self.contourf.enabled; kde2dComputationNeeded = true; end
                        if self.type.isContour3 && self.contour3.enabled; kde2dComputationNeeded = true; end
                        if self.type.isContour && self.contour.enabled; kde2dComputationNeeded = true; end
                    end
                    if (self.type.isHistogram2 && self.histogram2.enabled) ...
                    || (self.type.isContourf && self.contourf.enabled) ...
                    || (self.type.isContour3 && self.contour3.enabled) ...
                    || (self.type.isContour && self.contour.enabled)
                        if ycolindexlen>1
                            ydata = self.dfref{self.rowsindex,self.ycolindex(icol)};
                            if self.type.isContourf && self.contourf.enabled; kde2dComputationNeeded = true; end
                            if self.type.isContour3 && self.contour3.enabled; kde2dComputationNeeded = true; end
                            if self.type.isContour && self.contour.enabled; kde2dComputationNeeded = true; end
                        end
                    end
                    if kde2dComputationNeeded
                        [ self.currentFig.kde2d.bandwidth ...
                        , self.currentFig.kde2d.density ...
                        , self.currentFig.kde2d.crdx ...
                        , self.currentFig.kde2d.crdy ...
                        ] = kde2d( [xdata(:),ydata(:)], self.gridSize );
                        self.currentFig.kde2d.density(self.currentFig.kde2d.density<self.noiseLevel) = NaN;
                    end

                    if self.legend.enabled
                        if (self.type.isHistogram2 && self.histogram2.enabled) ...
                        ||  (self.type.isContourf && self.contourf.enabled) ...
                        ||  (self.type.isContour3 && self.contour3.enabled) ...
                        ||  (self.type.isContour && self.contour.enabled)
                            if xcolindexlen<2 && ycolindexlen>=1
                                lglabels = [ lglabels , self.ycolnames(icol) ];
                            elseif xcolindexlen>1 && ycolindexlen<2
                                lglabels = [ lglabels , self.xcolnames(icol) ];
                            else
                                lglabels = [ lglabels , self.xcolnames(icol)+"-"+self.ycolnames(icol) ];
                            end
                        else
                            lglabels = [ lglabels , self.xcolnames(icol) ];
                        end
                    end

                    % add plot

                    if self.type.isHistogram && self.histogram.enabled

                        counter.histogram = counter.histogram + 1;

                        self.currentFig.histogram{counter.histogram} = histogram( xdata ...
                                                                                , histogram_kws_cell{:} ...
                                                                                );
                        hold on;

                    elseif self.type.isHistogram2 && self.histogram2.enabled

                        counter.histogram2 = counter.histogram2 + 1;

                        self.currentFig.histogram2{counter.histogram2} = histogram2 ( xdata ...
                                                                                    , ydata ...
                                                                                    , histogram2_kws_cell{:} ...
                                                                                    );
                        hold on;

                    elseif self.type.isContourf && self.contourf.enabled

                        counter.contourf = counter.contourf + 1;

                        self.currentFig.contourf{counter.contourf} = struct();
                        self.currentFig.contourf{counter.contourf}.matrix = [];
                        self.currentFig.contourf{counter.contourf}.object = [];
                        [ self.currentFig.contourf{counter.contourf}.matrix ...
                        , self.currentFig.contourf{counter.contourf}.object ] = contourf( self.currentFig.kde2d.crdx ...
                                                                                        , self.currentFig.kde2d.crdy ...
                                                                                        , self.currentFig.kde2d.density ...
                                                                                        , self.contourf.levels ...
                                                                                        , kws_cell.(self.type.name){:} ...
                                                                                        );
                        hold on;

                    elseif self.type.isContour3 && self.contour3.enabled

                        counter.contour3 = counter.contour3 + 1;

                        self.currentFig.contour3{counter.contour3} = struct();
                        self.currentFig.contour3{counter.contour3}.matrix = [];
                        self.currentFig.contour3{counter.contour3}.object = [];
                        [ self.currentFig.contour3{counter.contour3}.matrix ...
                        , self.currentFig.contour3{counter.contour3}.object ] = contour3( self.currentFig.kde2d.crdx ...
                                                                                        , self.currentFig.kde2d.crdy ...
                                                                                        , self.currentFig.kde2d.density ...
                                                                                        , self.contour3.levels ...
                                                                                        , kws_cell.(self.type.name){:} ...
                                                                                        );
                        hold on;

                    elseif self.type.isContour && self.contour.enabled

                        counter.contour = counter.contour + 1;

                        self.currentFig.contour{counter.contour} = struct();
                        self.currentFig.contour{counter.contour}.matrix = [];
                        self.currentFig.contour{counter.contour}.object = [];
                        [ self.currentFig.contour{counter.contour}.matrix ...
                        , self.currentFig.contour{counter.contour}.object ] = contour   ( self.currentFig.kde2d.crdx ...
                                                                                        , self.currentFig.kde2d.crdy ...
                                                                                        , self.currentFig.kde2d.density ...
                                                                                        , self.contour.levels ...
                                                                                        , kws_cell.(self.type.name){:} ...
                                                                                        );
                        hold on;

                    elseif self.type.isHistfit && self.histfit.enabled

                        counter.histfit = counter.histfit + 1;

                        self.currentFig.histfit{counter.histfit} = histfit  ( xdata ...
                                                                            , self.histfit.nbins ...
                                                                            , self.histfit.dist ...
                                                                            , histfit_kws_cell{:} ...
                                                                            );
                        self.currentFig.histfit{counter.histfit}(1).EdgeAlpha = 0; % make histogram lines invisible
                        self.currentFig.histfit{counter.histfit}(1).FaceAlpha = 0.6; % make histogram face transparent 60%
                        self.currentFig.histfit{counter.histfit}(2).Color = "black"; % make the fitted line's color black
                        hold on;
                    end

                end % loop plot

            end

            %box on; grid on;

            self.currentFig.axes = gca;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% add colormap
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if self.type.isHistogram2 || self.type.isContourf || self.type.isContour3 || self.type.isContour
                if self.colormap.enabled && ~getVecLen(self.colormap.values)
                    if self.type.isContour3 || self.type.isContour
                        self.colormap.values = "parula";
                    else
                        self.colormap.values = "parula"; % flipud(cold());
                    end
                end
                colormap(self.currentFig.axes,self.colormap.values);
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% add axes labels
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if xcolindexlen>1
                self.currentFig.xlabel = xlabel("Variable Values", "Interpreter", "none");
            else
                self.currentFig.xlabel = xlabel(self.xcolnames(1), "Interpreter", "none");
            end

            if self.type.isHistogram2 || self.type.isContourf || self.type.isContour3 || self.type.isContour
                if ycolindexlen>1
                    self.currentFig.ylabel = ylabel("Variable Values", "Interpreter", "none");
                else
                    self.currentFig.ylabel = ylabel(self.ycolnames(1), "Interpreter", "none");
                end
                self.currentFig.zlabel = zlabel("Density of Points", "Interpreter", "none");
            else
                self.currentFig.ylabel = ylabel("Density of Points", "Interpreter", "none");
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% add colorbar
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            colorbarEnabled = (self.type.isHistogram2 || self.type.isContourf || self.type.isContour3 || self.type.isContour) ...
                            && self.colorbar.enabled && self.colormap.enabled;
            if colorbarEnabled
                if isempty(self.colorbar.kws.fontSize) || ~isa(self.colorbar.kws.fontSize,"numeric")
                    self.colorbar.kws.fontSize = self.currentFig.ylabel.FontSize;
                end
                colorbar_kws_cell = convertStruct2Cell(self.colorbar.kws,{"enabled","singleOptions"});
                self.currentFig.colorbar = colorbar(colorbar_kws_cell{:});
                ylabel(self.currentFig.colorbar,"Density of Points","fontsize",self.colorbar.kws.fontSize, "Interpreter", "none");
            else
                colorbar('off');
                self.currentFig.colorbar = [];
            end

            if self.legend.enabled && (~isfield(self.legend,"labels") || isempty(self.legend.labels)); self.legend.labels = lglabels; end
            self.doBasePlotStuff();

            hold off;

        end % function plot

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % classdef DensityPlot
