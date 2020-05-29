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
%   HistPlot(dataFrame, plotType)
%
%   This is the HistPlot class for generating instances of
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
%       ycolumns (available only in histogram2/contourf/contour3/contour objects)
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
%               4.  ycolumns = 7:9      # every column in the data frame starting from column #7
%               5.  ycolumns = 7:2:20   # every other column in the data frame starting from column #7
%
%           WARNING: In all cases, ycolumns must have a length that is either 0, or 1, or equal
%           WARNING: to the length of xcolumns. If the length is 1, then ycolumns will be
%           WARNING: plotted against data corresponding to each element of ycolumns.
%           WARNING: If it is an empty object having length 0, then the default value will be used.
%
%           The default value is the names of all columns of the input dataFrame.
%
%       colormap (available only in histogram2/contourf/contour3/contour objects) 
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
%           If colormap is not provided or is empty, the default will be "winter".
%
%       colorbar_kws (available only in histogram2/contourf/contour3/contour objects) 
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
%       histfit_kws (available only in histfit objects)
%
%           A MATLAB struct() whose fields (with the exception of few, e.g., enabled, singleOptions, ...)
%           are directly passed to the `histfit()` function of MATLAB.
%           This property exists only if the object is instantiated as a histogram object.
%
%           Example usage:
%
%               histfit_kws.enabled = true; % add histfit()
%
%           If a desired property is missing among the struct fields, simply add the field
%           and its value to histfit_kws.
%
%           WARNING: keep in mind that MATLAB keyword arguments are case-INsensitive.
%           WARNING: therefore make sure you do not add the keyword as multiple different fields.
%           WARNING: For example, histfit_kws.Enabled and histfit_kws.enabled are the same,
%           WARNING: and only one of the two will be processed.
%
%       histogram_kws (available only in histogram objects)
%
%           A MATLAB struct() whose fields (with the exception of few, e.g., enabled, singleOptions, ...)
%           are directly passed to the `histogram()` function of MATLAB.
%           This property exists only if the object is instantiated as a histogram object.
%
%           Example usage:
%
%               histogram_kws.enabled = true; % add histogram()
%               histogram_kws.edgecolor = "none";
%
%           If a desired property is missing among the struct fields, simply add the field
%           and its value to histogram_kws.
%
%           WARNING: keep in mind that MATLAB keyword arguments are case-INsensitive.
%           WARNING: therefore make sure you do not add the keyword as multiple different fields.
%           WARNING: histogram_kws.edgecolor and histogram_kws.Edgecolor are the same,
%           WARNING: and only one of the two will be processed.
%
%       histogram2_kws (available only in histogram2 objects)
%
%           A MATLAB struct() whose fields (with the exception of few, e.g., enabled, singleOptions, ...)
%           are directly passed to the `histogram2()` function of MATLAB.
%           This property exists only if the object is instantiated as a histogram object.
%
%           Example usage:
%
%               histogram2_kws.enabled = true; % add histogram2()
%               histogram2_kws.edgecolor = "none";
%
%           If a desired property is missing among the struct fields, simply add the field
%           and its value to histogram2_kws.
%
%           WARNING: keep in mind that MATLAB keyword arguments are case-INsensitive.
%           WARNING: therefore make sure you do not add the keyword as multiple different fields.
%           WARNING: histogram2_kws.edgecolor and histogram2_kws.Edgecolor are the same,
%           WARNING: and only one of the two will be processed.
%
%       contourf_kws (available only in contourf objects)
%
%           A MATLAB struct() whose fields (with the exception of few, e.g., enabled, singleOptions, ...)
%           are directly passed to the `contourf()` function of MATLAB.
%           This property exists only if the object is instantiated as a histogram object.
%
%           Example usage:
%
%               contourf_kws.enabled = true; % add contourf()
%               contourf_kws.linecolor = "none";
%
%           If a desired property is missing among the struct fields, simply add the field
%           and its value to contourf_kws.
%
%           WARNING: keep in mind that MATLAB keyword arguments are case-INsensitive.
%           WARNING: therefore make sure you do not add the keyword as multiple different fields.
%           WARNING: contourf_kws.linecolor and contourf_kws.LineColor are the same,
%           WARNING: and only one of the two will be processed.
%
%       contour3_kws (available only in contour3 objects)
%
%           A MATLAB struct() whose fields (with the exception of few, e.g., enabled, singleOptions, ...)
%           are directly passed to the `contour3()` function of MATLAB.
%           This property exists only if the object is instantiated as a histogram object.
%
%           Example usage:
%
%               contour3_kws.enabled = true; % add contour3()
%               contour3_kws.linewidth = 1;
%
%           If a desired property is missing among the struct fields, simply add the field
%           and its value to contour3_kws.
%
%           WARNING: keep in mind that MATLAB keyword arguments are case-INsensitive.
%           WARNING: therefore make sure you do not add the keyword as multiple different fields.
%           WARNING: contour3_kws.linewidth and contour3_kws.LineWidth are the same,
%           WARNING: and only one of the two will be processed.
%
%       contour_kws (available only in contour objects)
%
%           A MATLAB struct() whose fields (with the exception of few, e.g., enabled, singleOptions, ...)
%           are directly passed to the `contour()` function of MATLAB.
%           This property exists only if the object is instantiated as a histogram object.
%
%           Example usage:
%
%               contour_kws.enabled = true; % add contour()
%               contour_kws.linewidth = 1;
%
%           If a desired property is missing among the struct fields, simply add the field
%           and its value to contour_kws.
%
%           WARNING: keep in mind that MATLAB keyword arguments are case-INsensitive.
%           WARNING: therefore make sure you do not add the keyword as multiple different fields.
%           WARNING: contour_kws.linewidth and contour_kws.LineWidth are the same,
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
%       An object of HistPlot class
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
classdef HistPlot < BasePlot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Access = public)
        xcolumns
        target
    end

    properties (Access = protected, Hidden)
        plotType
        isHistogram = false;
        isHistogram2 = false;
        isContourf = false;
        isContour3 = false;
        isContour = false;
        isHistfit = false;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function reset(self)

            reset@BasePlot(self);
            self.xcolumns = {};

            if self.isHistogram

                self.histogram_kws = struct();
                self.histogram_kws.enabled = true;
                self.histogram_kws.edgecolor = {};
                self.histogram_kws.singleOptions = {};

            elseif self.isHistogram2 || self.isContourf || self.isContour3 || self.isContour

                self.ycolumns = {};
                self.colormap = {};

                if self.isHistogram2
                    self.histogram2_kws = struct();
                    self.histogram2_kws.enabled = true;
                    self.histogram2_kws.edgecolor = {};
                    self.histogram2_kws.facecolor = {};
                    self.histogram2_kws.displaystyle = {};
                    self.histogram2_kws.singleOptions = {};
                end

                if self.isContourf
                    self.contourf_kws = struct();
                    self.contourf_kws.enabled = true;
                    self.contourf_kws.levels = 8;
                    self.contourf_kws.linecolor = "none";
                    self.contourf_kws.singleOptions = {};
                end

                if self.isContour3
                    self.contour3_kws = struct();
                    self.contour3_kws.enabled = true;
                    self.contour3_kws.levels = 50;
                    self.contour3_kws.singleOptions = {};
                end

                if self.isContour
                    self.contour_kws = struct();
                    self.contour_kws.enabled = true;
                    self.contour_kws.levels = 8;
                    self.contour_kws.singleOptions = {};
                end

                self.colorbar_kws = struct();
                self.colorbar_kws.enabled = true;
                self.colorbar_kws.fontsize = [];
                self.colorbar_kws.singleOptions = {};

                self.legend_kws.enabled = false;

            elseif self.isHistfit

                self.histfit_kws = struct();
                self.histfit_kws.enabled = true;
                self.histfit_kws.singleOptions = {};

            end

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

        function self = HistPlot(varargin) % expected input arguments: dataFrame, plotType
            self = self@BasePlot(varargin{1});
            self.plotType = lower(string(varargin{2}));
            if strcmpi(self.plotType,"histogram")
                self.isHistogram = true;
                prop="histogram_kws"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            elseif strcmpi(self.plotType,"histfit")
                self.isHistfit = true;
                prop="histfit_kws"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            else
                if strcmpi(self.plotType,"histogram2"); self.isHistogram2 = true; dummy_kws="histogram2_kws"; end
                if strcmpi(self.plotType,"contourf"); self.isContourf = true; dummy_kws="contourf_kws"; end
                if strcmpi(self.plotType,"contour3"); self.isContour3 = true; dummy_kws="contour3_kws"; end
                if strcmpi(self.plotType,"contour"); self.isContour = true; dummy_kws="contour_kws"; end
                if self.isHistogram2 || self.isContourf || self.isContour3 || self.isContour
                    prop=dummy_kws; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                    prop="colorbar_kws"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                    prop="ycolumns"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                    prop="colormap"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                else
                    error("The requested plot type " + string(self.plotType) + " is not reqcognized.");
                end
            end
            self.reset();
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
            %       plot("xcolumns",[8])
            %       plot("xcolumns",7:10)
            %       plot( "gcf_kws", {"enabled",true,"color","none"} )
            %
            %

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% parse arguments
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.parseArgs(varargin{:});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% set what to plot
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if self.isHistogram && self.histogram_kws.enabled
                fname = "histogram_kws";
                key = "edgecolor"; val = "none"; if isfield(self.(fname),key) && isempty(self.(fname).(key)); self.(fname).(key) = val; end
            end

            if self.isHistogram2 && self.histogram2_kws.enabled
                fname = "histogram2_kws";
                key = "edgecolor"; val = "none"; if isfield(self.(fname),key) && isempty(self.(fname).(key)); self.(fname).(key) = val; end
                key = "facecolor"; val = "flat"; if isfield(self.(fname),key) && isempty(self.(fname).(key)); self.(fname).(key) = val; end
                key = "displaystyle"; val = "bar3"; if isfield(self.(fname),key) && isempty(self.(fname).(key)); self.(fname).(key) = val; end
            end

            if self.isContourf && self.contourf_kws.enabled
                fname = "contourf_kws";
            end

            if self.isContour3 && self.contour3_kws.enabled
                fname = "contour3_kws";
            end
            if self.isContour && self.contour_kws.enabled
                fname = "contour_kws";
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

            % check columns presence

            if getVecLen(self.xcolumns)
                [xcolnames, xcolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.xcolumns);
            else
                xcolindex = [];
                xcolnames = "Count";
                xdata = rowindex;
            end
            xcolindexlen = length(xcolindex);
            maxLenColumns = xcolindexlen;
            if xcolindexlen==1
                xdata = self.dfref{rowindex,xcolindex};
            end

            if self.isHistogram2 || self.isContourf || self.isContour3 || self.isContour
                [ycolnames, ycolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.ycolumns);
                ycolindexlen = length(ycolindex);
                maxLenColumns = max(maxLenColumns,ycolindexlen);
                if ycolindexlen~=maxLenColumns && ycolindexlen>1; error("length of ycolumns must be either 1 or equal to the maximum of the lengths of xcolumns."); end
                if ycolindexlen==1
                    ydata = self.dfref{rowindex,ycolindex};
                end
            end

            if xcolindexlen~=maxLenColumns && xcolindexlen>1; error("length of xcolumns must be either 1 or equal to the maximum of the lengths of ycolumns."); end

            %%%%%%%%%%%%%%%%%%%%%%%
            % get keyword arguments
            %%%%%%%%%%%%%%%%%%%%%%%

            if self.isHistogram
                histogram_kws_cell = convertStruct2Cell(self.histogram_kws,{"enabled","singleOptions"});
            end
            if self.isHistogram2
                histogram2_kws_cell = convertStruct2Cell(self.histogram2_kws,{"enabled","singleOptions"});
            end
            if self.isContourf
                contourf_kws_cell = convertStruct2Cell(self.contourf_kws,{"enabled","levels","singleOptions"});
            end
            if self.isContour3
                contour3_kws_cell = convertStruct2Cell(self.contour3_kws,{"enabled","levels","singleOptions"});
            end
            if self.isContour
                contour_kws_cell = convertStruct2Cell(self.contour_kws,{"enabled","levels","singleOptions"});
            end
            if self.isHistfit
                histfit_kws_cell = convertStruct2Cell(self.histfit_kws,{"enabled","singleOptions"});
            end

            %%%%%%%%%%%%%%%
            % add line plot
            %%%%%%%%%%%%%%%

            box on; grid on;

            lglabels = [];

            if (self.isContourf && self.contourf_kws.enabled) || (self.isContour3 && self.contour3_kws.enabled) || (self.isContour && self.contour_kws.enabled)
                gridSize = 2^9;
                [bandwidth,density,crdx,crdy] = kde2d([xdata(:),ydata(:)],gridSize);
            end

            if (self.isHistogram && self.histogram_kws.enabled) ...
            || (self.isHistogram2 && self.histogram2_kws.enabled) ...
            || (self.isHistfit && self.histfit_kws.enabled) ...
            || (self.isContourf && self.contourf_kws.enabled) ...
            || (self.isContour3 && self.contour3_kws.enabled) ...
            || (self.isContour && self.contour_kws.enabled)

                for i = 1:maxLenColumns

                    kde2dComputationNeeded = false;
                    if xcolindexlen>1
                        xdata = self.dfref{rowindex,xcolindex(i)};
                        if self.isContourf && self.contourf_kws.enabled; kde2dComputationNeeded = true; end
                        if self.isContour3 && self.contour3_kws.enabled; kde2dComputationNeeded = true; end
                        if self.isContour && self.contour_kws.enabled; kde2dComputationNeeded = true; end
                    end
                    if (self.isHistogram2 && self.histogram2_kws.enabled) ...
                    || (self.isContourf && self.contourf_kws.enabled) ...
                    || (self.isContour3 && self.contour3_kws.enabled) ...
                    || (self.isContour && self.contour_kws.enabled)
                        if ycolindexlen>1
                            ydata = self.dfref{rowindex,ycolindex(i)};
                            if self.isContourf && self.contourf_kws.enabled; kde2dComputationNeeded = true; end
                            if self.isContour3 && self.contour3_kws.enabled; kde2dComputationNeeded = true; end
                            if self.isContour && self.contour_kws.enabled; kde2dComputationNeeded = true; end
                        end
                    end
                    if kde2dComputationNeeded
                        [bandwidth,density,crdx,crdy] = kde2d([xdata(:),ydata(:)],gridSize);
                    end

                    if self.legend_kws.enabled
                        if (self.isHistogram2 && self.histogram2_kws.enabled) ...
                        ||  (self.isContourf && self.contourf_kws.enabled) ...
                        ||  (self.isContour3 && self.contour3_kws.enabled) ...
                        ||  (self.isContour && self.contour_kws.enabled)
                            if xcolindexlen<2 && ycolindexlen>=1
                                lglabels = [ lglabels , ycolnames(i) ];
                            elseif xcolindexlen>1 && ycolindexlen<2
                                lglabels = [ lglabels , xcolnames(i) ];
                            else
                                lglabels = [ lglabels , xcolnames(i)+"-"+ycolnames(i) ];
                            end
                        else
                            lglabels = [ lglabels , xcolnames(i) ];
                        end
                    end

                    % add plot

                    if self.isHistogram && self.histogram_kws.enabled
                        self.currentFig.histogram = histogram   ( xdata ...
                                                                , histogram_kws_cell{:} ...
                                                                );
                        hold on;
                    elseif self.isHistogram2 && self.histogram2_kws.enabled
                        self.currentFig.histogram2 = histogram2 ( xdata ...
                                                                , ydata ...
                                                                , histogram2_kws_cell{:} ...
                                                                );
                        hold on;
                    elseif self.isContourf && self.contourf_kws.enabled
                        self.currentFig.contourf = struct();
                        self.currentFig.contourf.matrix = [];
                        self.currentFig.contourf.object = [];
                        [ self.currentFig.contourf.matrix ...
                        , self.currentFig.contourf.object ] = contourf  ( crdx ...
                                                                        , crdy ...
                                                                        , density ...
                                                                        , self.contourf_kws.levels ...
                                                                        , contourf_kws_cell{:} ...
                                                                        );
                        hold on;
                    elseif self.isContour3 && self.contour3_kws.enabled
                        self.currentFig.contour3 = struct();
                        self.currentFig.contour3.matrix = [];
                        self.currentFig.contour3.object = [];
                        [ self.currentFig.contour3.matrix ...
                        , self.currentFig.contour3.object ] = contour3  ( crdx ...
                                                                        , crdy ...
                                                                        , density ...
                                                                        , self.contour3_kws.levels ...
                                                                        , contour3_kws_cell{:} ...
                                                                        );
                        hold on;
                    elseif self.isContour && self.contour_kws.enabled
                        self.currentFig.contour = struct();
                        self.currentFig.contour.matrix = [];
                        self.currentFig.contour.object = [];
                        [ self.currentFig.contour.matrix ...
                        , self.currentFig.contour.object ] = contour( crdx ...
                                                                    , crdy ...
                                                                    , density ...
                                                                    , self.contour_kws.levels ...
                                                                    , contour_kws_cell{:} ...
                                                                    );
                        hold on;
                    elseif self.isHistfit && self.histfit_kws.enabled
                        self.currentFig.histogram = histfit ( xdata ...
                                                            , histfit_kws_cell{:} ...
                                                            );
                        hold on;
                    end

                end % loop plot

                self.currentFig.gca = gca;

            end

            if self.isHistogram2 || self.isContourf || self.isContour3 || self.isContour
                if ~getVecLen(self.colormap)
                    if self.isContour3 || self.isContour
                        self.colormap = "parula";
                    else
                        self.colormap = flipud(cold());
                    end
                end
                colormap(self.currentFig.gca,self.colormap);
            end

            % add axis labels

            if xcolindexlen>1
                self.currentFig.xlabel = xlabel("Variable Values", "Interpreter", "none");
            else
                self.currentFig.xlabel = xlabel(xcolnames(1), "Interpreter", "none");
            end

            if self.isHistogram2 || self.isContourf || self.isContour3 || self.isContour
                if ycolindexlen>1
                    self.currentFig.ylabel = ylabel("Variable Values", "Interpreter", "none");
                else
                    self.currentFig.ylabel = ylabel(ycolnames(1), "Interpreter", "none");
                end
                self.currentFig.zlabel = zlabel("Count", "Interpreter", "none");
            else
                self.currentFig.ylabel = ylabel("Count", "Interpreter", "none");
            end

            % add line colorbar

            colorbarEnabled = (self.isHistogram2 || self.isContourf || self.isContour3 || self.isContour) && self.colorbar_kws.enabled;
            if colorbarEnabled
                if isempty(self.colorbar_kws.fontsize) || ~isa(self.colorbar_kws.fontsize,"numeric")
                    self.colorbar_kws.fontsize = self.currentFig.ylabel.FontSize;
                end
                colorbar_kws_cell = convertStruct2Cell(self.colorbar_kws,{"enabled","singleOptions"});
                self.currentFig.colorbar = colorbar(colorbar_kws_cell{:});
                ylabel(self.currentFig.colorbar,"Density of Points","fontsize",self.colorbar_kws.fontsize, "Interpreter", "none");
            else
                colorbar('off');
                self.currentFig.colorbar = [];
            end

            self.doBasePlotStuff(self.legend_kws.enabled,lglabels);

            hold off;

        end % function plot

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % classdef HistPlot