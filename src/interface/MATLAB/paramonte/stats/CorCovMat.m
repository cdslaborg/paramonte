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
%   CorCovMat   ( dataFrame ...
%               , columns ...
%               , method ...
%               , rows ...
%               , Err ...
%               , reportEnabled ...
%               )
%
%   This is the CorCovMat class for generating objects containing the
%   computed correlation or covariance matrices of the input data.
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
%           A MATLAB data Table from which the selected data will be plotted.
%           This is a low-level internal argument and is not meant
%           to be manipulated or be provided by the user.
%
%       columns
%
%           An array of strings or numbers corresponding to column names or indices
%           of the input dataFrame for which the correlation/covariance matrix is
%           computed. This is a low-level internal argument and is not meant
%           to be accessed or be provided by the user.
%
%           If the input value is empty, the default will be the names of all
%           columns of the input dataFrame.
%
%           Example usage:
%
%               1.  columns = [7,8,9]
%               2.  columns = ["SampleLogFunc","SampleVariable1"]
%               3.  columns = {"SampleLogFunc",9,"SampleVariable1"}
%               4.  columns = 7:9    # every column in dataFrame from column #7 to #9
%               5.  columns = 7:2:20 # every other column in dataFrame from column #7 to #20
%
%       method
%
%           A string or char vector with one of the following possible values:
%
%               "pearson"   : compute the Pearson's correlation matrix of the input data
%               "kendall"   : compute the Kendall's correlation matrix of the input data
%               "spearman"  : compute the Spearman's correlation matrix of the input data
%
%           If an empty object is provided as input, the covariance matrix
%           of the input dataFrame will be computed for the selected columns.
%           This is a low-level internal argument and is not meant
%           to be accessed or be provided by the user.
%
%       rows
%
%           A numeric vector that represents the rows of the dataFrame that have been used
%           or will be used to compute the correlation/covariance matrix. It can be either:
%
%               1.  a numeric range, or,
%               2.  a list of row indices of the dataFrame.
%
%           Example usage:
%
%               1.  rows = 10000:-2:3
%               2.  rows = [12,46,7,8,9,4,7,163]
%
%           If not provided, the default includes all rows of the input dataFrame.
%
%       Err
%
%           An object of class Err_class for error reporting and warnings.
%
%       reportEnabled
%
%           An boolean indicating whether any descriptive messages 
%           should be output or not.
%
%   Attributes
%   ----------
%
%       df
%
%           A MATLAB data Table that contains the computed correlation/covariance matrix of the input
%           dataFrame (MATLAB Table). This is a low-level internal argument and is not meant
%           to be manipulated or be provided by the user.
%
%       columns
%
%           Optional property that determines the columns of the dataFrame for which the
%           correlation/covariance matrix must be computed. It can have multiple forms:
%
%               1.  a numeric or cell array of column indices in the input dataFrame.
%               2.  a string or cell array of column names in dataFrame.Properties.VariableNames.
%               3.  a cell array of a mix of the above two.
%               4.  a numeric range.
%
%           Example usage:
%
%               1.  columns = [7,8,9]
%               2.  columns = ["SampleLogFunc","SampleVariable1"]
%               3.  columns = {"SampleLogFunc",9,"SampleVariable1"}
%               4.  columns = 7:9      # every column in the data frame starting from column #7 to #9
%               5.  columns = 7:2:20   # every other column in the data frame starting from column #7 to #20
%
%           The default value is the names of all columns of the input dataFrame.
%
%       method (available only in correlation matrix objects)
%
%           A string or char vector with one of the following possible values:
%
%               "pearson"   : compute the Pearson's correlation matrix of the input data
%               "kendall"   : compute the Kendall's correlation matrix of the input data
%               "spearman"  : compute the Spearman's correlation matrix of the input data
%
%       rows
%
%           A numeric vector that represents the rows of the dataFrame that have been used
%           or will be used to compute the correlation/covariance matrix. It can be either:
%
%               1.  a numeric range, or,
%               2.  a list of row indices of the dataFrame.
%
%           Example usage:
%
%               1.  rows = 15:-2:8
%               2.  rows = [12,46,7,8,9,4,7,163]
%
%           If not provided, the default includes all rows of the input dataFrame.
%
%       plot
%
%           A structure containing several plotting tools for visualization of the
%           computed correlation/covariance matrix as reported in the component ``df``.
%
%   Returns
%   -------
%
%       An object of ``CorCovMat`` class.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
classdef CorCovMat < dynamicprops

    properties(Access = public)
        columns
        plot
        rows
        df
    end

    properties(Hidden)
        cormatPrecision
        reportEnabled
        plotTypeList
        matrixType
        rowsindex
        isCorMat
        dfref
        Err
        title
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access=public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = CorCovMat   ( dataFrame ...
                                    , columns ...
                                    , method ...
                                    , rows ...
                                    , Err ...
                                    , reportEnabled ...
                                    )
            self.reportEnabled = reportEnabled;
            self.Err = Err;
            self.dfref = dataFrame;
            self.columns = columns;
            self.cormatPrecision = 2;
            if ~isempty(method)
                prop="method"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                self.method = method;
            end
            self.columns = columns;
            self.title = struct();
            self.title.text = "";
            self.title.subtext = "";
            self.rows = []; if ~isempty(rows); self.rows = rows; end
            self.df = [];

            %prop="plot"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            self.plot = struct();
            self.plot.helpme = @self.helpme;
            self.get();

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
                if contains(varargin{1},"reset")
                    cmd = "doc self.resetPlot";
                    methodNotFound = false;
                else
                    methodList = ["plot","helpme","get"];
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

        function get(self,varargin)
            %
            %   Compute the correlation matrix of the selected columns of the object's dataFrame.
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
            %       get("columns",[8,9,13]) % compute only for the column numbers 8, 9, 13 in the input dataFrame
            %       get("columns",7:10)     % compute only for the column numbers 7, 8, 9, 10 in the input dataFrame
            %       get("rows", 1:5:10000)  % refine the input data every other 5 points, then compute the quantities.
            %
            parseArgs(self,varargin{:});

            if isprop(self,"method")
                self.isCorMat = true;
                self.matrixType = "correlation";
                if ~any(strcmp(["pearson","kendall","spearman"],self.method))
                    error   ( newline ...
                            + "The requested correlation type must be one of the following string values, " + newline + newline ...
                            + "    pearson  : standard correlation coefficient" + newline ...
                            + "    kendall  : Kendall Tau correlation coefficient" + newline ...
                            + "    spearman : Spearman rank correlation." + newline ...
                            + newline ...
                            );
                end
            else
                self.isCorMat = false;
                self.matrixType = "covariance";
            end

            % check columns presence

            if getVecLen(self.columns)
                [colnames, ~] = getColNamesIndex(self.dfref.Properties.VariableNames,self.columns); % colindex
            else
                %colindex = 1:length(self.dfref.Properties.VariableNames);
                colnames = string(self.dfref.Properties.VariableNames);
            end

            % check rows presence

            if getVecLen(self.rows)
                self.rowsindex = self.rows;
            else
                self.rowsindex = 1:1:length(self.dfref{:,1});
            end

            % construct the matrix dataframe

            rowindexLen = length(self.rowsindex);
            colnamesLen = length(colnames);
            if  rowindexLen > colnamesLen
                if  self.isCorMat
                    self.df = array2table( corr( self.dfref{self.rowsindex,colnames} , "type" , self.method ) );
                else
                    self.df = array2table( cov( self.dfref{self.rowsindex,colnames} ) );
                end
            else
                warning ( "The number of columns (" + string(colnamesLen) + ") is more than the number of rows (" + string(rowindexLen) + ") in the data-frame (Table). " ...
                        + "The " + self.matrixType + " matrix cannot be computed." ...
                        ...+ newline ...
                        );
                self.df = array2table( NaN(colnamesLen , colnamesLen) );
            end

            self.df.Properties.VariableNames = colnames;
            self.df.Properties.RowNames = colnames;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% graphics
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.plotTypeList = ["heatmap"];

            self.resetPlot("hard");
            self.plot.reset = @self.resetPlot;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end % get

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function resetPlot(self,resetType,plotNames)
            %
            %   Reset the properties of the plot to the original default settings.
            %   Use this method when you change many attributes of the plot and
            %   you want to clean up and go back to the default settings.
            %
            %   Parameters
            %   ----------
            %
            %       resetType
            %
            %           An optional string with possible value of "hard" or "soft".
            %           If set to "hard", the plot object(s) will be regenerated from scratch.
            %           This includes reading the original data frame again and resetting everything.
            %           If set to "soft", then only the parameters of the plot objects will be reset
            %           to the default values. The default is "soft".
            %
            %       plotNames
            %
            %           An optional string or array of string values representing the names of plots to reset.
            %           If no value is provided, then all plots will be reset. Note that if ``plotNames`` is
            %           present, then ``resetType`` must be also given as input argument.
            %
            %   Returns
            %   -------
            %
            %       None.
            %
            %   Example
            %   -------
            %
            %       reset() % reset all plots to the default settings
            %       reset("soft","line") % reset the line plot from scratch.
            %       reset("hard",["line","scatter"]) % regenerate line and scatter plots from scratch
            %       reset("hard") % regenerate all plots from scratch
            %
            if nargin<3 || isempty(plotNames); plotNames = "all"; end
            if nargin<2 || isempty(resetType); resetType = "soft"; end

            requestedPlotTypeList = [];
            if isstring(plotNames) || ischar(plotNames)
                plotTypeLower = lower(string(plotNames));
                if strcmp(plotTypeLower,"all")
                    requestedPlotTypeList = self.plotTypeList;
                elseif any(contains(self.plotTypeList,plotNames))
                    requestedPlotTypeList = [plotNames];
                else
                    self.reportWrongPlotName(plotNames);
                end
            elseif getVecLen(plotNames)
                for plotName = plotNames
                    if ~any(contains(self.plotTypeList,plotName))
                        self.reportWrongPlotName(plotName);
                    end
                end
            else
                self.reportWrongPlotName("a none-string none-list object.")
            end

            resetTypeIsHard = false;
            if isstring(resetType) || ischar(resetType)
                resetTypeIsHard = strcmp(lower(resetType),"hard");
            else
                self.Err.marginTop = 1;
                self.Err.marginBot = 1;
                self.Err.msg    = "The input argument ``resetType`` must be a string representing" + newline ...
                                + "the type of the reset to be performed on the plots." + newline ...
                                + "A list of possible plots includes: ""hard"", ""soft""" + newline ...
                                + "Here is the help for the ``reset()`` method: " + newline ...
                                + newline ...
                                + string(help("self.resetPlot")) ...
                                ;
                self.Err.abort();
            end

            lenVariableNames = length(self.df.Properties.VariableNames);

            if resetTypeIsHard
                msgPrefix = "creating the ";
                msgSuffix = " plot object from scratch...";
            else
                msgPrefix = "resetting the properties of the ";
                msgSuffix = " plot...";
            end

            self.Err.marginTop = 0;
            self.Err.marginBot = 0;

            for requestedPlotTypeCell = requestedPlotTypeList

                requestedPlotType = string(requestedPlotTypeCell);
                requestedPlotTypeLower = lower(requestedPlotType);

                isLine          = contains(requestedPlotTypeLower,"line");
                isScatter       = contains(requestedPlotTypeLower,"scatter");
                isHistfit       = contains(requestedPlotTypeLower,"histfit");
                isHistogram2    = contains(requestedPlotTypeLower,"histogram2");
                isHistogram     = contains(requestedPlotTypeLower,"histogram") && ~isHistogram2;
                isContourf      = contains(requestedPlotTypeLower,"contourf");
                isContour3      = contains(requestedPlotTypeLower,"contour3");
                isContour       = contains(requestedPlotTypeLower,"contour") && ~(isContourf || isContour3);
                isGridPlot      = contains(requestedPlotTypeLower,"grid");
                isHeatmap       = contains(requestedPlotTypeLower,"heatmap");
                is3d            = contains(requestedPlotTypeLower,"3") || isHistogram2;

                isLineScatterPlot = isLine || isScatter;
                isDensityPlot = isHistfit || isHistogram2 || isHistogram || isContourf || isContour3 || isContour;

                if self.reportEnabled
                    self.Err.msg = msgPrefix + requestedPlotType + msgSuffix;
                    self.Err.note();
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%% reset heatmap
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if isHeatmap

                    if resetTypeIsHard
                        self.plot.(requestedPlotType) = HeatmapPlot( requestedPlotType, self.df, @self.resetPlot );
                    else
                        self.plot.(requestedPlotType).resetInternal();
                    end
                    if self.isCorMat
                        self.plot.(requestedPlotType).heatmap.ColorLimits = [-1 1];
                    end
                    self.plot.(requestedPlotType).xcolumns = self.df.Properties.VariableNames;
                    self.plot.(requestedPlotType).ycolumns = self.df.Properties.RowNames;

                    if self.isCorMat; self.plot.heatmap.precision = self.cormatPrecision; end

                    % set title for the heatmap

                    if ~( (isstring(self.title.text) || ischar(self.title.text)) && getVecLen(self.title.text) )
                        if self.isCorMat
                            self.title.text = "The " + self.method + "'s Correlation Strength Matrix";
                        else
                            self.title.text = "The Covariance Matrix";
                        end
                    end
                    self.plot.heatmap.heatmap.kws.title = self.title.text;

                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end