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
%   AutoCorr_class  ( dataFrame ...
%                   , columns ...
%                   , rows ...
%                   , Err ...
%                   )
%
%   This is the AutoCorr_class class for generating objects
%   containing information about autocorrelation of the input data.
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
%           to be manipulated or be provided by the user.
%
%       columns
%
%           an array of strings or numbers corresponding to column names or indices 
%           of the input dataFrame for which the autocorrelation is computed.
%           This is a low-level internal argument and is not meant
%           to be accessed or be provided by the user.
%
%           If the input value is empty, the default will be the names of all columns 
%           of the input dataFrame.
%
%           Example usage:
%
%               1.  columns = [7,8,9]
%               2.  columns = ["SampleLogFunc","SampleVariable1"]
%               3.  columns = {"SampleLogFunc",9,"SampleVariable1"}
%               4.  columns = 7:9      # every column in the data frame starting from column #7 to #9
%               5.  columns = 7:2:20   # every other column in the data frame starting from column #7 to #20
%
%       rows
%
%           a numeric vector that represents the rows of the dataFrame that have been used or will be used 
%           to compute the autocorrelations. It can be either: 
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
%           an object of class Err_class for error reporting and warnings.
%
%   Attributes
%   ----------
%
%       df
%
%           a MATLAB data Table that contains the computed autocorrelations of the input 
%           dataFrame (MATLAB Table). This is a low-level internal argument and is not meant
%           to be manipulated or be provided by the user.
%
%       columns
%
%           optional property that determines the columns of the dataFrame for which the 
%           autocorrelation must be computed. It can have multiple forms:
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
%       bounds
%
%           a MATLAB data Table with two rows (`lowerLimit`, `upperLimit`) and with as many 
%           columns as specified by the property `columns`. The values of `bounds` are determined 
%           upon computing the autocorrelations of the input chain. The values represent the 1-sigma 
%           lower and upper standard errors on the computed autocorrelations. Any autocorrelation 
%           value bound within these limits can be considered random fluctuations at 1-sigma confidence level.
%
%       rows
%
%           a numeric vector that represents the rows of the dataFrame that have been used or will be used 
%           to compute the autocorrelations. It can be either: 
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
%           a structure containing several plotting tools for visualization of the 
%           computed autocorrelations as reported in the component `df`.
%
%   Superclass Attributes
%   ----------------------
%
%       See the documentation for the BasePlot class
%
%   Returns
%   -------
%
%       an object of AutoCorr_class class
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
classdef AutoCorr_class < dynamicprops

    properties(Access = public)
        columns
        bounds
        rows
        df
    end

    properties(Hidden)
        plotTypeList = ["line","scatter","lineScatter"];
        offset
        dfref
        Err
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access=public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = AutoCorr_class  ( dataFrame ...
                                        , columns ...
                                        , rows ...
                                        , Err ...
                                        )
            self.Err = Err;
            self.dfref = dataFrame;
            self.columns = columns;
            self.rows = []; if ~isempty(rows); self.rows = rows; end
            self.df = [];
            self.offset = 2;

            prop="plot"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            self.plot = struct();

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
                    methodList = ["plot","helpme","get","addBounds"];
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
            %   Computes the autocorrelation (ACF) for the selected columns of the object's dataFrame.
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
            %       get("columns",[8,9,13]) % compute the ACF only for the column numbers 8, 9, 13 in the input dataFrame
            %       get("columns",7:10)     % compute the ACF only for the column numbers 7, 8, 9, 10 in the input dataFrame
            %       get("rows", 1:5:10000)  % refine the input data every other 5 points, then compute the ACF
            %
            parseArgs(self,varargin{:});

            % check columns presence

            if getVecLen(self.columns)
                [colnames, colindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.columns);
            else
                colindex = 1:length(self.dfref.Properties.VariableNames);
                colnames = string(self.dfref.Properties.VariableNames);
            end

            % check rows presence

            if getVecLen(self.rows)
                rowindex = self.rows;
            else
                rowindex = 1:1:length(self.dfref{:,1});
            end

            % compute the autocorrelations

            nvar = length(colindex);
            nlag = length(rowindex);
            acf = zeros(nlag,nvar);
            bounds_array = zeros(2,nvar);

            if ~license('test','Econometrics_toolbox')
                warning ( newline ...
                        + "MATLAB Econometrics_toolbox is missing on your system. You can download the toolbox from: " + newline ...
                        + newline ...
                        + "    " + href("https://www.mathworks.com/products/econometrics.html") + newline ...
                        + newline ...
                        + "No autocorrelations will be computed." ...
                        + newline ...
                        );
                return
            end

            if nlag>1
                for i = 1:length(colindex)
                    [acf(:,i),lags,bounds_array(:,i)] = autocorr( self.dfref{rowindex,colindex(i)}, nlag-1 );
                end
                self.df = array2table([lags,acf]);
                colnames = "ACF_" + colnames;
                self.df.Properties.VariableNames = ["Lag", colnames];
                self.df.Properties.RowNames = string(lags);
                self.bounds = array2table(flipud(bounds_array));
                self.bounds.Properties.VariableNames = colnames;
                self.bounds.Properties.RowNames = ["lowerLimit","upperLimit"];
            else
                warning("The sample size is less than 2. No autocorrelations will be computed.");
                self.df = array2table(NaN(1,nvar));
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% graphics
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for plotType = self.plotTypeList
                self.resetPlot(plotType,"hard");
                %self.plot.(plotType).target.value = [1 ]
            end
            self.plot.reset = @self.resetPlot;
            self.plot.helpme = @self.helpme;
            self.plot.addBounds = @self.addBounds;

        end % get

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function addBounds(self,plotType)
            %
            %   Add standard error bounds of the computed autocorrelations to the user-input plot 
            %   (which must already exist and be active).
            %
            %   Parameters
            %   ----------
            %
            %       plotType
            %
            %           A string representing the name of the plot to which the bounds have to be added.
            %           Possible choices include:
            %
            %               "line", "scatter", "lineScatter"
            %
            %           Note that under the hood, this method is basically using the Target_class object of the requested 
            %           plotType to generate horizontal lines representing the standard error bounds of the autocorrelations.
            %           Therefore, if any visualization properties of the bounding lines needs to be changed, you can set
            %           them by changing the `hline_kws` and `scatter_kws` properties of `target` object of the input plotType.
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
            %       addBounds("line") % add standard error bounds to the already-existing line plot of the computed autocorrelations.
            %       addBounds("scatter") % add standard error bounds to the already-existing scatter plot of the computed autocorrelations.
            %
            currentFig = [];
            msg1= "The method addBounds() optionally takes an input argument ""plotType"" which is a string containing " + newline ...
                + "the type of the figure to which the bounds will be added. Possible plotType options are: " + newline + newline ...
                + "    " + strjoin(self.plotTypeList,", ") + newline;
            if nargin==1
                currentFig = get(groot,'CurrentFigure');
                if isempty(currentFig)
                    error   ( newline ...
                            + "Could not find any active figure to which the Autocorrelation Significance Bounds could be added." + newline ...
                            + msg1 ...
                            + newline ...
                            );
                end
            elseif nargin==2
                currentFig = self.plot.(plotType).currentFig.gcf;
                set(0, "CurrentFigure", currentFig);
            end
            currentAxes = currentFig.CurrentAxes;
            xLimits = currentAxes.XLim;
            boundListLen = length(self.bounds.Properties.VariableNames);
            xval = mean(xLimits);
            boundList = zeros(2*boundListLen,2);
            for i = 1:boundListLen
                indx = 2*(i-1)+1;
                boundList(indx,:) = [xval self.bounds{1,i}];
                boundList(indx+1,:) = [xval self.bounds{2,i}];
            end
            self.plot.(plotType).target.vline_kws.enabled = false;
            self.plot.(plotType).target.scatter_kws.enabled = false;
            self.plot.(plotType).target.hline_kws.enabled = true;
            self.plot.(plotType).target.hline_kws.color = "black";%[150 150 150 150]/255;
            self.plot.(plotType).target.xlimits = [xLimits(1)+0.001*xLimits(1), xLimits(2)];
            self.plot.(plotType).target.value = boundList;
            self.plot.(plotType).target.plot();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function resetPlot(self,varargin)
            %
            %   Reset the properties of the plot to the original default settings.
            %   Use this method when you change many attributes of the plot
            %   and you want to clean up
            %
            %   Parameters
            %   ----------
            %
            %       plotNames
            %
            %           An optional string or array of string values representing the names of plots to reset.
            %           If no value is provided, then all plots will be reset.
            %
            %       resetType
            %
            %           An optional string with possible value of "hard".
            %           If provided, the plot object will be regenerated from scratch.
            %           This includes reading the original data frame again and reseting everything.
            %
            %   Returns
            %   -------
            %
            %       None. 
            %
            %   Example
            %   -------
            %
            %       reset("line") % reset line plot to the default dettings
            %       reset("line","hard") % regenerate line plot from scratch
            %       reset(["line","scatter"],"hard") % regenerate line and scatter plots from scratch
            %       reset("hard") % regenrate all plots from scratch
            %
            requestedPlotTypeList = [];

            if nargin==1
                requestedPlotTypeList = self.plotTypeList;
            else
                for requestedPlotTypeCell = varargin{1}
                    if isa(requestedPlotTypeCell,"cell")
                        requestedPlotType = string(requestedPlotTypeCell{1});
                    else
                        requestedPlotType = string(requestedPlotTypeCell);
                    end
                    plotTypeNotFound = true;
                    for plotTypeCell = self.plotTypeList
                        plotType = string(plotTypeCell{1});
                        if strcmp(plotType,requestedPlotType)
                            requestedPlotTypeList = [ requestedPlotTypeList , plotType ];
                            plotTypeNotFound = false;
                            break;
                        end
                    end
                    if plotTypeNotFound
                        error   ( newline ...
                                + "The input plot-type argument, " + varargin{1} + ", to the resetPlot method" + newline ...
                                + "did not match any plot type. Possible plot types include:" + newline ...
                                + "    " + strjoin(self.plotTypeList,", ") + newline ...
                                );
                    end
                end
            end

            resetTypeIsHard = false;
            if nargin==3 && strcmpi(varargin{2},"hard")
                resetTypeIsHard = true;
                msgPrefix = "creating the ";
                msgSuffix = " plot object from scratch...";
            else
                msgPrefix = "reseting the properties of the ";
                msgSuffix = " plot...";
            end

            self.Err.marginTop = 0;
            self.Err.marginBot = 0;

            for requestedPlotTypeCell = requestedPlotTypeList

                self.Err.msg = msgPrefix + requestedPlotType + msgSuffix;
                self.Err.note();

                requestedPlotType = string(requestedPlotTypeCell);
                requestedPlotTypeLower = lower(requestedPlotType);
                plotName = "";

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                is3d = false;
                if contains(requestedPlotTypeLower,"3")
                    is3d = true;
                end

                % line

                if strcmp(requestedPlotTypeLower,"line") || strcmp(requestedPlotTypeLower,"line3")
                    plotName = "line"; if is3d; plotName = plotName + "3"; end
                    if resetTypeIsHard
                        self.plot.(plotName) = LineScatterPlot( self.df, plotName );
                    else
                        self.plot.(plotName).reset();
                    end
                    self.plot.(plotName).xcolumns = string(self.df.Properties.VariableNames(self.offset-1));
                    self.plot.(plotName).ycolumns = string(self.df.Properties.VariableNames(self.offset:end));
                    self.plot.(plotName).surface_kws.enabled = false;
                    self.plot.(plotName).ccolumns = "SampleLogFunc";
                    self.plot.(plotName).gca_kws.xscale = "linear";
                    self.plot.(plotName).plot_kws.linewidth = 1;
                    self.plot.(plotName).surface_kws.linewidth = 1;
                end

                % scatter / scatter3

                if strcmp(requestedPlotTypeLower,"scatter") || strcmp(requestedPlotTypeLower,"scatter3")
                    plotName = "scatter"; if is3d; plotName = plotName + "3"; end
                    if resetTypeIsHard
                        self.plot.(plotName) = LineScatterPlot( self.df, plotName );
                    else
                        self.plot.(plotName).reset();
                    end
                    self.plot.(plotName).ccolumns = string(self.df.Properties.VariableNames(self.offset));
                    self.plot.(plotName).xcolumns = string(self.df.Properties.VariableNames(self.offset-1));
                    self.plot.(plotName).ycolumns = string(self.df.Properties.VariableNames(self.offset:end));
                    self.plot.(plotName).gca_kws.xscale = "linear";
                    self.plot.(plotName).scatter_kws.size = 10;
                end

                % lineScatter / lineScatter3

                if strcmp(requestedPlotTypeLower,"linescatter") || strcmp(requestedPlotTypeLower,"linescatter3")
                    plotName = "lineScatter"; if is3d; plotName = plotName + "3"; end
                    if resetTypeIsHard
                        self.plot.(plotName) = LineScatterPlot( self.df, plotName );
                    else
                        self.plot.(plotName).reset();
                    end
                    self.plot.(plotName).surface_kws.enabled = false;
                    self.plot.(plotName).ccolumns = string(self.df.Properties.VariableNames(self.offset));
                    self.plot.(plotName).xcolumns = string(self.df.Properties.VariableNames(self.offset-1));
                    self.plot.(plotName).ycolumns = string(self.df.Properties.VariableNames(self.offset:end));
                    self.plot.(plotName).gca_kws.xscale = "linear";
                    if is3d
                        self.plot.(plotName).plot_kws.color = [200 200 200 75] / 255;
                    else
                        self.plot.(plotName).plot_kws.linewidth = 1;
                        self.plot.(plotName).plot_kws.color = uint8([200 200 200 100]);
                        self.plot.(plotName).scatter_kws.size = 20;
                    end
                end

                % 3d

                if is3d
                    self.plot.(plotName).xcolumns = string(self.df.Properties.VariableNames(self.offset));
                    self.plot.(plotName).ycolumns = string(self.df.Properties.VariableNames(self.offset+1));
                    self.plot.(plotName).zcolumns = string(self.df.Properties.VariableNames(self.offset));
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
