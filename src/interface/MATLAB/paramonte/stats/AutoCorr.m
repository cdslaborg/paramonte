classdef AutoCorr < dynamicprops
%   This is the class for generating object of type CorMat which, 
%   upon construction, will provide methods to compute and plot the 
%   autocorrelations of the selected columns of the input dataFrame.
%
%   Parameters
%   ----------
%       dataFrame
%           a Pandas dataframe based upon the selected comlumns of which 
%           the autocorrelations will be computed.
%       columns
%           optional argument that determines the columns of the input dataFrame to be 
%           used in the computation of the autocorrelations. It can have three forms:
%               1. a list of column indices from the input dataFrame.
%               2. a list of column names from dataFrame.columns.
%               3. a range(start,stop,step), representing the column indices in dataFrame.
%           Examples:
%               1. columns = [0,1,4,3]
%               2. columns = ["SampleLogFunc","SampleVariable1"]
%               3. columns = range(17,7,-2)
%           If not provided, the default behavior includes all columns of the dataFrame.
%       rows
%           optional argument that determines the rows of the input dataFrame to be 
%           used in the computation of the autocorrelations. It can be either:
%               1. a range(start,stop,step), or, 
%               2. a list of row indices from dataFrame.index.
%           Examples:
%               1. rows = range(17,7,-2)
%               2. rows = [i for i in range(7,17)]
%           If not provided, the default behavior includes all rows of the dataFrame.
%
%   Attributes
%   ----------
%       All of the parameters described above, except dataFrame.
%           a reference to the dataFrame will be implicitly stored in the object.
%       df
%           a pandas dataframe containing the computed autocorrelations.
%       plot
%           a structure containing the following plotting tools:
%               heatmap
%                   a callable object of class HeatMap which will enable 
%                   plotting of the correlation matrix.
%
%   Returns
%   -------
%       corMat
%           an object of type class CorMat.
%
%   ----------------------------------------------------------------------


    properties(Access = public)
        columns
        bounds
        rows
        df
    end

    properties(Hidden)
        plotTypeList = ["line","scatter","lineScatter"];
        cormatPrecision
        matrixType
        isCorMat
        offset
        dfref
        title
        Err
    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods(Access=public)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function self = AutoCorr    ( dataFrame ...
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

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function get(self,varargin)
        %   compute the correlation matrix of the selected columns 
        %   of the input dataframe to the object's constructor.
        %   
        %   Parameters
        %   ----------
        %       None
        %   
        %   Returns
        %   -------
        %       None. However, this method causes side-effects by manipulating 
        %       the existing attributes of the object.

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
            nlag = length(rowindex)-1;
            acf = zeros(nlag,nvar);
            bounds_array = zeros(2,nvar);

            for i = 1:length(colindex)
                [acf(:,i),lags,bounds_array(:,i)] = autocorr( self.dfref{rowindex,colindex(i)}, nlag-1 );
            end
            self.df = array2table([lags,acf]);
            colnames = "ACF_" + colnames;
            self.df.Properties.VariableNames = ["Lag", colnames];
            self.df.Properties.RowNames = string(lags);
            self.bounds = array2table(bounds_array);
            self.bounds.Properties.VariableNames = colnames;
            self.bounds.Properties.RowNames = ["upperLimit","lowerLimit"];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% graphics
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for plotType = self.plotTypeList
                self.resetPlot(plotType,"hard");
                %self.plot.(plotType).target.value = [1 ]
            end
            self.plot.reset = @self.resetPlot;
            self.plot.addBounds = @self.addBounds;

        end % get

        %***************************************************************************************************************************
        %***************************************************************************************************************************

    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Access = public, Hidden)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function addBounds(self,plotType)
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

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function resetPlot(self,varargin)

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
                msgSuffix = " plot object from scrach...";
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

        %***************************************************************************************************************************
        %***************************************************************************************************************************

    end % methods

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

end
