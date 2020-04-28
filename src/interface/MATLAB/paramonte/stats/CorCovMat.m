classdef CorCovMat < dynamicprops

    properties(Access = public)
        columns
        rows
        df
    end

    properties(Hidden)
        cormatPrecision
        matrixType
        isCorMat
        dfref
        title
        Err
    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods(Access=public)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function self = CorCovMat   ( dataFrame ...
                                    , columns ...
                                    , method ...
                                    , rows ...
                                    , Err ...
                                    )
            self.Err = Err;
            self.dfref = dataFrame;
            self.columns = columns;
            self.cormatPrecision = 2;
            if ~isempty(method)
                prop="method"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                self.method = method;
            end
            self.columns = columns;
            self.rows = []; if ~isempty(rows); self.rows = rows; end
            self.df = [];

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
                rowindex = self.rows;
            else
                rowindex = 1:1:length(self.dfref{:,1});
            end

            % construct the matrix dataframe

            rowindexLen = length(rowindex);
            colnamesLen = length(colnames);
            if  rowindexLen > colnamesLen
                if  self.isCorMat
                    self.df = array2table( corr( self.dfref{rowindex,colnames} , "type" , self.method ) );
                else
                    self.df = array2table( cov( self.dfref{rowindex,colnames} ) );
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

            for plotType = ["heatmap"]
                self.resetPlot(plotType,"hard");
            end
            self.plot.reset = @self.resetPlot;
            self.plot.heatmap.title = self.title;
            if self.isCorMat; self.plot.heatmap.precision = self.cormatPrecision; end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % set title for the heatmap

            if isempty(self.title)
                if self.isCorMat
                    self.title = "The " + self.method + "'s Correlation Strength Matrix";
                else
                    self.title = "The Covariance Matrix";
                end
            end
            try
                self.plot.heatmap.title = self.title;
                if self.isCorMat; self.plot.heatmap.precision = self.cormatPrecision; end
            end

        end % get

        %***************************************************************************************************************************
        %***************************************************************************************************************************

    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Access = public, Hidden)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function resetPlot(self,varargin)

            requestedPlotTypeList = [];
            plotTypeList = ["heatmap"];

            if nargin==1
                requestedPlotTypeList = plotTypeList;
            else
                for requestedPlotTypeCell = varargin{1}
                    if isa(requestedPlotTypeCell,"cell")
                        requestedPlotType = string(requestedPlotTypeCell{1});
                    else
                        requestedPlotType = string(requestedPlotTypeCell);
                    end
                    plotTypeNotFound = true;
                    for plotTypeCell = plotTypeList
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
                                + "line, lineScatter." + newline ...
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

                % heatmap

                if strcmp(requestedPlotTypeLower,"heatmap")
                    plotName = "heatmap";
                    if resetTypeIsHard
                        self.plot.(plotName) = HeatmapPlot( self.df );
                    else
                        self.plot.(plotName).reset();
                    end
                    if self.isCorMat
                        self.plot.(plotName).heatmap_kws.ColorLimits = [-1 1];
                    end
                    self.plot.(plotName).xcolumns = self.df.Properties.VariableNames;
                    self.plot.(plotName).ycolumns = self.df.Properties.RowNames;
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