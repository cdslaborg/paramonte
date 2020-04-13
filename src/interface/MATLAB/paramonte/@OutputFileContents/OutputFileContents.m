classdef OutputFileContents < dynamicprops
%   This is the OutputFileContents class for generating instances 
%   of ParaMonte output file contents. The ParaMonte read* methods
%   return an object or a list of objects of class OutputFileContents.
%
%   Constructor Parameters
%   ----------------------
%       file
%           full path to the file containing the sample/chain.
%       delimiter
%           the delimiter used in the sample/chain file, which 
%           must be provided by the user.
%       fileType
%           The type of the file to be parsed:
%           "chain", "sample", "progress", "report", "restart"
%       mpiEnabled
%           boolean value indicating weather the code is being run in parallel
%
%   Attributes
%   ----------
%
%   The attributes of the constructed object depend on the type of the file being parsed:
%
%       Sample/Chain files
%       ------------------
%       file
%           full path to the file containing the sample/chain.
%       delimiter
%           the delimiter used in the sample/chain file, which 
%           must be provided by the user.
%       ndim
%           the inferred number of dimensions of the domain of the 
%           objective function for which the simulation has been performed.
%       count
%           number of points (states) in the sample/chain file. 
%           This is essentially, the number of rows in the file 
%           minus one (representing the header line).
%       [df]
%           if the input file contents is structured in a format that
%           could be read as a MATLAB Table, then the contents of the file
%           will be stored in the form of a MATLAB Table in this property, 
%           hence called 'df'.
%
%   Returns
%   -------
%       OutputFileContents
%
%   ----------------------------------------------------------------------
%   """

    properties(Access = public)
        file = [];
        delimiter = [];
        %count = [];
        %stats = [];
        %ndim = [];
        %df = [];
    end

    properties(Hidden)
        offset = [];
    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Access = public)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function self = OutputFileContents  ( file ...
                                            , fileType ...
                                            , delimiter ...
                                            , methodName ...
                                            , mpiEnabled ...
                                            , markovChainRequested ...
                                            , Err ...
                                            )

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.file = file;
            self.delimiter = delimiter;

            Err.marginTop = 0;
            Err.marginBot = 0;
            Err.msg = "reading file contents... "; Err.note();
            timer = Timer_class();
            timer.tic()

            self.file
            file
            d = importdata  ( self.file ...
                            ..., "delimiter", self.delimiter ...
                            );
            colheadersLen = length(d.colheaders);
            for icol = 1:colheadersLen
                if strcmp(d.colheaders{icol},"SampleLogFunc")
                    break
                end
            end
            self.offset = icol + 1; % index of the first variable
            self.addprop("ndim");
            self.ndim   = colheadersLen - self.offset + 1;
            self.addprop("count");
            self.count  = length(d.data(:,1));

            if markovChainRequested
                cumSumWeight = cumsum(d.data(:,self.offset-2));
                if cumSumWeight(end) ~= self.count % it is indeed a compact chain
                    dMarkov = zeros( cumSumWeight(end) , self.ndim + self.offset );
                    istart = 1;
                    for irow = 1:self.count
                        iend = cumSumWeight(irow);
                        dMarkov(istart:iend,:) = d.data(irow,:);
                        istart = iend + 1;
                    end
                    d.data = dMarkov;
                    self.count = cumSumWeight(end);
                end
            end

            self.addprop("df");
            self.df = array2table(d.data,'VariableNames',d.colheaders);

            updateUser([]);

            if mpiEnabled
                Err.marginTop = 0;
                Err.marginBot = 1;
                Err.msg = "ndim = " + string(self.ndim) + ", count = " + string(self.count);
                Err.note();
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% statistics
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.addprop("stats");
            self.stats = struct();

            % add chain cormat

%            self.stats.cormat = _pm.stats.CorMat( dataFrame     = self.df
%                                                , columns       = range(self.offset,self.offset+self.ndim)
%                                                , method        = "pearson"
%                                                )
%            
%            timer.tic( msg = "computing sample correlation matrix... " )
%            self.stats.cormat()
%            timer.toc()
%
%            # add chain covmat
%
%            self.stats.covmat = _pm.stats.CovMat( dataFrame     = self.df
%                                                , columns       = range(self.offset,self.offset+self.ndim)
%                                                )
%
%            timer.tic( msg = "computing sample covariance matrix... " )
%            self.stats.covmat()
%            timer.toc()
%
%            self.stats.maxLogFunc = _pm.stats.getMaxLogFunc(dataFrame = self.df)
%            #self.stats.maxLogFunc = _Struct()
%            #self.stats.maxLogFunc.idrow = self.df[["SampleLogFunc"]].idxmax().values[0]
%            #self.stats.maxLogFunc.value = self.df[["SampleLogFunc"]].iat[self.stats.maxLogFunc.idrow,0]
%            #self.stats.maxLogFunc.dfrow = self.df.iloc[self.stats.maxLogFunc.idrow,:]
%            #self.stats.maxLogFunc.state = self.df.iloc[self.stats.maxLogFunc.idrow,self.offset:]
%
%            # add chain autocorrelation
%
%            self.stats.acf = _pm.stats.AutoCorr ( dataFrame     = self.df
%                                                , columns       = range(self.offset-1,self.offset+self.ndim)
%                                                )
%
%            timer.tic( msg = "computing autocorrelations... " )
%            self.stats.acf()
%            timer.toc()
%
%            ############################################################################################################################
%            #### graphics
%            ############################################################################################################################
%
%            timer.tic( msg = "adding graphics tools... " )
%
%            # add HistPlot
%
%            self.plot = _Struct()
%
%            self.plot.hist = _pm.vis.HistPlot   ( dataFrame = self.df
%                                                , columns = self.df.columns[self.offset:]
%                                                )
%
%            # add LinePlot
%
%            self.plot.line = _pm.vis.LinePlot   ( dataFrame = self.df
%                                                , ycolumns = self.df.columns[self.offset:]
%                                                , ccolumns = "SampleLogFunc"
%                                                , lc_kws =  {
%                                                            #"linewidth":0.75,
%                                                            #"cmap":"viridis",
%                                                            "cmap":"autumn",
%                                                            #"alpha":0.5,
%                                                            }
%                                                , colorbar_kws =    {
%                                                                    "extend":"neither",
%                                                                    "orientation":"vertical",
%                                                                    #"spacing":"uniform",
%                                                                    }
%                                                #, legend_kws = None
%                                                )
%
%            # add ScatterPlot
%
%            self.plot.scatter = _pm.vis.ScatterPlot ( dataFrame = self.df
%                                                    , ycolumns = self.df.columns[self.offset:]
%                                                    , ccolumns = "SampleLogFunc"
%                                                    #, scatter_kws = {}
%                                                    , colorbar_kws =    {
%                                                                        "extend":"neither",
%                                                                        "orientation":"vertical",
%                                                                        #"spacing":"uniform",
%                                                                        }
%                                                    #, legend_kws = None
%                                                    )
%
%            # add DensityMapPlot
%
%            xindex = self.offset
%            yindex = self.offset + 1
%            if self.ndim==1: xindex, yindex = yindex-1, xindex-1
%                
%            self.plot.density = _pm.vis.DensityMapPlot  ( dataFrame = self.df
%                                                        , xcolumn = xindex
%                                                        , ycolumn = yindex
%                                                        )
%            #print(self.stats.maxLogFunc.idrow,xindex,yindex)
%            #print(self.df)
%            #print(self.df.iat[self.stats.maxLogFunc.idrow,xindex])
%            #print(self.df.iat[self.stats.maxLogFunc.idrow,yindex])
%            self.plot.density.target.__init__   ( value = [ self.df.iat[self.stats.maxLogFunc.idrow,xindex], self.df.iat[self.stats.maxLogFunc.idrow,yindex] ]
%                                                , scatter_kws = {"label":"maxLogFunc"}
%                                                )
%
%            # add GridPlot
%
%            endColindex = _np.min( [self.offset+3, self.offset+self.ndim] )
%            self.plot.grid = _pm.vis.GridPlot   ( dataFrame = self.df
%                                                , columns = self.df.columns[self.offset:endColindex]
%                                                , scatterplot_kws = {"ccolumns": "SampleLogFunc"}
%                                                , _methodName = _pm.names.paradram
%                                                )
%
%            # add ScatterLinePlot
%
%            self.plot._scatterline = _pm.vis.ScatterLinePlot( dataFrame = self.df
%                                                            , ycolumns = self.df.columns[self.offset:]
%                                                            , lccolumns = "SampleLogFunc"
%                                                            , lc_kws =  {
%                                                                        #"linewidth":0.75,
%                                                                        #"cmap":"viridis",
%                                                                        "cmap":"autumn",
%                                                                        #"alpha":0.5,
%                                                                        }
%                                                            #, scatter_kws = {}
%                                                            , colorbar_kws =    {
%                                                                                "extend":"neither",
%                                                                                "orientation":"vertical",
%                                                                                #"spacing":"uniform",
%                                                                                }
%                                                            #, legend_kws = None
%                                                            )
%
%            timer.toc()

            function updateUser(msg)
                if isempty(msg)
                    timer.toc();
                    Err.msg = "done in " + sprintf("%.6f",string(timer.delta)) + " seconds."; Err.note();
                else
                    Err.msg = msg; Err.note();
                end
            end

        end % constructor

        %***************************************************************************************************************************
        %***************************************************************************************************************************

    end % methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Hidden, Static)
        % These methods have been implemented to override the default 'handle' class methods,
        % so that they won't pop-up after pressing 'Tab' button.
        function addlistener    ();  end
        function delete         ();  end
        function findobj        ();  end
        function findprop       ();  end
        function valid          ();  end
        function listener       ();  end
        function notify         ();  end
    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

end % classdef OutputFileContents < handle