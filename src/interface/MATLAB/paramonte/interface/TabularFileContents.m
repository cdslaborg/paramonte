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
%   This is the TabularFileContents class for generating instances
%   of ParaMonte output file contents. For example, the ParaDRAM
%   readSample(), readChain(), readMarkovChain(), readProgress()
%   methods return an object or a list of objects of class
%   TabularFileContents.
%
%   This is an internal ParaMonte class and it is
%   not meant to be directly accessible to the users.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
classdef TabularFileContents < OutputFileContents

    properties(Access = public)
        delimiter = [];
        %count = [];
        %stats = [];
        %ndim = [];
        %df = [];
    end

    properties(Hidden)
        sampleLogFuncColName = "";
        isProgressFile = false;
        plotTypeList = [];
        offset = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = TabularFileContents ( file ...
                                            , fileType ...
                                            , delimiter ...
                                            , methodName ...
                                            , reportEnabled ...
                                            , markovChainRequested ...
                                            , Err ...
                                            )

            self = self@OutputFileContents(file,methodName,reportEnabled,Err);
            self.delimiter = delimiter;
            self.timer.tic();
            if strcmp(fileType,"progress")
                self.isProgressFile = true;
            else
                self.sampleLogFuncColName = "SampleLogFunc";
            end

            if ~strcmp(fileType,"sample") && ~strcmp(fileType,"chain") && ~(self.isProgressFile || markovChainRequested)
                self.Err.msg    = "Internal error occurred. The input fileType is not recognized. " + newline ...
                                + "Please report this error on the issues page of the ParaMonte " + newline ...
                                + "library's GitHub repository page." ...
                                ;
                self.Err.marginTop = 1;
                self.Err.marginBot = 1;
                self.Err.abort();
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            d = importdata  ( self.file ...
                            ..., "delimiter", self.delimiter ...
                            );

            if isfield(d,"colheaders")
                prop="ncol"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                self.ncol = length(d.colheaders);
            else
                self.Err.marginTop = 1;
                self.Err.marginBot = 1;
                self.Err.msg    = "The structure of the file """ + self.file + """ does not match a " + methodName + " " + fileType + " file. " ...
                                + "Verify the contents of this file before attempting to read this file.";
                self.Err.abort();
            end

            if self.isProgressFile
                self.offset = self.ncol;
            else
                for icol = 1:self.ncol
                    if strcmp(d.colheaders{icol},self.sampleLogFuncColName)
                        break;
                    end
                end
                self.offset = icol + 1; % index of the first variable
                prop="ndim"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                self.ndim = self.ncol - self.offset + 1;
            end

            prop="count"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            self.count = length(d.data(:,1));

            if markovChainRequested
                cumSumWeight = cumsum(d.data(:,self.offset-2));
                if cumSumWeight(end) > self.count % it is indeed a compact chain
                    dMarkov = zeros( cumSumWeight(end) , self.ncol );
                    istart = 1;
                    for irow = 1:self.count
                        iend = cumSumWeight(irow);
                        for i = istart:iend
                            dMarkov(i,:) = d.data(irow,:);
                        end
                        istart = iend + 1;
                    end
                    d.data = dMarkov;
                    self.count = cumSumWeight(end);
                elseif cumSumWeight(end) < self.count % there is something wrong about this chain output file
                    self.Err.marginTop = 1;
                    self.Err.marginBot = 1;
                    self.Err.msg    = "The internal contents of the file """ + self.file + """ does not match a " + methodName + " " + fileType + "file. " ...
                                    + "In particular, the SampleWeight data column in this file appears to contain values that are either non-positive or non-integer.";
                    self.Err.abort();
                end
            end

            prop="df"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            self.df = array2table(d.data,'VariableNames',d.colheaders);

            self.updateUser([]);

            if self.reportEnabled && ~self.isProgressFile
                self.Err.marginTop = 0;
                self.Err.marginBot = 1;
                self.Err.msg = "ndim = " + string(self.ndim) + ", count = " + string(self.count);
                self.Err.note();
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% statistics
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if ~self.isProgressFile

                prop="stats"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end

                self.stats = struct();

                % add chain cormat

                self.updateUser("computing the sample correlation matrix...");
                self.stats.cormat = CorCovMat   ( self.df ...
                                                , self.offset:self.ncol ...
                                                , "pearson" ... method
                                                , [] ... rows
                                                , self.Err ...
                                                , self.reportEnabled ...
                                                );
                self.updateUser([]);

                % add chain covmat

                self.updateUser("computing the sample covariance matrix...");
                self.stats.covmat = CorCovMat   ( self.df ...
                                                , self.offset:self.ncol ...
                                                , [] ... method
                                                , [] ... rows
                                                , self.Err ...
                                                , self.reportEnabled ...
                                                );
                self.updateUser([]);

                % add chain autocorr

                self.updateUser("computing the sample autocorrelation...");
                self.stats.autocorr = AutoCorr_class( self.df ...
                                                    , self.offset-1:self.ncol ...
                                                    , [] ... rows
                                                    , self.Err ...
                                                    , reportEnabled ...
                                                    );
                self.updateUser([]);

                self.stats.maxLogFunc = getMaxLogFunc(self.df,self.sampleLogFuncColName);

            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% graphics
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.plotTypeList = [ "line" ...
                                , "scatter" ...
                                , "lineScatter" ...
                                ];

            if ~self.isProgressFile
                self.plotTypeList = [ self.plotTypeList ...
                                    , "line3" ...
                                    , "scatter3" ...
                                    , "lineScatter3" ...
                                    , "histogram" ...
                                    , "histogram2" ...
                                    , "histfit" ...
                                    , "contour" ...
                                    , "contourf" ...
                                    , "contour3" ...
                                    , "grid" ...
                                    ];
            end

            self.updateUser("adding the graphics tools...");
            prop="plot"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            self.plot = struct();
            self.resetPlot("hard");
            self.plot.reset = @self.resetPlot;
            self.plot.helpme = @self.helpme;

        end % constructor

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
                if strcmpi(varargin{1},"reset")
                    cmd = "doc self.resetPlot";
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

    end % methods (Access = public)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function reportWrongPlotName(self,plotNames)
            self.Err.marginTop = 1;
            self.Err.marginBot = 1;
            self.Err.msg    = "The input argument `plotNames` must be a string representing" + newline ...
                            + "the name of a plot belonging to the TabularFileContents class or," + newline ...
                            + "a list of such plot names. You have entered: " + plotNames + newline ...
                            + "Possible plots are: " + newline ...
                            + newline ...
                            + join(self.plotTypeList,newline) + newline ...
                            + newline ...
                            + "Here is the help for the ``reset()`` method: " + newline ...
                            + newline ...
                            + string(help("self.resetPlot")) ...
                            ;
            self.Err.abort();
        end

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
                is3d            = contains(requestedPlotTypeLower,"3") || isHistogram2;

                isDensityPlot = isHistfit || isHistogram2 || isHistogram || isContourf || isContour3 || isContour;

                if self.reportEnabled
                    self.Err.msg = msgPrefix + requestedPlotType + msgSuffix;
                    self.Err.note();
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%% reset line / scatter
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if isLine || isScatter

                    if resetTypeIsHard
                        self.plot.(requestedPlotType) = LineScatterPlot( requestedPlotType, self.df, @self.resetPlot );
                    else
                        self.plot.(requestedPlotType).resetInternal();
                    end

                    self.plot.(requestedPlotType).ccolumns = self.sampleLogFuncColName;
                    self.plot.(requestedPlotType).axes.kws.xscale = "linear";

                    if isScatter
                        if is3d
                            self.plot.(requestedPlotType).scatter.size = 10;
                        else
                            self.plot.(requestedPlotType).scatter.size = 7;
                        end
                    end

                    if isScatter && isLine

                        self.plot.(requestedPlotType).surface.enabled = false;
                        self.plot.(requestedPlotType).plot.enabled = true;
                        self.plot.(requestedPlotType).plot.kws.lineWidth = 1;
                        if is3d
                            self.plot.(requestedPlotType).plot.kws.color = [200 200 200 75] / 255;
                        else
                            self.plot.(requestedPlotType).plot.kws.color = uint8([200 200 200 200]);
                        end

                    elseif isLine

                        self.plot.(requestedPlotType).plot.enabled = false;
                        self.plot.(requestedPlotType).plot.kws.lineWidth = 1;
                        self.plot.(requestedPlotType).surface.enabled = true;
                        self.plot.(requestedPlotType).surface.kws.lineWidth = 1;

                    end

                    self.plot.(requestedPlotType).ycolumns = string(self.df.Properties.VariableNames{self.offset});%:end);
                    if is3d
                        if self.ndim==1
                            self.plot.(requestedPlotType).xcolumns = {};
                        else
                            self.plot.(requestedPlotType).xcolumns = string(self.df.Properties.VariableNames{self.offset});
                            self.plot.(requestedPlotType).ycolumns = string(self.df.Properties.VariableNames{self.offset+1});
                        end
                        self.plot.(requestedPlotType).zcolumns = self.sampleLogFuncColName;
                    end

                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%% reset histogram / histogram2 / histfit / contour / contourf / contour3
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if isDensityPlot

                    if resetTypeIsHard
                        self.plot.(requestedPlotType) = DensityPlot( requestedPlotType, self.df, @self.resetPlot );
                    else
                        self.plot.(requestedPlotType).resetInternal();
                    end

                    if isHistogram
                        self.plot.(requestedPlotType).xcolumns = string(self.df.Properties.VariableNames{self.offset}); %:self.offset+2);
                        self.plot.(requestedPlotType).histogram.kws.faceAlpha = 0.6;
                        self.plot.(requestedPlotType).histogram.kws.faceColor = "auto";
                        self.plot.(requestedPlotType).histogram.kws.edgeColor = "none";
                    else
                        self.plot.(requestedPlotType).xcolumns = string(self.df.Properties.VariableNames{self.offset});
                        if isHistogram2 || isContour || isContourf || isContour3
                            if self.ndim==1
                                self.plot.(requestedPlotType).ycolumns = string(self.df.Properties.VariableNames{self.offset-1});
                            else
                                self.plot.(requestedPlotType).ycolumns = string(self.df.Properties.VariableNames{self.offset+1});
                            end
                        end
                    end

                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%% reset grid
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if isGridPlot
                    %self.plot.(requestedPlotType).columns = string(self.df.Properties.VariableNames(self.offset:end));
                    if resetTypeIsHard
                        endIndx = min(lenVariableNames,self.offset+4);
                        self.plot.(requestedPlotType) = GridPlot( self.df, self.df.Properties.VariableNames(self.offset-1:endIndx), @self.resetPlot );
                    else
                        self.plot.(requestedPlotType).resetInternal();
                    end
                    self.plot.(requestedPlotType).ccolumn = string(self.df.Properties.VariableNames{self.offset-1});
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%% reset target component
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if (isLine || isScatter || isDensityPlot) && ~(is3d || self.isProgressFile)

                    xtarget = 0; % dummy
                    if isDensityPlot
                        xtarget = self.df.(self.plot.(requestedPlotType).xcolumns)(self.stats.maxLogFunc.idrow);
                    end

                    if self.plot.(requestedPlotType).type.is1d
                        self.plot.(requestedPlotType).target.values = [ xtarget, 0 ];
                    end

                    if self.plot.(requestedPlotType).type.is2d
                        ytarget = self.df.(self.plot.(requestedPlotType).ycolumns)(self.stats.maxLogFunc.idrow);
                        self.plot.(requestedPlotType).target.values = [ xtarget, ytarget ];
                    end

                    if isDensityPlot && self.plot.(requestedPlotType).type.is1d
                        self.plot.(requestedPlotType).target.hline.enabled = false;
                        self.plot.(requestedPlotType).target.scatter.enabled = false;
                    end
                    if isLine || isScatter
                        self.plot.(requestedPlotType).target.vline.enabled = true;
                    end
                    self.plot.(requestedPlotType).target.label = "maxLogFunc";

                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods (Access = Hidden)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % classdef TabularFileContents < handle
