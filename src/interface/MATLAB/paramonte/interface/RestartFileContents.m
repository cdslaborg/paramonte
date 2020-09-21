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
%   This is the RestartFileContents class for generating instances
%   of ParaMonte restart output file contents. The ParaDRAM readRestart() 
%   method returns an object or a list of objects of class RestartFileContents.
%
%   This is an internal ParaMonte class and it is 
%   not meant to be directly accessible to the users.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
classdef RestartFileContents < OutputFileContents

    properties(Access = public)
        count = [];
        ndim = [];
        df = [];
    end

    properties(Hidden)
        pdFieldNames =  [ "meanAcceptanceRateSinceStart" ...
                        , "sampleSize" ...
                        , "logSqrtDeterminant" ...
                        , "adaptiveScaleFactorSquared" ...
                        , "meanVec" ...
                        , "covMat" ...
                        , "corMat" ...
                        ];
        plotTypeList =  [ "line" ...
                        , "scatter" ...
                        , "lineScatter" ...
                        , "covmat2" ...
                        , "covmat3" ...
                        , "cormat2" ...
                        , "cormat3" ...
                        ..., "line3" ...
                        ..., "scatter3" ...
                        ..., "lineScatter3" ...
                        ..., "histogram" ...
                        ..., "histogram2" ...
                        ..., "histfit" ...
                        ..., "contour" ...
                        ..., "contourf" ...
                        ..., "contour3" ...
                        ..., "grid" ...
                        ];
        fileType = "restart";
        contents = [];
        lineList = [];
        lineListLen = [];
        newline = char(10);
        offset = 1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = RestartFileContents ( file ... % this must be a full file path
                                            , methodName ...
                                            , reportEnabled ...
                                            , Err ...
                                            )

            self = self@OutputFileContents(file,methodName,reportEnabled,Err);
            self.timer.tic();

            if strcmpi(self.methodName,"ParaDRAM") || strcmpi(self.methodName,"MatDRAM")
                self.readRestartParaDRAM()
            else
                error("Internal error occurred. Unrecognized methodName in the constructor of RestartFileContents: " + self.methodName)
            end

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

        function readRestartParaDRAM(self)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.contents = strrep(fileread(self.file),char(13),'');
            self.lineList = strsplit(self.contents,self.newline);
            self.lineListLen = length(self.lineList);

            % find count of updates

            prop="count"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            self.count = count(self.contents,self.pdFieldNames(1));

            % find ndim

            prop="ndim"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end

            rowOffset = 1;
            while ~contains(lower(self.lineList(rowOffset)),"meanvec")
                rowOffset = rowOffset + 1;
                if rowOffset>self.lineListLen; self.reportCorruptFile(); end
            end
            rowOffset = rowOffset + 1; % the first numeric value

            self.ndim = 0;
            while isNumericString(self.lineList(rowOffset+self.ndim))
                self.ndim = self.ndim + 1;
            end
            if self.ndim==0; self.reportCorruptFile(); end

            % parse the restart file contents

            proposalUpdates.(self.pdFieldNames(1)) = zeros(self.count,1);
            proposalUpdates.(self.pdFieldNames(2)) = zeros(self.count,1);
            proposalUpdates.(self.pdFieldNames(3)) = zeros(self.count,1);
            proposalUpdates.(self.pdFieldNames(4)) = zeros(self.count,1);
            proposalUpdates.(self.pdFieldNames(5)) = zeros(self.count,self.ndim);
            proposalUpdates.(self.pdFieldNames(6)) = zeros(self.count,self.ndim,self.ndim);
            proposalUpdates.(self.pdFieldNames(7)) = zeros(self.count,self.ndim,self.ndim);
            skip = 10 + self.ndim * (self.ndim + 3) / 2;
            for icount = 1:self.count
                if mod(icount,10)==0; self.updateProgess(icount/self.count); end
                istart = (icount-1) * skip + 1;
                rowOffset = 1; proposalUpdates.(self.pdFieldNames(1))(icount) = str2double( self.lineList(istart+rowOffset) );
                rowOffset = 3; proposalUpdates.(self.pdFieldNames(2))(icount) = str2double( self.lineList(istart+rowOffset) );
                rowOffset = 5; proposalUpdates.(self.pdFieldNames(3))(icount) = str2double( self.lineList(istart+rowOffset) );
                rowOffset = 7; proposalUpdates.(self.pdFieldNames(4))(icount) = str2double( self.lineList(istart+rowOffset) );
                rowOffset = 9;
                iend = istart + rowOffset + self.ndim;
                proposalUpdates.(self.pdFieldNames(5))(icount,1:self.ndim) = str2double( self.lineList(istart+rowOffset:iend-1) );
                for i = 1:self.ndim % covmat
                    istart = iend + 1;
                    iend = iend + i;
                    proposalUpdates.(self.pdFieldNames(6))(icount,1:i,i) = str2double( self.lineList(istart:iend) );
                    proposalUpdates.(self.pdFieldNames(6))(icount,i,1:i-1) = proposalUpdates.(self.pdFieldNames(6))(icount,1:i-1,i);
                end
                proposalUpdates.(self.pdFieldNames(7))(icount,:,:) = corrcov(squeeze(proposalUpdates.(self.pdFieldNames(6))(icount,:,:)));
            end
            self.updateProgess(1);
            self.df = table ( proposalUpdates.(self.pdFieldNames(1)) ...
                            , proposalUpdates.(self.pdFieldNames(2)) ...
                            , proposalUpdates.(self.pdFieldNames(3)) ...
                            , proposalUpdates.(self.pdFieldNames(4)) ...
                            , proposalUpdates.(self.pdFieldNames(5)) ...
                            , proposalUpdates.(self.pdFieldNames(6)) ...
                            , proposalUpdates.(self.pdFieldNames(7)) ...
                            );
            self.df.Properties.VariableNames = self.pdFieldNames;

            self.updateUser([]);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% graphics
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            prop="plot"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            self.plot = struct();

            for plotType =  self.plotTypeList
                self.resetPlot("hard",plotType);
                self.updateUser([]);
            end
            self.plot.reset = @self.resetPlot;
            self.plot.helpme = @self.helpme;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function reportCorruptFile(self)
            self.Err.marginTop = 1;
            self.Err.marginBot = 1;
            self.Err.msg    = "The structure of the file """ + self.file + """ does not match a " + self.methodName + " " + self.fileType + " file. " ...
                            + "The contents of the file may have been compromised. Verify the integrity of the contents of this file before attempting to reread it.";
            self.Err.abort();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function reportWrongPlotName(self,plotNames)
            self.Err.marginTop = 1;
            self.Err.marginBot = 1;
            self.Err.msg    = "The input argument `plotNames` must be a string representing" + newline ...
                            + "the name of a plot belonging to the RestartFileContents class or," + newline ...
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
                isCorMatPlot    = contains(requestedPlotTypeLower,"cormat");
                isCovMatPlot    = contains(requestedPlotTypeLower,"covmat");
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

                    self.plot.(requestedPlotType).ccolumns = "sampleSize";
                    self.plot.(requestedPlotType).ycolumns = string(self.df.Properties.VariableNames{self.offset});%:end);
                    self.plot.(requestedPlotType).axes.kws.xscale = "linear";

                    if isScatter
                        if is3d
                            self.plot.(requestedPlotType).scatter.size = 12;
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

                    if is3d
                        if self.ndim==1
                            self.plot.(requestedPlotType).xcolumns = {};
                            self.plot.(requestedPlotType).ycolumns = string(self.df.Properties.VariableNames{self.offset});
                        else
                            self.plot.(requestedPlotType).xcolumns = string(self.df.Properties.VariableNames{self.offset});
                            self.plot.(requestedPlotType).ycolumns = string(self.df.Properties.VariableNames{self.offset+1});
                        end
                        self.plot.(requestedPlotType).zcolumns = "sampleSize";
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
                %%%% reset ellipsoid
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if isCorMatPlot || isCovMatPlot
                    if resetTypeIsHard
                        self.plot.(requestedPlotType) = EllipsoidPlot( requestedPlotType, self.df, @self.resetPlot );
                    else
                        self.plot.(requestedPlotType).resetInternal();
                    end
                    self.plot.(requestedPlotType).rows = self.plot.(requestedPlotType).getLogLinSpace   ( 1.05 ... base
                                                                                                        , 1 ... logskip
                                                                                                        , 1 ... lowerLim
                                                                                                        , self.count ... upperLim
                                                                                                        );
                    if isCorMatPlot
                        self.plot.(requestedPlotType).matrixColumn = self.df.Properties.VariableNames(end);
                    elseif isCovMatPlot
                        self.plot.(requestedPlotType).matrixColumn = self.df.Properties.VariableNames(end-1);
                    end
                    self.plot.(requestedPlotType).centerColumn = self.df.Properties.VariableNames(end-2);
                    %self.plot.(requestedPlotType).ccolumn = "sampleSize";
                    self.plot.(requestedPlotType).axes.kws.xscale = "linear";
                    self.plot.(requestedPlotType).axes.kws.yscale = "linear";
                    self.plot.(requestedPlotType).plot.kws.lineWidth = 1;
                    matrixType = "covariance"; if isCorMatPlot; matrixType = "correlation"; end
                    self.plot.(requestedPlotType).title.enabled = true;
                    self.plot.(requestedPlotType).title.text = "Evolution of the " + matrixType + " matrices of the proposal distribution";
                    self.plot.(requestedPlotType).title.kws.fontSize = 11;
                    self.plot.(requestedPlotType).title.kws.interpreter = "tex";

                    % 3d

                    if is3d
                        %if self.ndim==1
                        %    self.plot.(requestedPlotType).matrixColumn = {};
                        %    self.plot.(requestedPlotType).centerColumn = self.df.Properties.VariableNames(self.offset);
                        %else
                        %    self.plot.(requestedPlotType).matrixColumn = self.df.Properties.VariableNames(self.offset);
                        %    self.plot.(requestedPlotType).centerColumn = self.df.Properties.VariableNames(self.offset+1);
                        %end
                        %self.plot.(requestedPlotType).zcolumn = "sampleSize";
                        self.plot.(requestedPlotType).axes.kws.zscale = "log";
                    end

                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods (Access = Hidden)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % classdef
