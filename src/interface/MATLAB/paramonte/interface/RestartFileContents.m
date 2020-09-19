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
%   of ParaMonte output file contents. The ParaMonte read* methods
%   return an object or a list of objects of class RestartFileContents.
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
                self.resetPlot(plotType,"hard");
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

        function resetPlot(self,varargin)
            %
            %   Reset the properties of the plot to the original default settings.
            %   Use this method when you change many attributes of the plot
            %   and you want to clean up.
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
            %           This includes reading the original data frame again and resetting everything.
            %
            %   Returns
            %   -------
            %
            %       None. 
            %
            %   Example
            %   -------
            %
            %       reset("line3") % reset line3 plot to the default dettings
            %       reset("line3","hard") % regenerate line3 plot from scratch
            %       reset(["line","line3"],"hard") % regenerate line and line3 plots from scratch
            %       reset("hard") % regenrate all plots from scratch
            %
            resetTypeIsHard = false;
            requestedPlotTypeList = [];
            %plotTypeList =  [ "line" ...
            %                , "scatter" ...
            %                , "lineScatter" ...
            %                , "line3" ...
            %                , "scatter3" ...
            %                , "lineScatter3" ...
            %                , "histogram" ...
            %                , "histogram2" ...
            %                , "histfit" ...
            %                , "contour" ...
            %                , "contourf" ...
            %                , "contour3" ...
            %                , "grid" ...
            %                ];
            lenVariableNames = length(self.df.Properties.VariableNames);

            if nargin==1
                requestedPlotTypeList = self.plotTypeList;
            else
                argOne = string(varargin{1});
                if length(argOne)==1 && strcmpi(argOne,"hard")
                    requestedPlotTypeList = self.plotTypeList;
                    resetTypeIsHard = true;
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
                                    + "line, lineScatter." + newline ...
                                    );
                        end
                    end
                end
            end

            if nargin==3 && strcmpi(varargin{2},"hard"); resetTypeIsHard = true; end

            if resetTypeIsHard
                resetTypeIsHard = true;
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

                if self.reportEnabled
                    self.Err.msg = msgPrefix + requestedPlotType + msgSuffix;
                    self.Err.note();
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                isLineOrScatterPlot = contains(requestedPlotTypeLower,"scatter") || contains(requestedPlotTypeLower,"line");
                isEllipsoidPlot = contains(requestedPlotTypeLower,"covmat") || contains(requestedPlotTypeLower,"cormat");
                is3d = contains(requestedPlotTypeLower,"3") && ( isLineOrScatterPlot || isEllipsoidPlot );

                % line

                if strcmp(requestedPlotTypeLower,"line") || strcmp(requestedPlotTypeLower,"line3")
                    if resetTypeIsHard
                        self.plot.(requestedPlotType) = LineScatterPlot( requestedPlotType, self.df, @self.resetPlot );
                    else
                        self.plot.(requestedPlotType).resetInternal();
                    end
                    self.plot.(requestedPlotType).ycolumns = self.df.Properties.VariableNames(self.offset);%:end);
                    self.plot.(requestedPlotType).ccolumns = "sampleSize";
                    self.plot.(requestedPlotType).axes.kws.xscale = "linear";
                    self.plot.(requestedPlotType).plot.enabled = false;
                    self.plot.(requestedPlotType).plot.kws.lineWidth = 1;
                    self.plot.(requestedPlotType).surface.enabled = true;
                    self.plot.(requestedPlotType).surface.kws.lineWidth = 1;
                end

                % scatter / scatter3

                if strcmp(requestedPlotTypeLower,"scatter") || strcmp(requestedPlotTypeLower,"scatter3")
                    if resetTypeIsHard
                        self.plot.(requestedPlotType) = LineScatterPlot( requestedPlotType, self.df, @self.resetPlot );
                    else
                        self.plot.(requestedPlotType).resetInternal();
                    end
                    self.plot.(requestedPlotType).ccolumns = "sampleSize";
                    self.plot.(requestedPlotType).ycolumns = self.df.Properties.VariableNames(self.offset);%:end);
                    self.plot.(requestedPlotType).axes.kws.xscale = "linear";
                    self.plot.(requestedPlotType).scatter.size = 10;
                end

                % lineScatter / lineScatter3

                if strcmp(requestedPlotTypeLower,"linescatter") || strcmp(requestedPlotTypeLower,"linescatter3")
                    if resetTypeIsHard
                        self.plot.(requestedPlotType) = LineScatterPlot( requestedPlotType, self.df, @self.resetPlot );
                    else
                        self.plot.(requestedPlotType).resetInternal();
                    end
                    self.plot.(requestedPlotType).surface.enabled = false;
                    self.plot.(requestedPlotType).ccolumns = "sampleSize";
                    self.plot.(requestedPlotType).ycolumns = self.df.Properties.VariableNames(self.offset);%:end);
                    self.plot.(requestedPlotType).axes.kws.xscale = "linear";
                    if is3d
                        self.plot.(requestedPlotType).plot.kws.color = [200 200 200 75] / 255;
                    else
                        self.plot.(requestedPlotType).plot.kws.lineWidth = 1;
                        self.plot.(requestedPlotType).plot.kws.color = uint8([200 200 200 200]);
                        self.plot.(requestedPlotType).scatter.size = 20;
                    end
                end

                % 3d

                if isLineOrScatterPlot && is3d
                    if self.ndim==1
                        self.plot.(requestedPlotType).xcolumns = {};
                        self.plot.(requestedPlotType).ycolumns = self.df.Properties.VariableNames(self.offset);
                    else
                        self.plot.(requestedPlotType).xcolumns = self.df.Properties.VariableNames(self.offset);
                        self.plot.(requestedPlotType).ycolumns = self.df.Properties.VariableNames(self.offset+1);
                    end
                    self.plot.(requestedPlotType).zcolumns = "sampleSize";
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % hist / hist2 / histfit / contour / contourf / contour3

                isHist = strcmp(requestedPlotTypeLower,"histogram");
                isHist2 = strcmp(requestedPlotTypeLower,"histogram2");
                isHistfit = strcmp(requestedPlotTypeLower,"histfit");
                isContour = strcmp(requestedPlotTypeLower,"contour");
                isContourf = strcmp(requestedPlotTypeLower,"contourf");
                isContour3 = strcmp(requestedPlotTypeLower,"contour3");
                if isHist || isHist2 || isHistfit || isContour || isContourf || isContour3
                    if resetTypeIsHard
                        self.plot.(requestedPlotTypeLower) = HistPlot( self.df, requestedPlotTypeLower );
                    else
                        self.plot.(requestedPlotTypeLower).resetInternal();
                    end
                    if isHist
                        self.plot.(requestedPlotTypeLower).xcolumns = self.df.Properties.VariableNames(self.offset); %:self.offset+2);
                        self.plot.(requestedPlotTypeLower).histogram.kws.faceAlpha = 0.6;
                        self.plot.(requestedPlotTypeLower).histogram.kws.faceColor = "auto";
                        self.plot.(requestedPlotTypeLower).histogram.kws.edgeColor = "none";
                    else
                        self.plot.(requestedPlotTypeLower).xcolumns = self.df.Properties.VariableNames(self.offset);
                        if isHist2 || isContour || isContourf || isContour3
                            if self.ndim==1
                                self.plot.(requestedPlotTypeLower).ycolumns = self.df.Properties.VariableNames(self.offset-1);
                            else
                                self.plot.(requestedPlotTypeLower).ycolumns = self.df.Properties.VariableNames(self.offset+1);
                            end
                        end
                    end
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % grid

                isGrid = strcmp(requestedPlotTypeLower,"grid");
                if isGrid
                    %self.plot.(requestedPlotTypeLower).columns = string(self.df.Properties.VariableNames(self.offset:end));
                    if resetTypeIsHard
                        endIndex = min(lenVariableNames,self.offset+4);
                        self.plot.(requestedPlotTypeLower) = GridPlot( self.df, self.df.Properties.VariableNames(self.offset-1:endIndex));
                    else
                        self.plot.(requestedPlotTypeLower).resetInternal();
                    end
                    self.plot.(requestedPlotTypeLower).ccolumn = string(self.df.Properties.VariableNames(self.offset-1));
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % ellipsoid

                isCorMatPlot = strcmp(requestedPlotTypeLower,"cormat2") || strcmp(requestedPlotTypeLower,"cormat3");
                isCovMatPlot = strcmp(requestedPlotTypeLower,"covmat2") || strcmp(requestedPlotTypeLower,"covmat3");
                if isCorMatPlot || isCovMatPlot
                    plotType = "line"; if is3d; plotType = plotType + "3"; end
                    requestedPlotType = requestedPlotTypeLower;
                    if resetTypeIsHard
                        self.plot.(requestedPlotType) = EllipsoidPlot( plotType, self.df, @self.resetPlot );
                    else
                        self.plot.(requestedPlotType).resetInternal();
                    end
                    self.plot.(requestedPlotType).rows = self.plot.(requestedPlotType).getLogLinSpace ( 1.01 ... base
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
                    self.plot.(requestedPlotType).title.label = "Evolution of the " + matrixType + " matrices of the proposal distribution";
                    self.plot.(requestedPlotType).title.kws.fontsize = 11;
                    self.plot.(requestedPlotType).title.kws.interpreter = "tex";
                end

                % 3d

                if isEllipsoidPlot && is3d
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

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods (Access = Hidden)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % classdef RestartFileContents < handle

