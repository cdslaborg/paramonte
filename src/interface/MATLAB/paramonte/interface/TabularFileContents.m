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
%   This is the TabularFileContents class for generating instances 
%   of ParaMonte output file contents. The ParaMonte read* methods
%   return an object or a list of objects of class TabularFileContents.
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
        offset = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = TabularFileContents ( file ...
                                            , fileType ...
                                            , delimiter ...
                                            , methodName ...
                                            , mpiEnabled ...
                                            , markovChainRequested ...
                                            , Err ...
                                            )

            self = self@OutputFileContents(file,methodName,mpiEnabled,Err);
            self.delimiter = delimiter;
            self.timer.tic();

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            d = importdata  ( self.file ...
                            ..., "delimiter", self.delimiter ...
                            );
            if isfield(d,"colheaders")
                colheadersLen = length(d.colheaders);
            else
                self.Err.marginTop = 1;
                self.Err.marginBot = 1;
                self.Err.msg    = "The structure of the file """ + self.file + """ does not match a " + methodName + " " + fileType + " file. " ...
                                + "Verify the contents of this file before attempting to read this file.";
                self.Err.abort();
            end
            for icol = 1:colheadersLen
                if strcmp(d.colheaders{icol},"SampleLogFunc")
                    break
                end
            end
            self.offset = icol + 1; % index of the first variable
            prop="ndim"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            self.ndim   = colheadersLen - self.offset + 1;
            prop="count"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            self.count  = length(d.data(:,1));

            if markovChainRequested
                cumSumWeight = cumsum(d.data(:,self.offset-2));
                if cumSumWeight(end) > self.count % it is indeed a compact chain
                    dMarkov = zeros( cumSumWeight(end) , self.ndim + self.offset - 1 );
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

            if ~self.mpiEnabled
                self.Err.marginTop = 0;
                self.Err.marginBot = 1;
                self.Err.msg = "ndim = " + string(self.ndim) + ", count = " + string(self.count);
                self.Err.note();
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% statistics
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            prop="stats"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end

            self.stats = struct();

            % add chain cormat

            self.updateUser("computing the sample correlation matrix...");
            self.stats.cormat = CorCovMat   ( self.df ...
                                            , self.offset:self.offset+self.ndim-1 ...
                                            , "pearson" ... method
                                            , [] ... rows
                                            , self.Err ...
                                            );
            self.updateUser([]);

            % add chain covmat

            self.updateUser("computing the sample covariance matrix...");
            self.stats.covmat = CorCovMat   ( self.df ...
                                            , self.offset:self.offset+self.ndim-1 ...
                                            , [] ... method
                                            , [] ... rows
                                            , self.Err ...
                                            );
            self.updateUser([]);

            % add chain autocorr

            self.updateUser("computing the sample autocorrelation...");
            self.stats.autocorr = AutoCorr_class( self.df ...
                                                , self.offset-1:self.offset+self.ndim-1 ...
                                                , [] ... rows
                                                , self.Err ...
                                                );
            self.updateUser([]);

            %self.stats.maxLogFunc = getMax(self.df,"SampleLogFunc");

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% graphics
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            prop="plot"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            self.plot = struct();

            for plotType = ["line","line3","scatter","scatter3","lineScatter","lineScatter3","histogram","histogram2","histfit","contour","contourf","contour3","grid"]
                %self.updateUser("generating " + plotType + " plot...");
                self.resetPlot(plotType,"hard");
                self.updateUser([]);
            end
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

        function resetPlot(self,varargin)
            %
            %   Reset the properties of the plot to the original default settings.
            %   Use this method when you change many attributes of the plot and 
            %   you want to clean up and go back to the default settings.
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
            %       reset("line3") % reset line3 plot to the default settings
            %       reset("line3","hard") % regenerate line3 plot from scratch
            %       reset(["line","line3"],"hard") % regenerate line and line3 plots from scratch
            %       reset("hard") % regenerate all plots from scratch
            %
            resetTypeIsHard = false;
            requestedPlotTypeList = [];
            plotTypeList = ["line","scatter","lineScatter","line3","scatter3","lineScatter3","histogram","histogram2","histfit","contour","contourf","contour3","grid"];
            lenVariableNames = length(self.df.Properties.VariableNames);

            if nargin==1
                requestedPlotTypeList = plotTypeList;
            else
                argOne = string(varargin{1});
                if length(argOne)==1 && strcmpi(argOne,"hard")
                    requestedPlotTypeList = plotTypeList;
                    resetTypeIsHard = true;
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

                self.Err.msg = msgPrefix + requestedPlotType + msgSuffix;
                self.Err.note();

                plotName = "";

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                is3d = false;
                if contains(requestedPlotTypeLower,"3") && ( contains(requestedPlotTypeLower,"scatter") || contains(requestedPlotTypeLower,"line"))
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
                    self.plot.(plotName).ycolumns = self.df.Properties.VariableNames(self.offset);%:end);
                    self.plot.(plotName).ccolumns = "SampleLogFunc";
                    self.plot.(plotName).gca_kws.xscale = "linear";
                    self.plot.(plotName).plot_kws.enabled = false;
                    self.plot.(plotName).plot_kws.linewidth = 1;
                    self.plot.(plotName).surface_kws.enabled = true;
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
                    self.plot.(plotName).ccolumns = "SampleLogFunc";
                    self.plot.(plotName).ycolumns = self.df.Properties.VariableNames(self.offset);%:end);
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
                    self.plot.(plotName).ccolumns = "SampleLogFunc";
                    self.plot.(plotName).ycolumns = self.df.Properties.VariableNames(self.offset);%:end);
                    self.plot.(plotName).gca_kws.xscale = "linear";
                    if is3d
                        self.plot.(plotName).plot_kws.color = [200 200 200 75] / 255;
                    else
                        self.plot.(plotName).plot_kws.linewidth = 1;
                        self.plot.(plotName).plot_kws.color = uint8([200 200 200 200]);
                        self.plot.(plotName).scatter_kws.size = 20;
                    end
                end

                % 3d

                if is3d
                    if self.ndim==1
                        self.plot.(plotName).xcolumns = {};
                        self.plot.(plotName).ycolumns = self.df.Properties.VariableNames(self.offset);
                    else
                        self.plot.(plotName).xcolumns = self.df.Properties.VariableNames(self.offset);
                        self.plot.(plotName).ycolumns = self.df.Properties.VariableNames(self.offset+1);
                    end
                    self.plot.(plotName).zcolumns = "SampleLogFunc";
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
                        self.plot.(requestedPlotTypeLower).reset();
                    end
                    if isHist
                        self.plot.(requestedPlotTypeLower).xcolumns = self.df.Properties.VariableNames(self.offset); %:self.offset+2);
                        self.plot.(requestedPlotTypeLower).histogram_kws.facealpha = 0.6;
                        self.plot.(requestedPlotTypeLower).histogram_kws.facecolor = "auto";
                        self.plot.(requestedPlotTypeLower).histogram_kws.edgecolor = "none";
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
                        endIndx = min(lenVariableNames,self.offset+4);
                        self.plot.(requestedPlotTypeLower) = GridPlot( self.df, self.df.Properties.VariableNames(self.offset-1:endIndx));
                    else
                        self.plot.(requestedPlotTypeLower).reset();
                    end
                    self.plot.(requestedPlotTypeLower).ccolumn = string(self.df.Properties.VariableNames(self.offset-1));
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods (Access = Hidden)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % classdef TabularFileContents < handle