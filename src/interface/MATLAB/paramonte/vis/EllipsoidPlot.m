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
%   EllipsoidPlot(plotType, dataFrame, resetExternal)
%
%   This is the EllipsoidPlot class for generating instances of
%   ellipsoid-evolution plots via MATLAB's builtin function `plot()`.
%   It generates a plot of 2D ellipsoids corresponding to the input 
%   list of covariance/correlation matrices.
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
%           A MATLAB Table containing the covariance/correlation matrices 
%           that represent the characteristic covariance/correlation of ellipsoids. 
%           The covariance/correlation-matrix column of the input dataFrame must be a 3D 
%           matrix of the size (nrows,ndim,ndim) where count is the number of dataFrame nrows.
%
%   Attributes
%   ----------
%
%       dimensionPair
%
%           A pair of indices (vector of length 2) whose value determine the rows 
%           and columns from the covariance/correlation matrix which will be plotted.
%           The default value is dimensionPair = [1,2].
%
%           Example usage:
%
%               1.  dimensionPair = [1,2]
%               2.  dimensionPair = [3,1]
%
%           WARNING: In all cases, the indices must be distinct from each other and less 
%           WARNING: than ndim where ndim is the rank of the covariance/correlation matrix.
%
%       matrixColumn
%
%           The name of the column of the input dataFrame that represents the 
%           of the covariance/correlation or correlation matrices of the ellipsoids.
%
%       centerColumn
%
%           The name of the column of the input dataFrame that represents the 
%           corresponding center mean-vectors of the covariance/correlation matrices.
%
%       zcolumn (available only in 3D plotting objects)
%
%           Optional property that determines the column of dataFrame to serve as
%           the z-values. It can have multiple forms:
%
%               1.  a numeric or cell array of column indices in the input dataFrame.
%               2.  a string or cell array of column names in dataFrame.Properties.VariableNames.
%               3.  a cell array of a mix of the above two.
%               4.  a numeric range.
%
%           Example usage:
%
%               1.  zcolumn = 2
%               2.  zcolumn = "sampleSize"
%
%           WARNING: In all cases, zcolumn must have a length that is either 0, or 1
%           WARNING: If it is an empty object having length 0, then the default value will be used.
%
%           The default value is the index of the covariance/correlation matrices in the input data frame.
%
%       npoint
%
%           The number of points used to represent the ellipsoids. The higher the 
%           value of npoint is, the higher-resolution the ellipsoids would look like.
%
%           The default value is 100.
%
%       ccolumn (standing for color-columns)
%
%           Optional property that determines the columns of dataFrame to serve
%           as the color-mapping-values for each line/point element in the plot.
%           It can have multiple forms:
%
%               1.  a numeric or cell array of column indices in the input dataFrame.
%               2.  a string or cell array of column names in dataFrame.Properties.VariableNames.
%               3.  a cell array of a mix of the above two.
%               4.  a numeric range.
%
%           Example usage:
%
%               1.  ccolumn = 7
%               2.  ccolumn = "sampleSize"
%
%           WARNING: In all cases, ccolumn must have a length that is either 0, or 1.
%           WARNING: If the length is 1, then all data will be plotted with the same color 
%           WARNING: mapping determined by values specified by the elements of ccolumn. 
%           WARNING: If it is an empty object having length 0, then the default value will be used.
%
%           The default value is the indices of the rows of the input dataFrame.
%
%       colormap (available only in histogram2/contourf/contour3/contour objects) 
%
%           A MATLAB struct() with the following components:
%
%               enabled
%
%                   A logical value. If `true`, the colormap will be applied 
%                   to the plot. 
%
%               values
%
%                   A string or any other value that the colormap function 
%                   of MATLAB accepts as input.
%
%           Example usage:
%
%               1.  colormap.enabled = true;
%               1.  colormap.values = "winter";
%               1.  colormap.values = "autumn";
%
%       colorbar (available only in histogram2/contourf/contour3/contour objects) 
%
%           A MATLAB struct() with the following components:
%
%               enabled
%
%                   A logical value. If `true`, the colorbar will be applied to the plot.
%
%               kws
%
%                   A MATLAB struct() whose components' values are passed to 
%                   MATLAB's colorbar() function. If your desired attribute is 
%                   missing from the fieldnames of colorbar.kws, simply add 
%                   a new field named as the attribute and assign the desired 
%                   value to it.
%
%           Example usage:
%
%               colorbar.enabled = true; % add colorbar
%               colorbar.kws.location = "west";
%
%           WARNING
%
%               Keep in mind that MATLAB keyword arguments are case-INsensitive.
%               Hence, sure you do not add the keyword as multiple different fields.
%               For example, colorbar.kws.color and colorbar.kws.Color are the same,
%               and only one of the two will be processed.
%
%       plot
%
%           A MATLAB struct() with the following components:
%
%               enabled
%
%                   A logical value. If `true`, the plot() will be added to the plot.
%
%               kws
%
%                   A MATLAB struct() whose components' values are passed to 
%                   MATLAB's plot() function. If your desired attribute is 
%                   missing from the fieldnames of plot.kws, simply add 
%                   a new field named as the attribute and assign the desired 
%                   value to it.
%
%           Example usage:
%
%               plot.enabled = true; % add plot()
%               plot.kws.lineWidth = "none";
%
%           WARNING
%
%               Keep in mind that MATLAB keyword arguments are case-INsensitive.
%               Hence, sure you do not add the keyword as multiple different fields.
%               For example, plot.kws.lineWidth and plot.kws.LineWidth are the same,
%               and only one of the two will be processed.
%
%       scatter
%
%           A MATLAB struct() with the following components:
%
%               enabled
%
%                   A logical value. If `true`, the scatter() will be added to the plot.
%
%               kws
%
%                   A MATLAB struct() whose components' values are passed to 
%                   MATLAB's scatter() function. If your desired attribute is 
%                   missing from the fieldnames of scatter.kws, simply add 
%                   a new field named as the attribute and assign the desired 
%                   value to it.
%
%           Example usage:
%
%               scatter.enabled = true; % add scatter()
%
%           WARNING
%
%               Keep in mind that MATLAB keyword arguments are case-INsensitive.
%               Hence, sure you do not add the keyword as multiple different fields.
%               For example, scatter.kws.lineWidth and scatter.kws.LineWidth are the same,
%               and only one of the two will be processed.
%
%       target (available in 2D plots)
%
%           An object of class Target_class for adding target values to the plots.
%           For more information, see the documentation for the Target_class().
%
%   Superclass Attributes
%   ---------------------
%
%       See the documentation for the BasePlot class.
%
%   Returns
%   -------
%
%       An object of EllipsoidPlot class.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
classdef EllipsoidPlot < BasePlot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Access = public)
        dimensionPair
        matrixColumn
        centerColumn
        colorbar
        colormap
        ccolumn
        npoint
        title
    end

    properties (Hidden) % , Access = protected
        %dfref = [];
        %isdryrun = [];
        ndim
        is3d
        zdata
        isLinePlot = false;
        isScatterPlot = false;
        zcolnames
        zcolindex
        ccolnames
        ccolindex
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function resetInternal(self)

            resetInternal@BasePlot(self);
            self.dimensionPair = [];
            self.ccolumn = [];
            self.colormap = struct();
            self.colormap.enabled = true;
            self.colormap.values = [];
            self.npoint = [];

            self.title.kws = struct();
            self.title.enabled = false;
            self.title.text = [];
            self.title.kws.fontSize = 11;
            self.title.kws.interpreter = "tex";

            if self.type.is3d
                prop="zcolumn"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            end

            prop="plot"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            self.plot = struct();
            self.plot.enabled = true;
            self.plot.kws = struct();
            self.plot.kws.lineWidth = {};
            self.plot.kws.color = {};

            prop="scatter"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            self.scatter = struct();
            self.scatter.marker = [];
            self.scatter.filled = [];
            self.scatter.color = [];
            self.scatter.size = [];
            self.scatter.enabled = false;
            self.scatter.kws = struct();

            self.legend.enabled = false;

            self.colorbar = struct();
            self.colorbar.enabled = true;
            self.colorbar.kws = struct();
            self.colorbar.kws.fontSize = [];

            self.isdryrun = true;
            self.make();
            self.isdryrun = false;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = EllipsoidPlot(plotType, dataFrame, resetExternal)
            if nargin<3; resetExternal = []; end
            self = self@BasePlot(plotType, dataFrame, resetExternal);
            if nargin<3; self.resetExternal = @self.resetInternal; end
            self.resetInternal()
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

        function make(self,varargin)
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
            %       object's attributes, before calling the make() method.
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
            %       make( "gcf_kws", {"enabled",true,"color","none"} )
            %

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% parse arguments
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.parseArgs(varargin{:});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% set what to plot
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            dimPairLen = length(self.dimensionPair);
            if ~isnumeric(self.dimensionPair) || (dimPairLen~=0 && dimPairLen~=2)
                error   ( newline ...
                        + "dimensionPair must be a vector of length 2, representing the indices of the rows and the columns " ...
                        + "of the covariance/correlation matrix that will be used to form a 2-dimensional " ...
                        + "sub-covariance/correlation matrix. Example: dimensionPair = [1,2]" ...
                        + newline ...
                        );
            elseif dimPairLen==0
                self.dimensionPair = [1,2];
            end

            % activate at least one plot in the figure
            if ~self.plot.enabled
                warning ( newline ...
                        + "The line-plot type has been disabled by the user. There is nothing to display. " ...
                        + "To add at least one plot, set at least one the following components of the line-plot, " + newline ...
                        + newline ...
                        + "To add at least one plot, set at least one the following components of the line-plot, " + newline ...
                        + newline ...
                        + "    self.plot.enabled = true; % to generate single color monochromatic line plots" + newline ...
                        + newline ...
                        );
            end

            %%%% set scatter properties

            fname = "scatter";
            key = "marker"; val = "o"; if ~isfield(self.(fname),key) || isempty(self.(fname).(key)); self.(fname).(key) = val; end
            key = "filled"; val = true; if ~isfield(self.(fname),key) || isempty(self.(fname).(key)); self.(fname).(key) = val; end
            key = "size"; val = 5; if ~isfield(self.(fname),key) || isempty(self.(fname).(key)); self.(fname).(key) = val; end

            %%%% set line properties

            fname = "plot";
            key = "lineWidth"; val = 1; if isfield(self.(fname).kws,key) && isempty(self.(fname).kws.(key)); self.(fname).kws.(key) = val; end
            if self.scatter.enabled && self.plot.enabled
                key = "color"; val = uint8([200 200 200 150]); if ~isfield(self.(fname).kws,key) || isempty(self.(fname).kws.(key)); self.(fname).kws.(key) = val; end
            end

            if self.colormap.enabled && ~getVecLen(self.colormap.values)
                self.colormap.values = "winter";
            end

            %%%% set resolution

            if isempty(self.npoint)
                self.npoint = 100;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if self.isdryrun; return; end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% generate figure and axes if needed
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if self.figure.enabled
                figure_kws_cell = convertStruct2Cell(self.figure.kws,{"enabled","singleOptions"});
                % ATTN: DO NOT ERASE
                %if isfield(self.figure.kws,"singleOptions"); figure_kws_cell = { figure_kws_cell{:} , self.figure.kws.singleOptions{:} }; end
                self.currentFig.figure = figure( figure_kws_cell{:} );
            else
                set(0, "CurrentFigure", gcf);
                self.currentFig.figure = gcf;
                hold on;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% make 3d plot view, if requested
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if self.type.is3d; view(3); end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% check rows presence
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if getVecLen(self.rows)
                self.rowsindex = self.rows;
            else
                self.rowsindex = 1:1:length(self.dfref{:,1});
            end
            rowindexLen = length(self.rowsindex);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% check columns presence
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if getVecLen(self.matrixColumn)
                [covcolname, covcolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.matrixColumn); % m stands for (covariance/correlation) matrix.
                self.ndim = length(squeeze(self.dfref{1,covcolindex}));
            else
                error   ( newline ...
                        + "The column of the covariance/correlation matrices in the input dataFrame must be specified. " ...
                        + "No plots were made." ...
                        + newline ...
                        );
            end

            if getVecLen(self.centerColumn)
                [avgcolnames, avgcolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.centerColumn); % v stands for (the mean) vector.
            else
                avgcolnames = [];
                avgcolindex = [];
            end

            if self.type.is3d && getVecLen(self.zcolumn)
                [self.zcolnames, self.zcolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.zcolumn);
                self.zdata = self.dfref.(self.zcolnames)(self.rowsindex);
            else
                self.zcolindex = [];
                self.zcolnames = "Adaptive Update Count"; % "Count";
                self.zdata = self.rowsindex; % 1:1:rowindexLen;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% set color data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if self.colormap.enabled
                if getVecLen(self.ccolumn)
                    [self.ccolnames, self.ccolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.ccolumn);
                    cdata = self.dfref{self.rowsindex,self.ccolindex};
                elseif self.type.is3d
                    self.ccolindex = self.zcolindex;
                    self.ccolnames = self.zcolnames;
                    cdata = self.zdata;
                else
                    self.ccolindex = [];
                    self.ccolnames = "Adaptive Update Count"; % "Count";
                    cdata = self.rowsindex; % 1:1:rowindexLen;
                end
            else
                self.ccolindex = [];
                self.ccolnames = [];
                cdata = [];
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% check the lengths are consistent
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            zcolindexlen = length(self.zcolindex);
            ccolindexlen = length(self.ccolindex);
            covcolindexlen = length(covcolindex);
            avgcolindexlen = length(avgcolindex);
            maxcolindexlen = max(   [ covcolindexlen ...
                                    , avgcolindexlen ...
                                    , zcolindexlen ...
                                    , ccolindexlen ...
                                    ] ...
                                );

            if covcolindexlen~=maxcolindexlen && covcolindexlen>1; error("matrixColumn must be a unique name pointing to the column of the dataframe containing covariance/correlation matrix data."); end
            if avgcolindexlen~=maxcolindexlen && avgcolindexlen>1; error("centerColumn must be a unique name pointing to the column of the dataframe containing the centers of the covariance/correlation matrix data."); end
            if zcolindexlen~=maxcolindexlen && zcolindexlen>1; error("zcolumn must be a unique name pointing to the column of the dataframe that will be plotted on the z-axis."); end
            if ccolindexlen~=maxcolindexlen && ccolindexlen>1; error("ccolumn must be a unique name pointing to the column of the dataframe that will be used as the color map."); end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% get keyword arguments
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            excludes = {"enabled","singleOptions","color","cdata","size"};
            if self.colormap.enabled
                if isstring(self.colormap.values) || ischar(self.colormap.values)
                    try
                        cmap = eval(string(self.colormap.values)+"(rowindexLen)");
                    catch
                        error   ( "Failed to generate the color-mapping given the requested colormap: " + self.colormap.values ...
                                + "Please make sure the colormap string value is colormapping name recognized by the MATLAB's " ...
                                + "colormap() function." ...
                                );
                    end
                else
                    cmapsize = size(self.colormap.values);
                    if isnumeric(self.colormap.values) && length(cmapsize)==2 && cmapsize(1)==rowindexLen
                        cmap = self.colormap.values;
                    else
                        error   ( "A numeric value for the colormap must be given in the form of GRB triplet matrix, " ...
                                + "whose number of rows is the number of rows of the input dataframe to be visualized and, " ...
                                + "whose columns represent an RGB triplet." ...
                                );
                    end
                end
            end
            plot_kws_cell = convertStruct2Cell(self.plot.kws,excludes);
            scatter_kws_cell = convertStruct2Cell(self.scatter.kws,{"enabled","singleOptions","color","cdata","size"});
            if self.scatter.filled; scatter_kws_cell = {"filled", scatter_kws_cell{:}}; end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% add plot
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            lglabels = [];

            if ~self.type.is3d
                self.target.currentFig = struct();
                self.target.counter.xline = 0;
                self.target.counter.yline = 0;
                self.target.counter.scatter = 0;
            end

            if self.type.is3d
                zdataOnes = ones(self.npoint,1);
            end

            if getVecLen(self.plot.kws.color)
                lineColorKeyVal = {"color",self.plot.kws.color};
            else
                lineColorKeyVal = {};
            end

            if getVecLen(self.scatter.color)
                scatterColor = {self.scatter.color};
            else
                scatterColor = {};
            end

            counter.plot = 0;
            counter.plot3 = 0;
            counter.surface = 0;
            counter.scatter = 0;
            counter.scatter3 = 0;
            if self.plot.enabled && self.type.is3d
                self.currentFig.plot3 = cell(rowindexLen,1);
            else
                self.currentFig.plot = cell(rowindexLen,1);
            end
            if self.scatter.enabled && self.type.is3d
                self.currentFig.scatter3 = cell(rowindexLen,1);
            else
                self.currentFig.scatter = cell(rowindexLen,1);
            end

            meanVec = zeros(dimPairLen,1);
            for irow = 1:rowindexLen

                %%%% set colormapping if enabled

                if self.colormap.enabled
                    lineColorKeyVal = {"color",cmap(irow,:)};
                    scatterColor = {cmap(irow,:)};
                end

                %%%% compute the ellipsoid characteristics

                submatrix = squeeze(self.dfref{self.rowsindex(irow),covcolindex});
                submatrix = submatrix(self.dimensionPair,self.dimensionPair);
                if avgcolindexlen>0
                    meanVec = squeeze(self.dfref{self.rowsindex(irow),avgcolindex});
                    meanVec = meanVec(self.dimensionPair);
                end

                % get ellipsoid boundary

                bcrd = self.getEllipsoidBoundary( submatrix ...
                                                , meanVec ...
                                                , self.npoint ...
                                                );


                %%%% add plot

                if self.plot.enabled
                    if self.type.is3d
                        counter.plot3 = counter.plot3 + 1;
                        self.currentFig.plot3{counter.plot3}    = plot3 ( bcrd(:,1) ...
                                                                        , bcrd(:,2) ...
                                                                        , self.zdata(irow)*zdataOnes ...
                                                                        , plot_kws_cell{:} ...
                                                                        , lineColorKeyVal{:} ...
                                                                        );
                    else
                        counter.plot = counter.plot + 1;
                        self.currentFig.plot{counter.plot}  = plot  ( bcrd(:,1) ...
                                                            , bcrd(:,2) ...
                                                            , plot_kws_cell{:} ...
                                                            , lineColorKeyVal{:} ...
                                                            );
                    end
                    hold on;
                end

                if self.scatter.enabled
                    if self.type.is3d
                        counter.scatter3 = counter.scatter3 + 1;
                        self.currentFig.scatter3{counter.scatter3} = scatter3   ( bcrd(:,1) ...
                                                                                , bcrd(:,2) ...
                                                                                , self.zdata(irow)*zdataOnes ...
                                                                                , self.scatter.size ...
                                                                                , scatterColor{:} ...
                                                                                , self.scatter.marker ...
                                                                                , scatter_kws_cell{:} ...
                                                                                );
                    else
                        counter.scatter = counter.scatter + 1;
                        self.currentFig.scatter{counter.scatter} = scatter  ( bcrd(:,1) ...
                                                                            , bcrd(:,2) ...
                                                                            , self.scatter.size ...
                                                                            , scatterColor{:} ...
                                                                            , self.scatter.marker ...
                                                                            , scatter_kws_cell{:} ...
                                                                            );
                    end
                    hold on;
                end

            end % loop plot

            self.currentFig.axes = gca;

            % set gca properties
            %gca_kws_cell = convertStruct2Cell(self.axes.kws,{"enabled","singleOptions"});
            %%if isfield(self.axes.kws,"singleOptions"); gca_kws_cell = { gca_kws_cell{:} , self.axes.kws.singleOptions{:} }; end
            %self.currentFig.axes = set( self.currentFig.axes, gca_kws_cell{:} );

            if self.colormap.enabled
                colormap(self.currentFig.axes,self.colormap.values);
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% add axis labels
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.currentFig.xlabel = xlabel("Dimension " + string(self.dimensionPair(1)), "Interpreter", "none");
            self.currentFig.ylabel = ylabel("Dimension " + string(self.dimensionPair(2)), "Interpreter", "none");

            if self.type.is3d
                self.currentFig.zlabel = zlabel(self.zcolnames(1), "Interpreter", "none");
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% set title properties
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % this is now done in the BasePlot
            %if self.title.enabled && getVecLen(self.title.text)
            %    if ~(isstring(self.title.text) || ischar(self.title.text))
            %        error("The title component of an EllipsoidPlot object must be a string or character vector.");
            %    end
            %    title_kws_cell = convertStruct2Cell(self.title.kws,{"enabled","text"});
            %    title(self.title.text,title_kws_cell{:});
            %end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% add line colorbar
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if self.colorbar.enabled && self.colormap.enabled
                if isempty(self.colorbar.kws.fontSize) || ~isa(self.colorbar.kws.fontSize,"numeric")
                    self.colorbar.kws.fontSize = self.currentFig.ylabel.FontSize;
                end
                caxis(self.currentFig.axes, [min(cdata), max(cdata)] );
                colorbar_kws_cell = convertStruct2Cell(self.colorbar.kws,{"enabled","singleOptions"});
                self.currentFig.colorbar = colorbar(colorbar_kws_cell{:});
                ylabel(self.currentFig.colorbar,self.ccolnames(1),"fontSize",self.colorbar.kws.fontSize, "interpreter", "none");
            else
                colorbar('off');
                self.currentFig.colorbar = [];
            end

            if self.legend.enabled && (~isfield(self.legend,"labels") || isempty(self.legend.labels)); self.legend.labels = lglabels; end
            self.doBasePlotStuff();

            hold off;
            %box on; grid on; hold off;

        end % function plot

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Static)

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function bcrd = getEllipsoidBoundary(covMat, meanVec, npoint) % returns the coordinates of the boundary of the ellipsoid
            if isempty(npoint); npoint = 50; end
            independentVariable = linspace(0, 2*pi, npoint)';
            xval = cos(independentVariable);
            yval = sin(independentVariable);
            ap = [xval(:) yval(:)]';
            [eigenVectors,eigenValues] = eig(covMat);
            eigenValues = sqrt(eigenValues); % convert variance to std
            bcrd = transpose( eigenVectors * eigenValues * ap + repmat(meanVec(:), 1, npoint) );
            %h = plot(bcrd(:,1), bcrd(:,2), '-');
        end

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % classdef EllipsoidPlot
