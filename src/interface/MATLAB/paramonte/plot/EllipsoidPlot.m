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
%   EllipsoidPlot(dataFrame)
%
%   This is the EllipsoidPlot class for generating instances of
%   ellipsoid-evolution plots via MATLAB's builtin function `plot()`.
%   It generates a plot of 2D ellipsoids corresponding to the input 
%   list of covariance (or correlation) matrices.
%
%   NOTE: This is a low-level ParaMonte class and is not meant
%   NOTE: to be directly instantiated by the user.
%
%   Parameters
%   ----------
%
%       dataFrame
%
%           A MATLAB Table containing the covariance (or the correlation) 
%           matrices that represent the characteristic covariance of ellipsoids.
%           The covariance-matrix column of the input dataFrame must be a 3D matrix 
%           of the size (nrows,ndim,ndim) where count is the number of dataFrame nrows.
%
%           NOTE: This is a low-level internal argument and is not meant
%           NOTE: to be accessed or be provided by the user.
%
%   Attributes
%   ----------
%
%       dimensionPair
%
%           A pair of indices (vector of length 2) whose value determine 
%           the rows and columns from the covariance matrix which will be plotted.
%
%           Example usage:
%
%               1.  dimensionPair = [1,2]
%               2.  dimensionPair = [3,1]
%
%           WARNING: In all cases, the indices must distinct from each other and 
%           WARNING: <ndim where ndim is the rank of the covariance matrix.
%           The default value is dimensionPair = [1,2].
%
%       covMatColumn
%
%           The name of the column of the input dataFrame that represents the 
%           of the covariance matrices.
%
%       centerColumn
%
%           The name of the column of the input dataFrame that represents the 
%           corresponding center mean-vectors of the covariance matrices.
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
%           The default value is the index of the covariance matrices in the input data frame.
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
%           optional property that determines the columns of dataFrame to serve
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
%               1.  ccolumn = [7,8,9]
%               2.  ccolumn = ["SampleLogFunc","SampleVariable1"]
%               3.  ccolumn = {"SampleLogFunc",9,"SampleVariable1"}
%               4.  ccolumn = 7:9      # every column in the data frame starting from column #7 to #9
%               5.  ccolumn = 7:2:20   # every other column in the data frame starting from column #7 to #20
%
%           WARNING: In all cases, ccolumn must have a length that is either 0, or 1, or equal
%           WARNING: to the maximum lengths of (covMatColumns,zcolumn). If the length is 1, then all data
%           WARNING: will be plotted with the same color mapping determined by values specified by the elements
%           WARNING: of ccolumn. If it is an empty object having length 0, then the default value will be used.
%
%           The default value is the indices of the rows of the input dataFrame.
%
%       colormap
%
%           A MATLAB struct() property with two components:
%
%               1. enabled: logical value. If `true`, the colormap will be applied to the plot
%               1. values: a string or any other value that the colormap function of MATLAB accepts as input.
%
%           Example usage:
%
%               1.  colormap.values = "autumn"
%               1.  colormap.values = "winter"
%
%           If colormap is not provided or is empty, the default will be "winter".
%
%       colorbar_kws
%
%           A MATLAB struct() whose components' values are passed to MATLAB's colorbar function.
%           If your desired attribute is missing from the fieldnames of colorbar_kws, simply add
%           a new field named as the attribute and assign the desired value to it.
%
%           Example usage:
%
%               colorbar_kws.enabled = true % add colorbar
%               colorbar_kws.location = "west"
%
%           If a desired property is missing among the struct fields, simply add the field
%           and its value to colorbar_kws.
%
%           WARNING: keep in mind that MATLAB keyword arguments are case-INsensitive.
%           WARNING: therefore make sure you do not add the keyword as multiple different fields.
%           WARNING: For example, colorbar_kws.color and colorbar_kws.Color are the same,
%           WARNING: and only one of the two will be processed.
%
%       plot_kws (available only in line/line3 objects)
%
%           A MATLAB struct() whose fields (with the exception of few, e.g., enabled, singleOptions, ...)
%           are directly passed to the `plot()` function of MATLAB (or plot3() if plot is 3d).
%           This property exists only if the object is instantiated as a line/line3 object.
%
%           Example usage:
%
%               plot_kws.enabled = true; % add plot()
%               plot_kws.linewidth = 2;
%
%           If a desired property is missing among the struct fields, simply add the field
%           and its value to plot_kws.
%
%           WARNING: keep in mind that MATLAB keyword arguments are case-INsensitive.
%           WARNING: therefore make sure you do not add the keyword as multiple different fields.
%           WARNING: For example, plot_kws.color and plot_kws.Color are the same,
%           WARNING: and only one of the two will be processed.
%
%       target
%
%           an object of class Target_class for adding target values to the plots.
%           For more information, see the documentation for Target.
%
%   Superclass Attributes
%   ----------------------
%
%       See the documentation for the BasePlot class
%
%   Returns
%   -------
%
%       an object of EllipsoidPlot class
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
classdef EllipsoidPlot < BasePlot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Access = public)

        dimensionPair
        covMatColumn
        centerColumn
        ccolumn
        colorbar_kws
        colormap
        target
        npoint

    end

    properties (Access = protected, Hidden)
        %dfref = [];
        %isdryrun = [];
        ndim
        is3d
        plotType
        isLinePlot = false;
        isScatterPlot = false;
        surface_kws = []; % dummy variable
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function reset(self)

            reset@BasePlot(self);
            self.dimensionPair = {};
            self.ccolumn = {};
            self.colormap = struct();
            self.colormap.enabled = true;
            self.colormap.values = [];
            self.npoint = [];

            self.isLinePlot = false;
            if contains(self.plotType,"line")
                self.isLinePlot = true;
                prop="plot_kws"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            end
            if contains(self.plotType,"3")
                prop="zcolumn"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                self.is3d = true;
            else
                self.is3d = false;
            end

            if self.isLinePlot

                self.plot_kws = struct();
                self.plot_kws.enabled = true;
                self.plot_kws.linewidth = {};
                self.plot_kws.singleOptions = {};
                self.plot_kws.color = {};

            end

            self.legend_kws.enabled = false;

            self.colorbar_kws = struct();
            self.colorbar_kws.enabled = true;
            self.colorbar_kws.fontsize = [];

            self.isdryrun = true;
            self.plot();
            self.isdryrun = false;

            self.target = Target_class();

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = EllipsoidPlot(varargin) % expected input arguments: dataFrame, plotType

            self = self@BasePlot(varargin{1});
            self.plotType = lower(string(varargin{2}));
            self.reset()
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

        function plot(self,varargin)
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
            %       plot( "gcf_kws", {"enabled",true,"color","none"} )
            %

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% parse arguments
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.parseArgs(varargin{:});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% set what to plot
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if self.isLinePlot
                % activate at least one plot in the figure
                if ~self.plot_kws.enabled
                    warning ( newline ...
                            + "The line plot() type has been disabled by the user. There is nothing to display. " ...
                            + "To add at least one plot, set at least one the following components of the line-plot, " + newline  + newline  ...
                            + "To add at least one plot, set at least one the following components of the line-plot, " + newline  + newline  ...
                            + "    self.plot_kws.enabled = true; % to generate single color monochromatic line plots" + newline ...
                            + "You can also pass these arguments at the time of calling the plot() method:" + newline  + newline  ...
                            + "    self.plot(""plot_kws"",{""enabled"",true}); % to generate single color monochromatic line plots" + newline ...
                            + newline ...
                            );
                end
            end

            if self.isLinePlot && self.plot_kws.enabled
                key = "linewidth"; val = 1;
                if isfield(self.plot_kws,key) && isempty(self.plot_kws.(key))
                    self.plot_kws.(key) = val;
                end
            end

            if self.colormap.enabled && ~getVecLen(self.colormap.values)
                self.colormap.values = "winter";
                %if self.is3d
                %    self.colormap = "winter";
                %else
                %    self.colormap = "autumn";
                %end
            end

            if isempty(self.npoint)
                self.npoint = 100;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if self.isdryrun; return; end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % generate figure and axes if needed

            if self.gcf_kws.enabled
                gcf_kws_cell = convertStruct2Cell(self.gcf_kws,{"enabled","singleOptions"});
                if isfield(self.gcf_kws,"singleOptions"); gcf_kws_cell = { gcf_kws_cell{:} , self.gcf_kws.singleOptions{:} }; end
                self.currentFig.gcf = figure( gcf_kws_cell{:} );
            else
                set(0, "CurrentFigure", gcf);
                self.currentFig.gcf = gcf;
                hold on;
            end

            % check rows presence

            if getVecLen(self.rows)
                rowindex = self.rows;
            else
                rowindex = 1:1:length(self.dfref{:,1});
            end
            rowindexLen = length(rowindex);

            % check columns presence

            if getVecLen(self.covMatColumn)
                [mcolnames, mcolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.covMatColumn); % m stands for (covariance) matrix.
                self.ndim = length(squeeze(self.dfref{1,mcolindex}));
            else
                error   ( newline ...
                        + "The column of the covariance matrices in the input dataFrame must be specified. " ...
                        + "No plots were made." ...
                        + newline ...
                        );
            end

            if getVecLen(self.centerColumn)
                [vcolnames, vcolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.centerColumn); % v stands for (the mean) vector.
            else
                vcolnames = [];
                vcolindex = [];
            end

            % set color data

            if self.colormap.enabled
                if getVecLen(self.ccolumn)
                    [ccolnames, ccolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.ccolumn);
                else
                    ccolindex = [];
                    ccolnames = "Count";
                    cdata = 1:1:rowindexLen;
                end
            else
                ccolindex = [];
                ccolnames = [];
                cdata = [];
            end

            if self.is3d && getVecLen(self.zcolumn)
                [zcolnames, zcolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.zcolumn);
            else
                zcolindex = [];
                zcolnames = "Count";
                zdata = rowindex;
            end

            % check the lengths are consistent

            mcolindexlen = length(mcolindex);
            vcolindexlen = length(vcolindex);
            zcolindexlen = length(zcolindex);
            ccolindexlen = length(ccolindex);
            maxLenColumns = max (   [ mcolindexlen ...
                                    , vcolindexlen ...
                                    , zcolindexlen ...
                                    , ccolindexlen ...
                                    ] ...
                                );

            if mcolindexlen~=maxLenColumns && mcolindexlen>1; error("length of mcolumns must be either 1 or equal to the maximum of the lengths of vcolumns, zcolumn, ccolumn."); end
            if vcolindexlen~=maxLenColumns && vcolindexlen>1; error("length of vcolumns must be either 1 or equal to the maximum of the lengths of mcolumns, zcolumn, ccolumn."); end
            if zcolindexlen~=maxLenColumns && zcolindexlen>1; error("length of zcolumn must be either 1 or equal to the maximum of the lengths of mcolumns, vcolumns, ccolumn."); end
            if ccolindexlen~=maxLenColumns && ccolindexlen>1; error("length of ccolumn must be either 1 or equal to the maximum of the lengths of mcolumns, vcolumns, zcolumn."); end

            % assign data in case of single column assignments

            %mdata = [];
            %if mcolindexlen==1
            %    mdata = self.dfref{rowindex,mcolindex}(:,:);
            %end
            %if vcolindexlen==1
            %    vdata = self.dfref{rowindex,vcolindex};
            %end
            %if zcolindexlen==1
            %    zdata = self.dfref{rowindex,zcolindex};
            %end
            %if ccolindexlen==1
            %    cdata = self.dfref{rowindex,ccolindex};
            %end

            %%%%%%%%%%%%%%%%%%%%%%%
            % get keyword arguments
            %%%%%%%%%%%%%%%%%%%%%%%

            if self.isLinePlot
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
                    if numeric(self.colormap.values) && length(cmapsize)==2 && cmapsize(1)==rowindexLen
                        cmap = self.colormap.values;
                    else
                        error   ( "A numeric value for the colormap must be given in the form of GRB triplet matrix, " ...
                                + "whose number of rows is the number of rows of the input dataframe to be visualized and, " ...
                                + "whose columns represent an RGB triplet." ...
                                );
                    end
                end
                excludes = {"enabled","singleOptions"};
                if self.colormap.enabled; excludes = {excludes{:},"color"}; end
                plot_kws_cell = convertStruct2Cell(self.plot_kws,excludes);
            end
            %if self.isScatterPlot
            %    scatter_kws_cell = convertStruct2Cell(self.scatter_kws,{"enabled","singleOptions","cdata","size"});
            %end
            %if self.isLinePlot
            %    surface_kws_cell = convertStruct2Cell(self.surface_kws,{"enabled","singleOptions","color","size"});
            %end

            %%%%%%%%%%%%%%%%%%%%%%%
            % add line/scatter plot
            %%%%%%%%%%%%%%%%%%%%%%%

            lglabels = [];

            if self.isLinePlot && self.plot_kws.enabled

                if ~self.is3d
                    zdata = zeros(rowindexLen,1);
                end

                if ~getVecLen(self.plot_kws.color)
                    colorKeyVal = {"color",self.plot_kws.color};
                end

                meanVec = zeros(self.ndim,1);
                for irow = 1:rowindexLen

                    %if mcolindexlen>1
                    %    mdata = self.dfref{rowindex,mcolindex(irow)};
                    %end
                    %if vcolindexlen>1
                    %    vdata = self.dfref{rowindex,vcolindex(irow)};
                    %end
                    %if zcolindexlen>1
                    %    zdata = self.dfref{rowindex,zcolindex(irow)};
                    %end
                    %if ccolindexlen>1
                    %    cdata = self.dfref{rowindex,ccolindex(irow)};
                    %end
                    %if self.legend_kws.enabled && ~self.is3d
                    %    if mcolindexlen<2 && vcolindexlen>=1
                    %        lglabels = [ lglabels , vcolnames(irow) ];
                    %    elseif mcolindexlen>1 && vcolindexlen<2
                    %        lglabels = [ lglabels , mcolnames(irow) ];
                    %    else
                    %        lglabels = [ lglabels , mcolnames(irow)+"-"+vcolnames(irow) ];
                    %    end
                    %end

                    %if isMultiColorScatterPlot
                    %    currentScatterMarkerColor = scatterMultiColorList(irow,:);
                    %end

                    % add plot

                    if self.colormap.enabled
                        colorKeyVal = {"color",cmap(irow,:)};
                    end
                    covMat = squeeze(self.dfref{rowindex(irow),mcolindex});
                    if vcolindexlen>0; meanVec = squeeze(self.dfref{rowindex(irow),vcolindex}); end
                    bcrd = self.makeEllipsoid   ( covMat ...
                                                , meanVec ...
                                                , self.npoint ...
                                                );

                    if self.isLinePlot
                        if self.plot_kws.enabled
                            if self.is3d
                                error("plot3 not implemented")
                                %self.currentFig.plot3 = plot3   ( mdata ...
                                %                                , vdata ...
                                %                                , zdata ...
                                %                                , plot_kws_cell{:} ...
                                %                                );
                            else
                                self.currentFig.plot = plot ( bcrd(:,1) ...
                                                            , bcrd(:,2) ...
                                                            , plot_kws_cell{:} ...
                                                            , colorKeyVal{:} ...
                                                            );
                            end
                            hold on;
                        end
                        %if self.surface_kws.enabled && getVecLen(self.colormap.values)
                        %    self.currentFig.surface = surface   ( "XData",[mdata(:) mdata(:)] ...
                        %                                        , "YData",[vdata(:) vdata(:)] ...
                        %                                        , "ZData",[zdata(:) zdata(:)] ...
                        %                                        , "CData",[cdata(:) cdata(:)] ...
                        %                                        , surface_kws_cell{:} ...
                        %                                        );
                        %    if self.is3d; view(3); end
                        %    hold on;
                        %end
                    end

                    %if self.isScatterPlot && self.scatter_kws.enabled
                    %    if self.colormap.enabled
                    %        if self.is3d
                    %            self.currentFig.scatter3 = scatter3 ( mdata ...
                    %                                                , vdata ...
                    %                                                , zdata ...
                    %                                                , self.scatter_kws.size ...
                    %                                                , cdata ...
                    %                                                , scatter_kws_cell{:} ...
                    %                                                );
                    %        else
                    %            self.currentFig.scatter = scatter   ( mdata ...
                    %                                                , vdata ...
                    %                                                , self.scatter_kws.size ...
                    %                                                , cdata ...
                    %                                                , scatter_kws_cell{:} ...
                    %                                                );
                    %        end
                    %    else
                    %        if self.is3d
                    %            %plot_kws = {};
                    %            %if ~isa(self.plot_kws,"cell"); plot_kws = self.plot_kws;
                    %            self.currentFig.scatter3 = scatter3 ( mdata ...
                    %                                                , vdata ...
                    %                                                , zdata ...
                    %                                                , self.scatter_kws.size ...
                    %                                                , currentScatterMarkerColor ...
                    %                                                , scatter_kws_cell{:} ...
                    %                                                );
                    %           %self.currentFig.plot = plot3( mdata ...
                    %           %                            , vdata ...
                    %           %                            , zdata ...
                    %           %                            , scatter_kws_cell{:} ...
                    %           %                            );
                    %        else
                    %            %plot_kws = {};
                    %            %if ~isa(self.plot_kws,"cell"); plot_kws = self.plot_kws;
                    %            self.currentFig.scatter = scatter   ( mdata ...
                    %                                                , vdata ...
                    %                                                , self.scatter_kws.size ...
                    %                                                , currentScatterMarkerColor ...
                    %                                                , scatter_kws_cell{:} ...
                    %                                                );
                    %           %self.currentFig.plot = plot ( mdata ...
                    %           %                            , vdata ...
                    %           %                            , scatter_kws_cell{:} ...
                    %           %                            );
                    %        end
                    %    end
                    %    hold on;
                    %end

                end % loop plot

                self.currentFig.gca = gca;
                if self.colormap.enabled
                    colormap(self.currentFig.gca,self.colormap.values);
                end

            end

            % add axis labels

            if mcolindexlen>1
                self.currentFig.xlabel = xlabel("Variable Values", "Interpreter", "none");
            else
                self.currentFig.xlabel = xlabel(mcolnames(1), "Interpreter", "none");
            end

            if vcolindexlen>1
                self.currentFig.ylabel = ylabel("Variable Values", "Interpreter", "none");
            else
                self.currentFig.ylabel = ylabel(vcolnames(1), "Interpreter", "none");
            end

            if self.is3d
            if zcolindexlen>1
                self.currentFig.zlabel = zlabel("Variable Values", "Interpreter", "none");
            else
                self.currentFig.zlabel = zlabel(zcolnames(1), "Interpreter", "none");
            end
            end

            % add line colorbar

            if self.colorbar_kws.enabled && self.colormap.enabled
                if isempty(self.colorbar_kws.fontsize) || ~isa(self.colorbar_kws.fontsize,"numeric")
                    self.colorbar_kws.fontsize = self.currentFig.ylabel.FontSize;
                end
                colorbar_kws_cell = convertStruct2Cell(self.colorbar_kws,{"enabled","singleOptions"});
                self.currentFig.colorbar = colorbar(colorbar_kws_cell{:});
                ylabel(self.currentFig.colorbar,ccolnames(1),"fontsize",self.colorbar_kws.fontsize, "Interpreter", "none");
            else
                colorbar('off');
                self.currentFig.colorbar = [];
            end

            if ~self.is3d || (self.is3d && self.legend_kws.enabled)
                self.doBasePlotStuff(self.legend_kws.enabled,lglabels)
            end

            box on; grid on; hold off;

        end % function plot

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Static)

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function bcrd = makeEllipsoid(covMat, meanVec, npoint) % returns the coordinates of the boundary of the ellipsoid
            if isempty(npoint); npoint = 50; end
            independentVariable = linspace(0, 2*pi, npoint)';
            xval = cos(independentVariable);
            yval = sin(independentVariable);
            ap = [xval(:) yval(:)]';
            [eigenVectors,eigenValues] = eig(covMat);
            eigenValues = sqrt(eigenValues); % convert variance to std
            bcrd = transpose(eigenVectors * eigenValues * ap + repmat(meanVec(:), 1, size(ap,2)));
            %h = plot(bcrd(:,1), bcrd(:,2), '-');
        end

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % classdef EllipsoidPlot