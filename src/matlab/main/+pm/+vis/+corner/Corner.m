%
%   This is the abstract class for generating instances
%   of figures containing a symmetric grid (corners) of subplots.
%
%       varargin
%
%           Any ``property, value`` pair of the parent object.
%           If the property is a ``struct()``, then its value must be given as a cell array,
%           with consecutive elements representing the struct ``property-name, property-value`` pairs.
%           Note that all of these property-value pairs can be also directly set via the
%           parent object attributes, before calling the ``make()`` method.
%
%>  \return
%       self
%
%           The output scalar object of class ``pm.vis.corner.Corner``.
%
%   Interface
%   ---------
%
%       plot = pm.vis.corner.Corner(subplot);
%
%   Attributes
%   ----------
%
%       See the list of class attributes below,
%       also those of the superclass ``pm.vis.figure.Tiling``.
%
classdef Corner < pm.vis.figure.Tiling

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = public, Hidden)
        % %
        % %   df
        % %
        % %       A scalar object of class ``pm.data.DataFrame``
        % %       containing the user-specified data to visualize.
        % %
        % df = [];
        %
        %   rows
        %
        %       A numeric vector that serves as a storage for an arbitrary subset of indices
        %       of the rows of the input dataframe reference ``dfref`` to the class constructor .
        %
        %       It can be either:
        %
        %           1.  a numeric range, or,
        %           2.  a list of row indices of the ``dfref``.
        %
        %       Example usage:
        %
        %           1.  rows = 15 : -2 : 8
        %           2.  rows = [12, 46, 7, 8, 9, 4, 7, 163]
        %
        %       If ``rows`` is empty, the default will be set by the template subplots
        %       specified in the components ``diag``, ``lower``, and ``upper`` of the
        %       parent corner object.
        %
        %rows = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = public)
        %
        %   cols (standing for columns)
        %
        %       The vector of MATLAB integers or names of of the input dataframe columns
        %       that determines the columns of the specified dataframe to visualize.
        %       It can have multiple forms:
        %
        %           1.  a numeric or cell array of column indices in the input ``dfref``.
        %           2.  a string or cell array of column names in ``dfref.Properties.VariableNames``.
        %           3.  a cell array of a mix of the above two.
        %           4.  a numeric range.
        %
        %       Example usage:
        %
        %           1.  self.cols = [7, 8, 9]
        %           2.  self.cols = ["sampleLogFunc", "sampleVariable1"]
        %           3.  self.cols = {"sampleLogFunc", 9, "sampleVariable1"}
        %           4.  self.cols = 7:9      # every column in the data frame starting from column #7 to #9
        %           5.  self.cols = 7:2:20   # every other column in the data frame starting from column #7 to #20
        %
        %       The ``i``th element of the specified ``cols`` will be used as value for:
        %
        %           1.  the ``i``th column of data to visualize on the x-axis in subplot ``subplot{i, i}.
        %           1.  the ``i``th column of data to visualize on the x-axis in subplot ``subplot{i, :}.
        %           1.  the ``i``th column of data to visualize on the y-axis in subplot ``subplot{:, i}.
        %
        %       If ``cols`` is empty, no visualizations will be made.
        %
        cols = [];
        %
        %       diag
        %
        %           The scalar object of superclass ``pm.vis.subplot.Subplot``
        %           representing the template of the diagonal subplots to display.
        %           Note that only the visualization properties of the template are used.
        %           The data properties of the template are set by the ``make()`` method.
        %           (**optional**. The default value is ``pm.vis.subplot.Histogram([])``.)
        %
        diag = [];
        %
        %       lower
        %
        %           The scalar object of superclass ``pm.vis.subplot.Subplot``
        %           representing the template of the lower-triangle subplots to display.
        %           The data properties of the template are set by the ``make()`` method.
        %           (**optional**. The default value is ``pm.vis.subplot.Contour([])``.)
        %
        lower = [];
        %
        %       upper
        %
        %           The scalar object of superclass ``pm.vis.subplot.Subplot``
        %           representing the template of the upper-triangle subplots to display.
        %           The data properties of the template are set by the ``make()`` method.
        %           (**optional**. The default value is ``pm.vis.subplot.LineScatter([])``.)
        %
        upper = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = Corner(diag, lower, upper, varargin)
            if  nargin < 3
                upper = [];
            end
            if  nargin < 2
                lower = [];
            end
            if  nargin < 1
                diag = [];
            end
            varargin = {"diag", diag, "lower", lower, "upper", upper, varargin{:}};
            self = self@pm.vis.figure.Tiling(cell(0, 0), varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %
        %   Configure the figure settings and specifications,
        %   make the figure and the subplots, and return nothing.
        %   The subplots are made by calling their ``make()`` methods.
        %
        %   \warning
        %
        %       This method has side-effects by manipulating
        %       the existing attributes of the parent object.
        %
        %   Parameters
        %   ----------
        %
        %       varargin
        %
        %           Any ``property, value`` pair of the parent object.
        %           If the property is a ``struct()``, then its value must be given as a cell array,
        %           with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %           Note that all of these property-value pairs can be also directly set via the
        %           parent object attributes, before calling the ``make()`` method.
        %
        %   Returns
        %   -------
        %
        %       None
        %
        %   Interface
        %   ---------
        %
        %       f = pm.vis.corner.Corner.make(varargin);
        %
        %   Example
        %   -------
        %
        %       f = pm.vis.corner.Corner();
        %       f.make()
        %
        %   LICENSE
        %   -------
        %
        %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
        %
        function make(self, varargin)

            self.premake(varargin{:});

            %%%% First set the tile dimensions.

            if  isempty(self.cols)
                warning ( "The ``cols`` component of the parent corner object is empty." + newline ...
                        + "There is nothing to visualize." + newline ...
                        );
                return;
            else
                rank = length(self.cols);
            end

            %%%% Define the subplot cell matrix.
            %%%% The following code block may be improved in
            %%%% the future to avoid full data copy to subplots.

            iplt = 0;
            self.subplot = cell(rank, rank);
            for irow = 1 : rank
                for icol = 1 : rank
                    iplt = iplt + 1;
                    if  icol == irow
                        copyStream = getByteStreamFromArray(self.diag);
                    elseif  icol < irow
                        copyStream = getByteStreamFromArray(self.lower);
                    elseif  irow < icol
                        copyStream = getByteStreamFromArray(self.upper);
                    end
                    self.subplot{irow, icol} = getArrayFromByteStream(copyStream);
                    isEllipsoid = isa(self.subplot{irow, icol}, "pm.vis.subplot.Ellipse") || isa(self.subplot{irow, icol}, "pm.vis.subplot.Ellipse3");

                    if  isEllipsoid
                        self.subplot{irow, icol}.dims = [self.cols(irow), self.cols(icol)];
                    else
                        if  isprop(self.subplot{irow, icol}, "colx")
                            self.subplot{irow, icol}.colx = self.cols(icol);
                        end
                        if  isprop(self.subplot{irow, icol}, "coly")
                            self.subplot{irow, icol}.coly = self.cols(irow);
                        end

                        if  1 < lencolc
                            self.subplot{irow, icol}.colc = self.template.colc(iplt);
                        elseif lencolc == 1 && isprop(self.template, "colorbar") && isempty(self.template.colorbar.enabled)
                            self.subplot{irow, icol}.colorbar.enabled = false;
                        elseif lencolc == 0 && isprop(self.template, "colormap") && isempty(self.template.colormap.enabled)
                            self.subplot{irow, icol}.colormap.enabled = true;
                        end
                    end
                end
            end

            make@pm.vis.figure.Tiling(self);

            %%%% Define a single colorbar.

            if ~isEllipsoid

                unicbar = lencolc == 1;
                if  unicbar && ~isempty(self.template.colormap.enabled)
                    unicbar = self.template.colormap.enabled;
                end
                if  unicbar && ~isempty(self.template.colorbar.enabled)
                    unicbar = self.template.colorbar.enabled;
                end
                if  unicbar

                    %%%% Get the start and end positions of the leftmost, lowest, rightmost, and highest axes.

                    % iplt = 0;
                    % positions = zeros(4, nplt);
                    % for icol = 1 : ncol
                    %     for irow = 1 : nrow
                    %         if  pm.introspection.istype(self.subplot{irow, icol}, "pm.vis.subplot.Subplot")
                    %             iplt = iplt + 1;
                    %             positions(:, iplt) = self.subplot{irow, icol}.fout.axes.Position;
                    %         end
                    %     end
                    % end
                    %
                    % for i = 2 : -1 : 1
                    %     start(i) = min(positions(i, :));
                    %     finit(i) = max(positions(i, :) + positions(i + 2, :));
                    % end

                    %%%% Create new invisible axes.

                    kws = struct();
                    for prop =  [ "colorbar" ...
                                ]
                        if  isprop(self.template, prop)
                            kws.(prop) = self.template.comp2hash(prop);
                        end
                    end
                    %ax = axes("position", [start, finit], "visible", "off");
                    self.fout.colorbar = colorbar(kws.colorbar{:});
                    self.fout.colorbar.Layout.Tile = 'east';

                    dfcopy = self.template.df.copy();
                    [~, colnamc] = pm.str.locname(dfcopy.Properties.VariableNames, self.template.colc);
                    ylabel(self.fout.colorbar, colnamc(1));

                end

            end

        end % function

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %
        %   Reset the properties of the figure to the original default settings.
        %   Use this method when you change many attributes of the plot and
        %   you want to clean up and go back to the default settings.
        %
        %   Parameters
        %   ----------
        %
        %       varargin
        %
        %           Any ``property, value`` pair of the parent object.
        %           If the property is a ``struct()``, then its value must be given as a cell array,
        %           with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %           Note that all of these property-value pairs can be also directly set via the
        %           parent object attributes, before calling the ``make()`` method.
        %
        %   Returns
        %   -------
        %
        %       None
        %
        %   Interface
        %   ---------
        %
        %       pm.vis.corner.Corner.reset() # reset all object properties to the default settings.
        %
        %   LICENSE
        %   -------
        %
        %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
        %
        function reset(self, varargin)

            self.dfref = [];
            self.cols = [];
            self.rows = [];

            reset@pm.vis.figure.Tiling(self, varargin{:});

            if  isempty(self.diag)
                self.diag = pm.vis.subplot.Histogram([]);
            end
            if  isempty(self.lower)
                self.diag = pm.vis.subplot.Contour([]);
            end
            if  isempty(self.upper)
                self.diag = pm.vis.subplot.LineScatter([]);
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: Any non-MATLAB-default setting must be preferably set in the make() method to override user null values.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %self.premake(varargin{:}); % This is the subclass method!
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %
        %   Prepare the corner plot specs for subsequent visualization by ``make()``.
        %
        %   \warning
        %
        %       This method causes side-effects by manipulating
        %       the existing attributes of the object.
        %
        %   Parameters
        %   ----------
        %
        %       varargin
        %
        %           Any ``property, value`` pair of the parent object.
        %           If the property is a ``struct()``, then its value must be given as a cell array,
        %           with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %           Note that all of these property-value pairs can be also directly set via the
        %           parent object attributes, before calling the ``make()`` method.
        %
        %   Returns
        %   -------
        %
        %       None.
        %
        %   Example
        %   -------
        %
        %       For example, change the left/bottom margin of the main
        %       axis of the figure to provide room for lengthy variable names.
        %       Then call the ``self.update()`` method to reflect the changes.
        %
        function premake(self, varargin)

            premake@pm.vis.figure.Tiling(self, varargin{:});
            % if ~isa(self.diag, "self.type.is.histfit") && ~isa(self.diag, "self.type.is.histogram")
            %     help("pm.vis.corner.Corner");
            %     error   ( newline ...
            %             + "The component ``diag`` of the parent object must be of" + newline ...
            %             + "superclass ``pm.vis.subplot.Histfit`` or ``pm.vis.subplot.Histogram``." + newline ...
            %             + "For more information, see the class documentation displayed above." + newline ...
            %             + newline ...
            %             );
            % end
            %
            % if  self.lower.type.is.d1
            %     help("pm.vis.corner.Corner");
            %     error   ( newline ...
            %             + "The component ``lower`` of the parent object must be of superclass" + newline ...
            %             + "    1.  ``pm.vis.subplot.Line``" + newline ...
            %             + "    1.  ``pm.vis.subplot.Scatter``" + newline ...
            %             + "    1.  ``pm.vis.subplot.LineScatter``" + newline ...
            %             + "    1.  ``pm.vis.subplot.LineScatter``" + newline ...
            %             + "For more information, see the class documentation displayed above." + newline ...
            %             + newline ...
            %             );
            % end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end