%
%   This is the abstract class for generating instances
%   of figures containing a tile of subplots.
%
%       subplot
%
%           The input cell matrix of MATLAB objects of superclass ``pm.vis.subplot.Subplot``.
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
%           The output scalar object of class ``pm.vis.figure.Tiling``.
%>
%>  \interface{}
%>  \code{.m}
%>  \endcode
%>
%       plot = pm.vis.figure.Tiling(subplot);
%
%   \note
%       See the list of class attributes below,
%       also those of the superclass ``pm.vis.figure.Figure``.
%
classdef Tiling < pm.vis.figure.Figure

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = public)
        %
        %   tiledlayout
        %
        %       A MATLAB ``struct`` whose fields and values are passed
        %       as keyword arguments to the MATLAB intrinsic ``tiledlayout()``.
        %
        tiledlayout = [];
        %
        %   subplot
        %
        %       The MATLAB cell matrix containing objects of superclass ``pm.vis.subplot.Subplot``
        %       each of which represents one subplot axes to display in the figure.
        %
        subplot = cell(0, 0);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %properties(Access = public, Hidden)
    %    %
    %    %       ncol
    %    %
    %    %           The MATLAB scalar whole-number whose value represents ``size(self.suplot, 2)``.
    %    %
    %    ncol = [];
    %    %
    %    %       nrow
    %    %
    %    %           The MATLAB scalar whole-number whose value represents ``size(self.suplot, 1)``.
    %    %
    %    nrow = [];
    %end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = Tiling(subplot, varargin)
            if  nargin < 1
                subplot = cell(0, 0);
            end
            failed = ~iscell(subplot);
            if ~failed
                for irow = 1 : size(subplot, 1)
                    for icol = 1 : size(subplot, 2)
                        failed = ~isempty(subplot{irow, icol}) && ~pm.introspection.istype(subplot{irow, icol}, "pm.vis.subplot.Subplot");
                        if  failed
                            break;
                        end
                    end
                    if  failed
                        break;
                    end
                end
            end
            if ~failed
                varargin = {"subplot", subplot, varargin{:}};
            else
                help("pm.vis.figure.Tiling");
                error   ( newline ...
                        + "The input argument ``subplot`` must be a MATLAB cell matrix of " + newline ...
                        + "empty objects or objects of superclass ``pm.vis.subplot.Subplot``." + newline ...
                        + "For more information, see the class documentation displayed above." + newline ...
                        + newline ...
                        );
            end
            self = self@pm.vis.figure.Figure(varargin{:});
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
        %       varargin
        %
        %           Any ``property, value`` pair of the parent object.
        %           If the property is a ``struct()``, then its value must be given as a cell array,
        %           with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %           Note that all of these property-value pairs can be also directly set via the
        %           parent object attributes, before calling the ``make()`` method.
        %>
        %>  \return
        %>
        %       None
        %
        %>  \interface{copy}
        %>  \code{.m}
        %>  \endcode
        %
        %       f = pm.vis.figure.Tiling.make(varargin);
        %
        %   Example
        %   -------
        %
        %       f = pm.vis.figure.Tiling();
        %       f.make()
        %
        %>  \final{}
        %>
        %>  \author
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function make(self, varargin)

            make@pm.vis.figure.Figure(self, varargin{:});

            %%%% Resize the figure to allow good default visualization.

            if  isempty(self.figure.outerPosition) && 0 < numel(self.subplot)
                fullSize = get(0, 'ScreenSize');
                maxSize = fullSize;
                maxSize(1:2) = .05 * maxSize(3:4);
                maxSize(3:4) = maxSize(3:4) - maxSize(1:2);
                figSize = self.fout.figure.OuterPosition;
                maxScale = maxSize(3:4) ./ figSize(3:4);
                newWidth = figSize(3) * min(maxScale(1), size(self.subplot, 2));
                newHeight = figSize(4) * min(maxScale(2), size(self.subplot, 1));
                figStart = [(fullSize(3) - newWidth) / 2, (fullSize(4) - newHeight) / 2];
                set(self.fout.figure, "OuterPosition", [figStart(1:2), newWidth, newHeight]);
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: No component of ``self`` is allowed to appear to the left of assignment operator, except ``fout``.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            kws = struct();
            for prop =  [ "tiledlayout" ...
                        ]
                if  isprop(self, prop)
                    kws.(prop) = self.comp2hash(prop);
                end
            end

            try
                self.fout.tiledlayout = tiledlayout(size(self.subplot, 1), size(self.subplot, 2), kws.tiledlayout{:}); % requires MATLAB R2019b.
            end
            iplt = 0;
            timer = pm.timing.Timer();
            spinner = pm.timing.Spinner();
            for irow = 1 : size(self.subplot, 1)
                for icol = 1 : size(self.subplot, 2)
                    iplt = iplt + 1;
                    spinner.spin(iplt / numel(self.subplot));
                    if  pm.introspection.istype(self.subplot{irow, icol}, "pm.vis.subplot.Subplot")
                        try
                            nexttile;
                        catch
                            subplot(size(self.subplot, 1), size(self.subplot, 2), iplt);
                        end
                        self.subplot{irow, icol}.make();
                    end
                end
            end
            if  0 < iplt
                disp("done in " + sprintf("%.6f", string(timer.toc())) + " seconds.");
            end

        end % function

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %
        %   Reset the properties of the figure to the original default settings.
        %   Use this method when you change many attributes of the plot and
        %   you want to clean up and go back to the default settings.
        %
        %       varargin
        %
        %           Any ``property, value`` pair of the parent object.
        %           If the property is a ``struct()``, then its value must be given as a cell array,
        %           with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %           Note that all of these property-value pairs can be also directly set via the
        %           parent object attributes, before calling the ``make()`` method.
        %>
        %>  \return
        %>
        %       None
        %
        %>  \interface{copy}
        %>  \code{.m}
        %>  \endcode
        %
        %       pm.vis.figure.Tiling.reset() # reset all object properties to the default settings.
        %
        %>  \final{}
        %>
        %>  \author
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function reset(self, varargin)

            self.tiledlayout.innerPosition = [];
            self.tiledlayout.outerPosition = [];
            self.tiledlayout.position = [];
            self.tiledlayout.positionConstraint = [];
            self.tiledlayout.padding = [];
            self.tiledlayout.tileSpacing = [];
            self.tiledlayout.units = [];

            reset@pm.vis.figure.Figure(self, varargin{:}); % this will automatically call the ``premake()`` method of the object.

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: Any non-MATLAB-default setting must be preferably set in the make() method to override user null values.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %self.premake(varargin{:}); % This is the subclass method!

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %
        %   Update the layout of the tile with the new axes position changes.
        %
        %   \warning
        %
        %       This method causes side-effects by manipulating
        %       the existing attributes of the object.
        %
        %       varargin
        %
        %           Any ``property, value`` pair of the parent object.
        %           If the property is a ``struct()``, then its value must be given as a cell array,
        %           with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %           Note that all of these property-value pairs can be also directly set via the
        %           parent object attributes, before calling the ``make()`` method.
        %>
        %>  \return
        %>
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

            premake@pm.vis.figure.Figure(self, varargin{:});

            %self.nrow = size(self.subplot, 1);
            %self.ncol = size(self.subplot, 2);

            %%%% Set the default margins.

            self.setKeyVal("tiledlayout", "tileSpacing", "tight");
            %self.setKeyVal("tiledlayout", "padding", "tight");

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end