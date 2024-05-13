classdef Cascade < pm.matlab.Handle
    %
    %   This is the abstract class for generating instances of objects
    %   that contain the specifications of a cascade of plots.
    %
    %   This is a generic class for generating multiple separate
    %   cascading figures each containing a single plot (axes).
    %
    %   Parameters
    %   ----------
    %
    %       template
    %
    %           The input scalar object of superclass ``pm.vis.cascade.Cascade``.
    %           It serves as the template based upon which
    %           the cascade of plots are constructed.
    %
    %       varargin
    %
    %           Any ``property, value`` pair of the parent object.
    %           If the property is a ``struct()``, then its value must be given as a cell array,
    %           with consecutive elements representing the struct ``property-name, property-value`` pairs.
    %           Note that all of these property-value pairs can be also directly set via the
    %           parent object attributes, before calling the ``make()`` method.
    %
    %           \note
    %
    %               The input ``varargin`` can also contain the components
    %               of the ``plot`` component of the parent object.
    %
    %   Returns
    %   -------
    %
    %       self
    %
    %           The output scalar object of class ``pm.vis.cascade.Cascade``.
    %
    %   Interface
    %   ---------
    %
    %       plot = pm.vis.cascade.Cascade(template);
    %
    %   Attributes
    %   ----------
    %
    %       See the list of class attributes below,
    %       also those of the superclass ``pm.vis.figure.Handle``.
    %
    properties(Access = public)
        %
        %       window
        %
        %           The scalar object of superclass ``pm.vis.cascade.Cascade``
        %           representing the cascade of plots to display.
        %
        window = [];
        %
        %       template
        %
        %           The scalar object of superclass ``pm.vis.cascade.Cascade``
        %           representing the template of the cascade of plots to display.
        %
        template = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = Cascade(template, varargin)
            self.template = template;
            self.reset(varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function reset(self, varargin)
            %
            %   Reset the properties of the cascade of plots to the original default settings.
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
            %       \note
            %
            %           The input ``varargin`` can also contain the components
            %           of the ``template`` component of the parent object.
            %
            %   Returns
            %   -------
            %
            %       None
            %
            %   Interface
            %   ---------
            %
            %       pm.vis.cascade.Cascade.reset() # reset the plot to the default settings.
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            if ~isempty(varargin)
                [varpar, varsub] = pm.matlab.hashmap.popKeyVal(["template", "window"], varargin);
                if ~isempty(varpar)
                    self.hash2comp(varpar);
                end
                if ~isempty(varsub)
                    self.template.reset(varsub{:});
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: Any non-MATLAB-default setting must be preferably set in premake() method to override user null values.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function make(self, varargin)
            %
            %   Configure the cascade settings and specifications,
            %   make the cascade of plots, and return nothing.
            %
            %   In making the figure, this method we call the ``make()``
            %   methods of each of the plot objects stored in the
            %   ``window`` component.
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
            %       \note
            %
            %           The input ``varargin`` can also contain the components
            %           of the ``template`` component of the parent object.
            %
            %   Returns
            %   -------
            %
            %       None
            %
            %   Interface
            %   ---------
            %
            %       p = pm.vis.cascade.Cascade.make(varargin);
            %
            %   Example
            %   -------
            %
            %       p = pm.vis.cascade.Cascade("coly", 1:3);
            %       p.make()
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            self.premake(varargin{:});

            %%%% First set the cascade dimensions.

            isEllipsoid = isa(self.template, "pm.vis.plot.Ellipse") || isa(self.template, "pm.vis.plot.Ellipse3");

            if  isEllipsoid
                nplt = size(self.template.subplot.dims, 1);
                if  size(self.template.subplot.dims, 2) ~= 2
                    help("pm.vis.plot.Ellipse3");
                    error   ( newline ...
                            + "The condition ``size(self.template.dims, 2) == 2`` must hold." + newline ...
                            + "For more information, see the documentation displayed above." + newline ...
                            + newline ...
                            );
                end

            else

                lencolx = 0;
                if  isprop(self.template.subplot, "colx")
                    lencolx = pm.array.len(self.template.subplot.colx);
                end

                lencoly = 0;
                if  isprop(self.template.subplot, "coly")
                    lencoly = pm.array.len(self.template.subplot.coly);
                end

                lencolz = 0;
                if  isprop(self.template.subplot, "colz")
                    lencolz = pm.array.len(self.template.subplot.colz);
                end

                lencolc = 0;
                if  isprop(self.template.subplot, "colc")
                    lencolc = pm.array.len(self.template.subplot.colc);
                end

                nplt = max([lencolx, lencoly, lencolz, lencolc]);

            end

            %%%% Define the cascade.
            %%%% The following code block may be improved in
            %%%% the future to avoid full data copy to subplots.

            self.window = cell(nplt);
            for iplt = 1 : nplt
                copyStream = getByteStreamFromArray(self.template);
                self.window{iplt} = getArrayFromByteStream(copyStream);
                if  isEllipsoid
                    self.window{iplt}.subplot.dims = self.template.subplot.dims(iplt, :);
                else
                    if  1 < lencolx
                        self.window{iplt}.subplot.colx = self.template.subplot.colx(iplt);
                    end
                    if  1 < lencoly
                        self.window{iplt}.subplot.coly = self.template.subplot.coly(iplt);
                    end
                    if  1 < lencolz
                        self.window{iplt}.subplot.colz = self.template.subplot.colz(iplt);
                    end
                    if  1 < lencolc
                        self.window{iplt}.subplot.colc = self.template.subplot.colc(iplt);
                    end
                end
                self.window{iplt}.make();
            end

        end % function

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function savefig(self, file, varargin)
            %
            %   Export the current cascade of figures to the specified external files.
            %
            %   \note
            %
            %       This method is merely a wrapper around
            %       all ``savefig`` methods of the figure objects in
            %       the ``window`` component of the parent cascade object.
            %
            %   Parameters
            %   ----------
            %
            %       file
            %
            %           The input vector of MATLAB strings or cell array of char vectors that must be
            %           of the same length as the length of the ``window`` component of the parent object.
            %           containing the paths to the external files to contain the visualizations.
            %           For more information, see the corresponding argument of the ``savefig``
            %           method of class ``pm.vis.figure.Figure``.
            %           (**optional**.  If ``file`` or any elements of it are is missing or empty,
            %           the default will be set by the ``savefig`` method of the corresponding
            %           cascade figure in the ``window`` component.)
            %
            %       varargin
            %
            %           Any optional key-val pair that is accepted by
            %           the ``savefig`` method of class ``pm.vis.figure.Figure``.
            %           For more information, see the ``savefig`` method of class ``pm.vis.figure.Figure``.
            %
            %   Returns
            %   -------
            %
            %       None.
            %
            %   Interface
            %   ---------
            %
            %       c.savefig();
            %       c.savefig(file);
            %       c.savefig(file, varargin{:});
            %
            %   Example
            %   -------
            %
            %       c = pm.vis.cascade.Cascade(pm.vis.plot.Line());
            %       c.savefig(); % export the current figure with the default name.
            %       c.savefig("gridplot.pdf") % export figure to the specified PDF file.
            %       c.savefig("gridplot.png", "-m4 -transparent") % export a large png plot of magnitude 4 with transparency.
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            if  nargin < 2 || isempty(file)
                file = strings(length(self.window), 1);
            end
            if  length(file) ~= length(self.window)
                help("pm.vis.cascade.Cascade");
                error   ( newline ...
                        + "The condition ``length(file) == length(self.window)`` must hold." + newline ...
                        + "For more information, see the documentation displayed above." + newline ...
                        + newline ...
                        );
            end
            for iwin = 1 : length(self.window)
                self.window{iwin}.savefig(file(iwin), varargin{:});
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function premake(self, varargin)
            %
            %   Configure the cascade of plots settings and specifications and return nothing.
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
            %           parent object attributes, before calling the ``premake()`` method.
            %
            %   Returns
            %   -------
            %
            %       None
            %
            %   Interface
            %   ---------
            %
            %       c = pm.vis.cascade.Cascade.premake(varargin);
            %
            %   Example
            %   -------
            %
            %       c = pm.vis.cascade.Cascade();
            %       c.premake("figure", {"color", "none"})
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            if ~isempty(varargin)
                [varpar, varsub] = pm.matlab.hashmap.popKeyVal(["template", "window"], varargin);
                if ~isempty(varpar)
                    self.hash2comp(varpar);
                end
                if ~isempty(varsub)
                    self.template.premake(varsub{:});
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% These settings must happen here so that they can be reset every time user nullifies the values.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end