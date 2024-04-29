classdef Plot < pm.vis.Figure
    %
    %   This is the abstract class for generating instances of objects
    %   that contain the specifications of various types of figures.
    %
    %   This is a generic class for generating figures
    %   containing a single subplot (axes).
    %
    %   Parameters
    %   ----------
    %
    %       subplot
    %
    %           The input cell matrix of MATLAB objects of superclass ``pm.vis.subplot.Subplot``.
    %           The input ``subplot`` object must minimally have the ``make()`` and ``reset()`` methods.
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
    %           of the ``subplot`` component of the parent object.
    %
    %   Returns
    %   -------
    %
    %       self
    %
    %           The output scalar object of class ``pm.vis.plot.Plot``.
    %
    %   Interface
    %   ---------
    %
    %       plot = pm.vis.plot.Plot(subplot);
    %
    %   Attributes
    %   ----------
    %
    %       See the list of class attributes below,
    %       also those of the superclass ``pm.vis.Figure``.
    %
    properties(Access = public)
        %
        %       subplot
        %
        %           The scalar object of superclass ``pm.vis.subplot.Subplot``
        %           representing the set of subplots to display in the figure.
        %
        subplot = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = Plot(subplot, varargin)
            varargin = {"subplot", subplot, varargin{:}};
            self = self@pm.vis.Figure(varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function reset(self, varargin)
            %
            %   Reset the properties of the plot to the original default settings.
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
            %           of the ``subplot`` component of the parent object.
            %
            %   Returns
            %   -------
            %
            %       None
            %
            %   Interface
            %   ---------
            %
            %       pm.vis.plot.Plot.reset() # reset the plot to the default settings.
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            [varfig, varsub] = pm.matlab.hashmap.popKeyVal(["figure", "subplot"], varargin);
            reset@pm.vis.Figure(self, varfig{:});
            self.subplot.reset(varsub{:});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: Any non-MATLAB-default setting must be preferably set in premake() method to override user null values.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function make(self, varargin)
            %
            %   Configure the plot settings and specifications,
            %   make the plot, and return nothing.
            %
            %   In making the figure, this method we call the ``make()``
            %   methods of each of the subplot objects stored in the
            %   ``subplot`` component.
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
            %           of the ``subplot`` component of the parent object.
            %
            %   Returns
            %   -------
            %
            %       None
            %
            %   Interface
            %   ---------
            %
            %       p = pm.vis.plot.Plot.make(varargin);
            %
            %   Example
            %   -------
            %
            %       p = pm.vis.plot.Plot(pm.vis.subplot.Line());
            %       p.make()
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            [varfig, varsub] = pm.matlab.hashmap.popKeyVal(["figure", "subplot"], varargin);
            make@pm.vis.Figure(self, varfig{:});
            self.subplot.make(varsub{:});
            hold off;

        end % function

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function premake(self, varargin)
            %
            %   Configure the plot settings and specifications and return nothing.
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
            %       f = pm.vis.plot.Plot.premake(varargin);
            %
            %   Example
            %   -------
            %
            %       f = pm.vis.plot.Plot(pm.vis.Line());
            %       f.premake("figure", {"color", "none"})
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            premake@pm.vis.Figure(self, varargin{:});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% These settings must happen here so that they can be reset every time user nullifies the values.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end