classdef Ellipse < pm.vis.subplot.Ellipse3
    %
    %   This is the Ellipse class for generating
    %   instances of 2-dimensional Ellipse plots
    %   based on the relevant MATLAB
    %   intrinsic functions.
    %
    %   In the following documentation,
    %
    %       1.  the variable ``ndim`` represents the number of ellipsoids in the input data.
    %           The value of ``ndim`` is inferred from the shapes
    %           of the input ``gramian`` and ``center`` arguments
    %           as ``max([2, size(gramian, 1), size(center, 1)])``.
    %
    %       2.  the variable ``nell`` represents the number of ellipsoids in the input data.
    %           The value of ``nell`` is inferred from the shapes
    %           of the input ``gramian``, ``center``, and ``cval`` arguments
    %           as ``max([size(gramian, 3), size(center, 2), size(cval, 2)])``.
    %           If the above expression yields zero, then ``nell`` is set to ``75``.
    %
    %       3.  the variable ``npnt`` represents the number of points used in visualizing each ellipsoid.
    %           The value of ``npnt`` is inferred from the shapes of the input argument ``cval`` as ``size(cval, 1)``.
    %           If the above expression yields zero or one, then ``npnt`` is set to ``100``.
    %
    %   \note
    %
    %       This class is merely a simple extension of
    %       the superclass ``pm.vis.subplot.Ellipse3``.
    %       The only difference is that all ``zval`` values for
    %       all points of all coordinates are set to default values.
    %       The only tangible action this class performs is
    %       to adjust the camera view to a 2D axis.
    %
    %   Parameters
    %   ----------
    %
    %       gramian
    %
    %           See the corresponding input argument to the constructor of the parent object.
    %
    %       center
    %
    %           See the corresponding input argument to the constructor of the parent object.
    %
    %       cval
    %
    %           See the corresponding input argument to the constructor of the parent object.
    %
    %       varargin
    %
    %           Any ``property, value`` pair of the object.
    %           If the property is a ``struct()``, then its value must be given as a cell array,
    %           with consecutive elements representing the struct ``property-name, property-value`` pairs.
    %           Note that all of these property-value pairs can be also directly set via the
    %           parent object attributes, before calling the ``make()`` method.
    %
    %   Attributes
    %   ----------
    %
    %       See below and also the documentation of the attributes
    %       of the superclass ``pm.vis.subplot.Ellipse3``.
    %
    %   Returns
    %   -------
    %
    %       An object of class ``pm.vis.subplot.Ellipse``.
    %
    %   Interface
    %   ---------
    %
    %       s = pm.vis.subplot.Ellipse();
    %       s = pm.vis.subplot.Ellipse(gramian);
    %       s = pm.vis.subplot.Ellipse(gramian, center);
    %       s = pm.vis.subplot.Ellipse(gramian, center, cval);
    %       s = pm.vis.subplot.Ellipse(gramian, center, cval, varargin);
    %
    %   Example
    %   -------
    %
    %       s = pm.vis.subplot.Ellipse();
    %       s.make("dims", [1, 2]);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    methods(Access = public)
        function self = Ellipse(gramian, center, cval, varargin)

            %%%% Define the missing optional values as empty with the right rank.

            if  nargin < 3
                cval = zeros(0, 0);
            end
            if  nargin < 2
                center = zeros(0, 0);
            end
            if  nargin < 1
                gramian = zeros(0, 0, 0);
            end
            self = self@pm.vis.subplot.Ellipse3(gramian, center, [], cval, varargin{:})

        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)
        function make(self, varargin)
            %
            %   Generate a plot from the selected dimensions
            %   of the input ellipsoid data to the object constructor.
            %
            %   \note
            %
            %       This method is pure, except for the changes in ``fout`` component.
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
            %       pm.vis.subplot.Ellipse.make(varargin);
            %
            %   Example
            %   -------
            %
            %       s = pm.vis.subplot.Ellipse();
            %       s.make("dims", [1, 2]);
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            make@pm.vis.subplot.Ellipse3(self, varargin{:});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: No component of ``self`` is allowed to appear to the left of assignment operator, except ``fout``.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            view(2);

        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end