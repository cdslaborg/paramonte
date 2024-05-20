%
%   This is the Heatmap class for generating
%   instances of 2-dimensional Heatmap plots
%   based on the relevant MATLAB
%   intrinsic functions.
%
%       dfref
%
%           See the documentation of the corresponding input
%           argument of the superclass ``pm.vis.subplot.Subplot``.
%
%   \note
%       See the documentation of the attributes
%       of the superclass ``pm.vis.subplot.Subplot``.
%
%>  \return
%       An object of ``pm.vis.subplot.Heatmap`` class.
%>
%>  \interface{}
%>  \code{.m}
%>  \endcode
%>
%       s = pm.vis.subplot.Heatmap(dfref);
%       s = pm.vis.subplot.Heatmap(dfref, varargin);
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef Heatmap < pm.vis.subplot.Subplot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = Heatmap(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.subplot.Subplot("Heatmap", dfref, varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %
        %   Reset the heatmap colormap limits to the user-specified limits ``lb, ub``.
        %   If either is empty or missing, keep the existing limit for the corresponding empty limit.
        %   If both are empty or missing, symmetrize the existing limits.
        %
        %   The specified input values will be used to set the ``ColorLimits``
        %   component of the Heatmap object as ``ColorLimits = [lb, ub]``.
        %
        %   If both input values are missing, the current colormap
        %   range of the Heatmap plots will be symmetrized.
        %
        %   \warning
        %
        %       This method has side-effects by manipulating
        %       the existing attributes of the parent object.
        %
        %       lb
        %
        %           The input MATLAB object that can be either:
        %
        %               1.  a scalar of type double representing the lower bound of the Heatmap colormap range.
        %               1.  a vector of type double of size ``2`` representing the
        %                   lower and upper bounds of the Heatmap colormap range.
        %
        %           (**optional**. If missing, the current value will remain intact.)
        %
        %       ub
        %
        %           The input MATLAB scalar double representing the upper bound of the Heatmap colormap limits.
        %           Its value is completely ignored if the input ``lb`` argument is a vector of size ``2``.
        %           (**optional**. If missing, the current value will remain intact.)
        %>
        %>  \return
        %>
        %       None
        %
        %>  \interface{copy}
        %>  \code{.m}
        %>  \endcode
        %
        %       h = pm.vis.subplot.Subplot.make();
        %       h = pm.vis.subplot.Subplot.make([]);
        %       h = pm.vis.subplot.Subplot.make([], []);
        %       h = pm.vis.subplot.Subplot.make(lb, []);
        %       h = pm.vis.subplot.Subplot.make([], ub);
        %       h = pm.vis.subplot.Subplot.make(lb, ub);
        %
        %   Example
        %   -------
        %
        %       h = pm.vis.subplot.Heatmap(dfref);
        %       h.make()
        %       h.setcl() % symmetrize the current range.
        %       h.setcl(1) % set the lower bound to 1.
        %       h.setcl([], 1) % set the upper bound to 1.
        %       h.setcl(1, 2) % set the lower and upper bounds to 1 and 2.
        %       h.setcl([1, 2]) % set the lower and upper bounds to 1 and 2.
        %
        %>  \final{}
        %>
        %>  \author
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function setcl(self, lb, ub)

            if  isfield(self.fout, "heatmap") && isprop(self.fout.heatmap, "ColorLimits")
                if  nargin < 3
                    ub = [];
                end
                if  nargin < 2
                    lb = [];
                elseif 1 < length(lb)
                    ub = lb(2);
                    lb = lb(1);
                end
                if  isempty(lb) && isempty(ub)
                    % symmetrize the existing range.
                    maxval = max(abs(self.fout.heatmap.ColorLimits));
                    limits = [-maxval, maxval];
                else
                    limits = self.fout.heatmap.ColorLimits;
                    if ~isempty(lb)
                        limits(1) = lb;
                    end
                    if ~isempty(ub)
                        limits(2) = ub;
                    end
                end
                self.fout.heatmap.ColorLimits = limits;
            else
                error   ( newline ...
                        + "There is no component ``fout.heatmap.ColorLimits`` for this object to adjust the limits." + newline ...
                        + newline ...
                        );
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % classdef