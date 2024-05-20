%
%   This is the Line class for generating
%   instances of 2-dimensional Line plots
%   based on the relevant MATLAB
%   intrinsic functions.
%
%       dfref
%
%           See the documentation of the corresponding input
%           argument of the superclass ``pm.vis.subplot.Subplot``.
%
%   Attributes
%   ----------
%
%       See the documentation of the attributes
%       of the superclass ``pm.vis.subplot.Subplot``.
%
%>  \return
%       An object of ``pm.vis.subplot.Line`` class.
%
%   Interface
%   ---------
%
%       s = pm.vis.subplot.Line(dfref);
%       s = pm.vis.subplot.Line(dfref, varargin);
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef Line < pm.vis.subplot.Subplot
    methods(Access = public)
        function self = Line(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.subplot.Subplot("Line", dfref, varargin{:});
        end
    end
end