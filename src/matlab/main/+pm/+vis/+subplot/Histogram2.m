%
%   This is the Histogram2 class for generating
%   instances of 3-dimensional Histogram2 plots
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
%       An object of ``pm.vis.subplot.Histogram2`` class.
%
%   Interface
%   ---------
%
%       s = pm.vis.subplot.Histogram2(dfref);
%       s = pm.vis.subplot.Histogram2(dfref, varargin);
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef Histogram2 < pm.vis.subplot.Subplot
    methods(Access = public)
        function self = Histogram2(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.subplot.Subplot("Histogram2", dfref, varargin{:});
        end
    end
end