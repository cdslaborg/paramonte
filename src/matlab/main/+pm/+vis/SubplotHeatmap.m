%>  \brief
%>  This is the SubplotHeatmap class for generating
%>  instances of 2-dimensional Heatmap [Subplot visualizations](@ref Subplot)
%>  based on the relevant MATLAB
%>  intrinsic functions.<br>
%>
%>  \note
%>  See the documentation of the constructor of the class
%>  [pm.vis.SubplotHeatmap](@ref SubplotHeatmap::SubplotHeatmap) for example usage.<br>
%>
%>  \see
%>  [pm.vis.Cascade](@ref Cascade)<br>
%>  [pm.vis.Subplot](@ref Subplot)<br>
%>  [pm.vis.Figure](@ref Figure)<br>
%>  [pm.vis.Corner](@ref Corner)<br>
%>  [pm.vis.Plot](@ref Plot)<br>
%>  [pm.vis.Tile](@ref Tile)<br>
%>
%>  \final
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 6:05 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef SubplotHeatmap < pm.vis.Subplot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Construct and return an object of class [pm.vis.SubplotHeatmap](@ref SubplotHeatmap).<br>
        %>
        %>  \details
        %>  This is the constructor of the class [pm.vis.SubplotHeatmap](@ref SubplotHeatmap).<br>
        %>
        %>  \param[in]      dfref       :   See the documentation of the corresponding input
        %>                                  argument of the superclass [pm.vis.Subplot](@ref Subplot).<br>
        %>  \param[in]      varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                                  If the property is a ``struct()``, then its value must be given as a cell array,
        %>                                  with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                                  Note that all of these property-value pairs can be also directly set via the
        %>                                  parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \return
        %>  ``self``                    :   The output object of class [pm.vis.SubplotHeatmap](@ref SubplotHeatmap).<br>
        %>
        %>  \interface{SubplotHeatmap}
        %>  \code{.m}
        %>
        %>      s = pm.vis.SubplotHeatmap(dfref);
        %>      s = pm.vis.SubplotHeatmap(dfref, varargin);
        %>
        %>  \endcode
        %>
        %>  \note
        %>  See also the documentation of the attributes
        %>  of the superclass [pm.vis.Subplot](@ref Subplot).<br>
        %>
        %>  \example{SubplotHeatmap}
        %>  \include{lineno} example/vis/SubplotHeatmap/main.m
        %>  \vis{Subplot}
        %>  <br>\image html example/vis/SubplotHeatmap/SubplotHeatmap.1.png width=700
        %>
        %>  \final{SubplotHeatmap}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:05 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = SubplotHeatmap(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.Subplot("Heatmap", dfref, varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Reset the heatmap colormap limits to the user-specified limits ``lb`` and ``ub``.<br>
        %>  If either is empty or missing, keep the existing limit for the corresponding empty limit.<br>
        %>  If both are empty or missing, symmetrize the existing limits.<br>
        %>  The specified input values will be used to set the ``ColorLimits``
        %>  component of the Heatmap object as ``ColorLimits = [lb, ub]``.<br>
        %>  If both input values are missing, the current colormap
        %>  range of the Heatmap subplots will be symmetrized.<br>
        %>
        %>  \warning
        %>  This method has side-effects by manipulating
        %>  the existing attributes of the parent object.<br>
        %>
        %>  \param[in]  lb  :   The input MATLAB object that can be either:<br>
        %>                      <ol>
        %>                          <li>    a scalar of type double representing the lower bound of the Heatmap colormap range.<br>
        %>                          <li>    a vector of type double of size ``2`` representing the lower
        %>                                  and upper bounds of the Heatmap colormap range.<br>
        %>                      </ol>
        %>                      (**optional**. If missing, the current value will remain intact.)
        %>
        %>  \param[in]  ub  :   The input MATLAB scalar double representing the upper bound of the Heatmap colormap limits.<br>
        %>                      Its value is completely ignored if the input ``lb`` argument is a vector of size ``2``.<br>
        %>                      (**optional**. If missing, the current value will remain intact.)
        %>
        %>  \interface{setColorLim}
        %>  \code{.m}
        %>
        %>      h = pm.vis.Subplot.make();
        %>      h = pm.vis.Subplot.make([]);
        %>      h = pm.vis.Subplot.make([], []);
        %>      h = pm.vis.Subplot.make(lb, []);
        %>      h = pm.vis.Subplot.make([], ub);
        %>      h = pm.vis.Subplot.make(lb, ub);
        %>
        %>  \endcode
        %>
        %>  \example{setColorLim}
        %>  \code{.m}
        %>
        %>      h = pm.vis.SubplotHeatmap(dfref);
        %>      h.make()
        %>
        %>      h.setColorLim()       % symmetrize the current range.
        %>      h.setColorLim(1)      % set the lower bound to 1.
        %>      h.setColorLim([], 1)  % set the upper bound to 1.
        %>      h.setColorLim(1, 2)   % set the lower and upper bounds to 1 and 2.
        %>      h.setColorLim([1, 2]) % set the lower and upper bounds to 1 and 2.
        %>
        %>  \endcode
        %>
        %>  \example{setColorLim}
        %>  \include{lineno} example/vis/SubplotHeatmap/main.m
        %>  \vis{setColorLim}
        %>  <br>\image html example/vis/SubplotHeatmap/SubplotHeatmap.1.png width=700
        %>
        %>  \final{setColorLim}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 5:54 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function setColorLim(self, lb, ub)

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