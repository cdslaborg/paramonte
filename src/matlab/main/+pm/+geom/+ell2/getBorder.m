%>  \brief
%>  Return a matrix of MATLAB doubles of shape ``[npnt, 2]``
%>  containing the coordinates of a set of points on the boundary
%>  of a 2D ellipsoid whose Gramian matrix is specified as input
%>  and whose center is also optionally specified.<br>
%>
%>  \param[in]  gramian :   The input square matrix of MATLAB doubles of shape ``[2, 2]`` containing the
%>                          Gramian of the target 2D ellipsoid whose boundary points are to be returned.<br>
%>                          (**optional**. If not present or empty, the default is ``eye(2, 2)``.)
%>  \param[in]  center  :   The input vector of MATLAB doubles of size ``2`` containing the 2D coordinates
%>                          of the center of the target 2D ellipsoid whose boundary points are to be returned.<br>
%>                          (**optional**. If not present or empty, the default is ``zeros(2, 1)``.)
%>  \param[in]  npnt    :   The input scalar MATLAB whole number containing the number of points to
%>                          return on the boundary of the target 2D ellipsoid.<br>
%>                          (**optional**, default = ``50``.)
%>
%>  \return
%>  `bcrd`              :   The output matrix of MATLAB doubles of shape ``[npnt, 2]``
%>                          containing the coordinates of a set of ``npnt`` points on
%>                          the boundary of the target 2D ellipsoid.<br>
%>
%>  \interface{getBorder}
%>  \code{.m}
%>
%>      bcrd = pm.geom.ell2.getBorder();
%>      bcrd = pm.geom.ell2.getBorder(gramian);
%>      bcrd = pm.geom.ell2.getBorder(gramian, center);
%>      bcrd = pm.geom.ell2.getBorder(gramian, center, npnt);
%>      bcrd = pm.geom.ell2.getBorder([], [], []);
%>
%>  \endcode
%>
%>  \see
%>  [pm.geom.ell2.getBorders](@ref getBorders)<br>
%>  [pm.vis.cascade.Ellipse](@ref Ellipse)<br>
%>  [pm.vis.subplot.Ellipse](@ref Ellipse)<br>
%>  [pm.vis.plot.Ellipse](@ref Ellipse)<br>
%>  [pm.vis.tile.Ellipse](@ref Ellipse)<br>
%>
%>  \example{getBorder}
%>  \include{lineno} example/geom/ell2/getBorder/main.m
%>  \vis{getBorder}
%>  \image html example/geom/ell2/getBorder/getBorder.2d.png width=700
%>  \image html example/geom/ell2/getBorder/getBorder.3d.png width=700
%>  \image html example/geom/ell2/getBorder/getBorder.wavy.png width=700
%>
%>  \final{getBorder}
%>
%>  \author
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function bcrd = getBorder(gramian, center, npnt)
    if  nargin < 3
        npnt = [];
    end
    if  nargin < 2
        center = [];
    end
    if  nargin < 1
        gramian = [];
    end
    if  isempty(npnt)
        npnt = 50;
    end
    if  isempty(center)
        center = zeros(2, 1);
    end
    if  isempty(gramian)
        gramian = eye(2);
    end
    independentVariable = linspace(0, 2 * pi, npnt)';
    xval = cos(independentVariable);
    yval = sin(independentVariable);
    ap = [xval(:) yval(:)]';
    [eigvec, eigval] = eig(gramian);
    eigval = sqrt(eigval); % convert variance to std.
    bcrd = transpose(eigvec * eigval * ap + repmat(center(:), 1, npnt));
end