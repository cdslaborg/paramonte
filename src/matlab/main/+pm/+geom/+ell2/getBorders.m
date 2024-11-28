%>  \brief
%>  Return a matrix of MATLAB doubles of shape ``[npnt, 2 * nell]``
%>  (or ``[npnt, 3 * nell]`` if the input ``zval`` is present)
%>  containing the 2D (or 3D) coordinates of a set of points on the boundary
%>  of a set of 2D ellipsoids whose Gramian matrices and centers are specified
%>  as input argument.<br>
%>
%>  \details
%>  Here ``nell`` refers to the number of ellipsoids which
%>  is equal to ``2 * max(1, size(gramian, 3), size(center, 2))``.
%>
%>  \param[in]  gramian :   The input matrix of MATLAB doubles of shape ``[2, 2, nell]``
%>                          containing the Gramian matrices of the ``nell`` 2D ellipsoids
%>                          whose boundary points are to be returned.<br>
%>                          (**optional**, default = ``eye(2, 2, 1)``)
%>  \param[in]  center  :   The input matrix of MATLAB doubles of shape ``[2, nell]`` containing
%>                          the 2D coordinates of the centers of the ``nell`` 2D ellipsoids
%>                          whose boundary points are to be returned.<br>
%>                          (**optional**, default = ``zeros(2, 1)``)
%>  \param[in]  zval    :   The input scalar (or matrix of shape ``[npnt, nell]`` of) MATLAB double(s)
%>                          containing the z-axis coordinates of the points on the borders of the ``nell`` 2D ellipsoids.<br>
%>                          If present, the specified value(s) will occupy the ``[1 : 2 : 3 * nell]`` columns of the output ``bcrd``.<br>
%>                          This argument is present for the sake of convenience and better performance of the algorithm to avoid further reallocations.<br>
%>                          (**optional**. If not present or empty, it will not be present in the output.)
%>  \param[in]  npnt    :   The input scalar MATLAB whole number containing the number of points to return on the boundary of the target 2D ellipsoid.<br>
%>                          (**optional**, default = ``50``)
%>
%>  \return
%>  ``bcrd``            :   The output matrix of MATLAB doubles of shape ``[npnt, 2 * nell]``
%>                          (or ``[npnt, 3 * nell]`` if the input argument ``zval`` is present)
%>                          containing the coordinates of a set of ``npnt`` points on
%>                          the boundary of the target 2D ellipsoid.
%>
%>  \interface{getBorders}
%>  \code{.m}
%>
%>      pm.geom.ell2.getBorders()
%>      pm.geom.ell2.getBorders(gramian)
%>      pm.geom.ell2.getBorders(gramian, center)
%>      pm.geom.ell2.getBorders(gramian, center, zval)
%>      pm.geom.ell2.getBorders(gramian, center, zval, npnt)
%>
%>  \endcode
%>
%>  \warning
%>  The condition ``npnt == size(zval, 1)`` must hold for the corresponding input arguments.
%>
%>  \example{getBorders}
%>  \include{lineno} example/geom/ell2/getBorders/main.m
%>  \vis{getBorders}
%>  \image html example/geom/ell2/getBorders/getBorders.2d.png width=700
%>  \image html example/geom/ell2/getBorders/getBorders.3d.png width=700
%>
%>  \final{getBorders}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 5:33 PM, University of Texas at Arlington<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function bcrd = getBorders(gramian, center, zval, npnt)

    if  nargin < 4
        npnt = [];
    end
    if  nargin < 3
        zval = [];
    end
    if  nargin < 2
        center = [];
    end
    if  nargin < 1
        gramian = [];
    end
    if  isempty(npnt)
        if  1 < size(zval, 1)
            npnt = size(zval, 1);
        else
            npnt = 50;
        end
    end
    if  isempty(center)
        center = zeros(2, 1);
    end
    if  isempty(gramian)
        gramian = zeros(2, 2, 1);
        gramian(:, :, 1) = eye(2);
    end

    nell = max(size(gramian, 3), size(center, 2));
    asserted = size(gramian, 3) == size(center, 2) || size(gramian, 3) == 1 || size(center, 2) == 1;
    if ~asserted
        help("pm.geom.ell2.getBorders")
        disp("size(gramian)")
        disp( size(gramian) )
        disp("size(center)")
        disp( size(center) )
        error   ( newline ...
                + "The condition ``size(gramian, 3) == size(center, 2) || size(gramian, 3) == 1 || size(center, 2) == 1`` must hold." + newline ...
                + "For more information, see the documentation displayed above." + newline ...
                + newline ...
                );
    end

    %%%%
    %%%% Make the 2D ellipsoid boundaries.
    %%%%

    if  isempty(zval)
        bcrd = zeros(npnt, 2);
        for iell = 1 : nell
            icol = (iell - 1) * 2 + 1;
            bcrd(:, icol : icol + 1) = pm.geom.ell2.getBorder(gramian(:, :, min(iell, size(gramian, 3))), center(:, min(iell, size(center, 2))), npnt);
        end
    else
        asserted = nell == size(zval, 2) || size(zval, 2) == 1 || (size(gramian, 3) == 1 || size(center, 2) == 1);
        if ~asserted
            help("pm.geom.ell2.getBorders")
            disp("size(gramian, 3)")
            disp( size(gramian, 3) )
            disp("size(canter, 2)")
            disp( size(canter, 2) )
            disp("size(zval)")
            disp( size(zval) )
            error   ( newline ...
                    + "The specified number of ellipsoids must be either the same or one in the input ``gramian``, ``center``, and ``zval``." + newline ...
                    + "For more information, see the documentation displayed above." + newline ...
                    + newline ...
                    );
        end
        zvalScalar = numel(zval) == 1;
        nell = max(nell, size(zval, 2));
        bcrd = zeros(npnt, 3 * nell);
        for iell = 1 : nell
            icol = (iell - 1) * 3 + 1;
            bcrd(:, icol : icol + 1) = pm.geom.ell2.getBorder(gramian(:, :, min(iell, size(gramian, 3))), center(:, min(iell, size(center, 2))), npnt);
            if  zvalScalar
                bcrd(:, icol + 2) = zval;
            else
                bcrd(:, icol + 2) = zval(:, min(iell, size(zval, 2)));
            end
        end
    end

end