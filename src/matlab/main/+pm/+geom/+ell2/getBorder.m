function bcrd = getBorder(gramian, center, zval, npnt)
    %
    %   Return a matrix of MATLAB doubles of shape ``[npnt, 2]``
    %   (or of shape ``[npnt, 3]`` if the input argument ``zval`` is present)
    %   containing the coordinates of a set of points on the boundary
    %   of a 2D ellipsoid whose Gramian matrix is specified as input
    %   and whose center is also optionally specified.
    %
    %   The third column of the output matrix will contain the
    %   specified input ``zval`` or otherwise, will not exit.
    %
    %   Parameters
    %   ----------
    %
    %       gramian
    %
    %           The input square matrix of MATLAB doubles of shape ``[2, 2]``
    %           containing the Gramian of the target 2D ellipsoid whose
    %           boundary points are to be returned.
    %           (**optional**. If not present or empty, the default is ``eye(2, 2)``.)
    %
    %       center
    %
    %           The input vector of MATLAB doubles of size ``2`` containing
    %           the 2D coordinates of the center of the target 2D ellipsoid
    %           whose boundary points are to be returned.
    %           (**optional**. If not present or empty, the default is ``zeros(2, 1)``.)
    %
    %       zval
    %
    %           The input scalar (or vector of size ``npnt`` of) MATLAB double(s) containing
    %           the z-axis coordinates of the points of the output 2D ellipsoid border.
    %           If present, the specified value(s) will occupy the third column of the output ``bcrd``.
    %           (**optional**. If not present or empty, it will not be present in the output.)
    %
    %           \note
    %
    %               This argument is present for the sake of convenience and
    %               better performance of the algorithm to avoid further reallocations.
    %
    %       npnt
    %
    %           The input scalar MATLAB whole number containing the number of
    %           points to return on the boundary of the target 2D ellipsoid.
    %           (**optional**, default = ``50`` or ``length(zval)`` if present.)
    %
    %   Returns
    %   -------
    %
    %       bcrd
    %
    %           The output matrix of MATLAB doubles of shape ``[npnt, 2]``
    %           (or ``[npnt, 3]`` if the input argument ``zval`` is present)
    %           containing the coordinates of a set of ``npnt`` points on
    %           the boundary of the target 2D ellipsoid.
    %
    %   Interface
    %   ---------
    %
    %       pm.geom.ell2.getBorder()
    %       pm.geom.ell2.getBorder(gramian)
    %       pm.geom.ell2.getBorder(gramian, center)
    %       pm.geom.ell2.getBorder(gramian, center, zval)
    %       pm.geom.ell2.getBorder(gramian, center, zval, npnt)
    %
    %   Example
    %   -------
    %
    %       bcrd = pm.geom.ell2.getBorder();
    %       figure; h = plot(bcrd(:, 1), bcrd(:, 2), '-');
    %
    %       bcrd = pm.geom.ell2.getBorder([], [], 1);
    %       figure; h = plot3(bcrd(:, 1), bcrd(:, 2), bcrd(:, 3), '-');
    %
    %       npnt = 500;
    %       range = 10 * pi * [-1 : 2 / (npnt - 1) : 1];
    %       bcrd = pm.geom.ell2.getBorder([], [], sin(range / 5) + sin(range) / 5, npnt);
    %       figure; plot3(bcrd(:, 1), bcrd(:, 2), bcrd(:, 3));
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
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
        if  1 < length(zval)
            npnt = length(zval);
        else
            npnt = 50;
        end
    end
    if  isempty(center)
        center = zeros(2, 1);
    end
    if  isempty(gramian)
        gramian = eye(2);
    end
    if ~isempty(zval)
        bcrd = zeros(npnt, 3);
        bcrd(:, 3) = zval;
    else
        bcrd = zeros(npnt, 2);
    end
    independentVariable = linspace(0, 2 * pi, npnt)';
    xval = cos(independentVariable);
    yval = sin(independentVariable);
    ap = [xval(:) yval(:)]';
    [eigvec, eigval] = eig(gramian);
    eigval = sqrt(eigval); % convert variance to std.
    bcrd(:, 1 : 2) = transpose(eigvec * eigval * ap + repmat(center(:), 1, npnt));
end