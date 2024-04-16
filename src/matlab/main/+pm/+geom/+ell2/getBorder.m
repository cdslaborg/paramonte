function bcrd = getBorder(gramian, center, npnt)
    %
    %   Return a matrix of MATLAB doubles of shape ``[npnt, 2]``
    %   containing the coordinates of a set of points on the boundary
    %   of a 2D ellipsoid whose Gramian matrix is specified as input
    %   and whose center is also optionally specified.
    %
    %   Parameters
    %   ----------
    %
    %       gramian
    %
    %           The input square matrix of MATLAB doubles of shape ``[2, 2]``
    %           containing the Gramian of the target 2D ellipsoid whose
    %           boundary points are to be returned.
    %           (**optional**, default = ``eye(2, 2)``)
    %
    %       center
    %
    %           The input vector of MATLAB doubles of size ``2`` containing
    %           the 2D coordinates of the center of the target 2D ellipsoid
    %           whose boundary points are to be returned.
    %           (**optional**, default = ``zeros(2, 1)``)
    %
    %       npnt
    %
    %           The input scalar MATLAB whole number containing the number of
    %           points to return on the boundary of the target 2D ellipsoid.
    %           (**optional**, default = ``50``)
    %
    %   Returns
    %   -------
    %
    %       bcrd
    %
    %           The output matrix of MATLAB doubles of shape ``[npnt, 2]``
    %           containing the set of ``npnt`` points on the boundary of
    %           the target 2D ellipsoid.
    %
    %   Interface
    %   ---------
    %
    %       pm.geom.ell2.getBorder()
    %       pm.geom.ell2.getBorder(gramian)
    %       pm.geom.ell2.getBorder(gramian, center)
    %       pm.geom.ell2.getBorder(gramian, center, npnt)
    %       pm.geom.ell2.getBorder(gramian, center, npnt)
    %
    %       pm.geom.ell2.getBorder([])
    %       pm.geom.ell2.getBorder([], center)
    %       pm.geom.ell2.getBorder([], center, npnt)
    %       pm.geom.ell2.getBorder(gramian, center, [])
    %       pm.geom.ell2.getBorder(gramian, [], npnt)
    %       pm.geom.ell2.getBorder(gramian, [])
    %
    %   Example
    %   -------
    %
    %       bcrd = pm.geom.ell2.getBorder();
    %       h = plot(bcrd(:, 1), bcrd(:, 2), '-');
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
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
        gramian = eye(2, 2);
    end
    independentVariable = linspace(0, 2 * pi, npnt)';
    xval = cos(independentVariable);
    yval = sin(independentVariable);
    ap = [xval(:) yval(:)]';
    [eigvec, eigval] = eig(gramian);
    eigval = sqrt(eigval); % convert variance to std.
    bcrd = transpose(eigvec * eigval * ap + repmat(center(:), 1, npnt));
end