function bcrd = getEllipsoidBoundary(covMat, meanVec, npoint) % returns the coordinates of the boundary of the ellipsoid
    if nargin<3 || isempty(npoint); npoint = 50; end
    if nargin<2 || isempty(meanVec); meanVec = zeros(length(covMat(:,1))); end
    independentVariable = linspace(0, 2*pi, npoint)';
    xval = cos(independentVariable);
    yval = sin(independentVariable);
    ap = [xval(:) yval(:)]';
    [eigenVectors,eigenValues] = eig(covMat);
    eigenValues = sqrt(eigenValues); % convert variance to std
    bcrd = transpose( eigenVectors * eigenValues * ap + repmat(meanVec(:), 1, npoint) );
    %h = plot(bcrd(:,1), bcrd(:,2), '-');
end
