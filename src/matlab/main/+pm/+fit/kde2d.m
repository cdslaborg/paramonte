%>  \brief
%>  Compute the Gaussian kernel density estimate of the input bivariate data.
%>
%>  \details
%>  This function offers a fast and accurate bivariate
%>  kernel density estimator with diagonal bandwidth matrix.<br>
%>  The kernel is assumed to be Gaussian.<br>
%>  The two bandwidth parameters are chosen optimally without ever
%>  using/assuming a parametric model for the data or any *rules of thumb*.<br>
%>  Particularly, this algorithm is immune to accuracy failures in the estimation
%>  of multimodal densities with widely separated modes (see examples below).<br>
%>
%>  \param[in]  data    :   The input matrix of shape ``[nobs, 2]`` containing the
%>                          continuous data whose density estimate is to be returned.<br>
%>  \param[in]  resol   :   The input scalar whole number representing the size of
%>                          the ``resol`` by ``resol`` grid over which the density
%>                          is computed where ``resol`` has to be a power of ``2``,
%>                          otherwise ``resol = 2 ^ ceil(log2(resol))`` will be used.<br>
%>                          (**optional**. If missing or empty, it will be set to ``2^8``.)
%>  \param[in]  xymin   :   The input vector of size ``2`` containing the lower limits
%>                          of the bounding box over which the density is computed.<br>
%>                          (**optional**. If missing or empty, it will be set to an appropriate value as prescribed in the note below.)
%>  \param[in]  xymax   :   The input vector of size ``2`` containing the upper limits
%>                          of the bounding box over which the density is computed.<br>
%>                          (**optional**. If missing or empty, it will be set to an appropriate value as prescribed in the note below.)
%>
%>  \return
%>  ``bandwidth``       :   The output row vector with the two optimal bandwidths for a bivaroate Gaussian kernel.<br>
%>                          The format of the vector is ``bandwidth = [bandwidth_x, bandwidth_y]``.<br>
%>                          <br>
%>  ``density``         :   The output ``resol`` by ``resol`` matrix containing
%>                          the density values over the ``resol`` by ``resol`` grid.<br>
%>                          This argument is not computed unless the output argument is present.<br>
%>                          <br>
%>  ``meshx``           :   The x-axis meshgrid over which the variable "density" has been computed.<br>
%>                          This output is particularly useful in combination with MATLAB intrinsic ``surf``.<br>
%>                          The intended usage is ``surf(meshx, meshy, density)``.<br>
%>                          <br>
%>  ``meshy``           :   The y-axis meshgrid over which the variable "density" has been computed.<br>
%>                          This output is particularly useful in combination with MATLAB intrinsic ``surf``.<br>
%>                          The intended usage is ``surf(meshx, meshy, density)``.<br>
%>
%>  \interface{kde2d}
%>  \code{.m}
%>
%>      [bandwidth, density, meshx, meshy] = pm.fit.kde2d(data);
%>      [bandwidth, density, meshx, meshy] = pm.fit.kde2d(data, resol);
%>      [bandwidth, density, meshx, meshy] = pm.fit.kde2d(data, resol, xymin);
%>      [bandwidth, density, meshx, meshy] = pm.fit.kde2d(data, resol, xymin, xymax);
%>
%>  \endcode
%>
%>  \note
%>  The procedure for computing the default values for the
%>  input arguments ``xymin`` and ``xymax`` are as follows:<br>
%>  \code{.m}
%>      dmax = max(data, [], 1);
%>      dmin = min(data, [], 1);
%>      dlim = dmax - dmin;
%>      xymax = dmax + dlim / 4;
%>      xymin = dmin - dlim / 4;
%>  \endcode
%>
%>  \see
%>  Kernel density estimation via diffusion, Z. I. Botev, J. F. Grotowski,
%>  and D. P. Kroese (2010) Annals of Statistics, Volume 38, Number 5, pages 2916-2957.<br>
%>  [pm.vis.subplot.Contour](@ref Contour)<br>
%>  [pm.vis.cascade.Contour](@ref Contour)<br>
%>  [pm.vis.plot.Contour](@ref Contour)<br>
%>  [pm.vis.tile.Contour](@ref Contour)<br>
%>  [pm.vis.subplot.Contourf](@ref Contourf)<br>
%>  [pm.vis.cascade.Contourf](@ref Contourf)<br>
%>  [pm.vis.plot.Contourf](@ref Contourf)<br>
%>  [pm.vis.tile.Contourf](@ref Contourf)<br>
%>  [pm.vis.subplot.Contour3](@ref Contour3)<br>
%>  [pm.vis.cascade.Contour3](@ref Contour3)<br>
%>  [pm.vis.plot.Contour3](@ref Contour3)<br>
%>  [pm.vis.tile.Contour3](@ref Contour3)<br>
%>
%>  \example{mvnmix}
%>  \include{lineno} example/fit/kde2d/mvnmix/main.m
%>  \vis{mvnmix}
%>  \image html example/fit/kde2d/mvnmix/pm.fit.kde2d.mvnmix.png width=700
%>
%>  \example{mvnmixfar}
%>  \include{lineno} example/fit/kde2d/mvnmixfar/main.m
%>  \vis{mvnmixfar}
%>  \image html example/fit/kde2d/mvnmixfar/pm.fit.kde2d.mvnmixfar.png width=700
%>
%>  \example{sincos}
%>  \include{lineno} example/fit/kde2d/sincos/main.m
%>  \vis{sincos}
%>  \image html example/fit/kde2d/sincos/pm.fit.kde2d.sincos.png width=700
%>
%>  \final{kde2d}
%>
%>  Copyright (c) 2015, Dr. Zdravko Botev
%>  All rights reserved.
%>
%>  Redistribution and use in source and binary forms, with or without
%>  modification, are permitted provided that the following conditions are
%>  met:
%>
%>  * Redistributions of source code must retain the above copyright
%>    notice, this list of conditions and the following disclaimer.
%>  * Redistributions in binary form must reproduce the above copyright
%>    notice, this list of conditions and the following disclaimer in
%>    the documentation and/or other materials provided with the distribution
%>  * Neither the name of the The University of New South Wales nor the names
%>    of its contributors may be used to endorse or promote products derived
%>    from this software without specific prior written permission.
%>
%>  This software is provided by the copyright holders and contributors "as is"
%>  and any express or implied warranties, including, but not limited to, the
%>  implied warranties of merchantability and fitness for a particular purpose
%>  are disclaimed. In no event shall the copyright owner or contributors be
%>  liable for any direct, indirect, incidental, special, exemplary, or
%>  consequential damages (including, but not limited to, procurement of
%>  substitute goods or services; loss of use, data, or profits; or business
%>  interruption) however caused and on any theory of liability, whether in
%>  contract, strict liability, or tort (including negligence or otherwise)
%>  arising in any way out of the use of this software, even if advised of the
%>  possibility of such damage.
%>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function [bandwidth, density, meshx, meshy] = kde2d(data, resol, xymin, xymax)

    %global nobs asq rangesq

    if  nargin < 2
        resol = 2 ^ 8;
    end

    %%%% round up resol to the next power of 2;

    resol = 2 ^ ceil(log2(resol));

    nobs = size(data, 1);
    if  nargin < 3
        dmax = max(data, [], 1);
        dmin = min(data, [], 1);
        dlim = dmax - dmin;
        xymax = dmax + dlim / 2;
        xymin = dmin - dlim / 2;
    end
    scaling = xymax - xymin;
    if  nobs <= size(data, 2)
        error('data has to be a ``nobs`` by ``2`` array where each row represents a two dimensional observation.');
    end
    transformed_data = (data - repmat(xymin, nobs, 1)) ./ repmat(scaling, nobs, 1);

    %%%% Bin the data uniformly using regular grid.

    initial_data = ndhist(transformed_data, resol);

    %%%% Discrete cosine transform of initial data.

    a =  dct2d(initial_data);

    %%%% Now compute the optimal bandwidth^2.

    asq = a .^ 2;
    rangesq = (0 : resol - 1) .^ 2;
    t_star = getRoot(@(t)(t - evolve(t, nobs, asq, rangesq)), nobs);
    p_02 = func([0,2], t_star, nobs, asq, rangesq);
    p_20 = func([2,0], t_star, nobs, asq, rangesq);
    p_11 = func([1,1], t_star, nobs, asq, rangesq);
    t_y = (p_02 ^ (3 / 4) / (4 * pi * nobs * p_20 ^ (3 / 4) * (p_11 + sqrt(p_20 * p_02)))) ^ (1 / 3);
    t_x = (p_20 ^ (3 / 4) / (4 * pi * nobs * p_02 ^ (3 / 4) * (p_11 + sqrt(p_20 * p_02)))) ^ (1 / 3);

    %%%% Smooth the discrete cosine transform of initial data using t_star.

    a_t = exp(-(0 : resol - 1)' .^ 2 * pi ^ 2 * t_x / 2) * exp(-(0 : resol - 1) .^ 2 * pi ^ 2 * t_y / 2) .* a;

    %%%% Now apply the inverse discrete cosine transform.

    if  nargout > 1
        density = idct2d(a_t) * (numel(a_t) / prod(scaling));
        %%%% remove any negative density values
        density(density < 0) = eps;
        [meshx,meshy] = meshgrid(xymin(1) : scaling(1) / (resol - 1) : xymax(1), xymin(2) : scaling(2) / (resol - 1) : xymax(2));
    end

    bandwidth = sqrt([t_x, t_y]) .* scaling;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>  \cond excluded
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out, time] = evolve(t, nobs, asq, rangesq)
    sumFunc = func([0, 2], t, nobs, asq, rangesq) ...
            + func([2, 0], t, nobs, asq, rangesq) ...
            + func([1, 1], t, nobs, asq, rangesq) * 2;
    time = (2 * pi * nobs * sumFunc) ^ (-1 / 3);
    out = (t - time) / time;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = func(s, t, nobs, asq, rangesq)
    if sum(s) <= 4
        sumFunc = func([s(1) + 1, s(2)], t, nobs, asq, rangesq) ...
                + func([s(1), s(2) + 1], t, nobs, asq, rangesq);
        const = (1 + 1 / 2 ^ (sum(s) + 1)) / 3;
        time = (-2 * const * K(s(1)) * K(s(2)) / nobs / sumFunc) ^ (1 / (2 + sum(s)));
        out = psi(s, time, asq, rangesq);
    else
        out = psi(s, t, asq, rangesq);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = psi(s, time, asq, rangesq)
    %global asq rangesq
    %%%% The input ``s`` is a vector.
    w = exp(-rangesq * pi ^ 2 * time) .* [1, .5 * ones(1, length(rangesq) - 1)];
    wx = w .* (rangesq .^ s(1));
    wy = w .* (rangesq .^ s(2));
    out = (-1) ^ sum(s) * (wy * asq * wx') * pi ^ (2 * sum(s));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = K(s)
    out = (-1) ^ s * prod((1 : 2 : 2 * s - 1)) / sqrt(2 * pi);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = dct2d(data)
    %%%%
    %%%% Compute the two-dimensional discrete cosine transform of data data is an ``ndim`` cube.
    %%%%
    [nrow, ncol] =  size(data);
    if  nrow ~= ncol
        error('data is not a square array!')
    end
    %%%% Compute weights to multiply DFT coefficients.
    w = [1; 2 * (exp(-i * (1 : nrow - 1) * pi / (2 * nrow))).'];
    weight = w(:, ones(1, ncol));
    data = dct1d(dct1d(data)')';
    function transform1d = dct1d(x)
        %%%% Re-order the elements of the columns of x.
        x = [x(1 : 2 : end, :); x(end : -2 : 2, :)];
        %%%% Multiply FFT by weights.
        transform1d = real(weight .* fft(x));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = idct2d(data)
    %%%% Computes the 2 dimensional inverse discrete cosine transform.
    [nrow, ncol] = size(data);
    % Compute weights.
    w = exp(i * (0 : nrow - 1) * pi / (2 * nrow)).';
    weights = w(:, ones(1, ncol));
    data = idct1d(idct1d(data)');
    function out = idct1d(x)
        y = real(ifft(weights .* x));
        out = zeros(nrow, ncol);
        out(1 : 2 : nrow, :) = y(1 : nrow / 2, :);
        out(2 : 2 : nrow, :) = y(nrow : -1 : nrow / 2 + 1, :);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function binned_data = ndhist(data, nbin)
    %%%%    This function computes the histogram of an ndim-dimensional data set.
    %%%%    'data' is ``ndat`` rows by ``ndim`` columns.
    %%%%    The input ``nbin`` is the number of bins used in each dimension.
    %%%%    The output ``binned_data`` is a hypercube of length ``nbin``.
    [ndat, ncol] = size(data);
    bins = zeros(ndat, ncol);
    for icol = 1 : ncol
        [dum, bins(:, icol)] = histc(data(:, icol), [0 : 1 / nbin : 1], 1);
        bins(:, icol) = min(bins(:, icol), nbin);
    end
    %%%% Combine the  vectors of 1D bin counts into a grid of ND bin counts.
    binned_data = accumarray(bins(all(bins > 0, 2), :), 1 / ndat, nbin(ones(1, ncol)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function root = getRoot(getFunc, nobs)
    %%%% Find smallest root whenever there is more than one.
    nobs = 50 * (nobs <= 50) + 1050 * (nobs >= 1050) + nobs * ((nobs < 1050) & (nobs > 50));
    tol = 10^ -12 + 0.01 * (nobs - 50) / 1000;
    flag = 0;
    while flag == 0
        try
            root = fzero(getFunc, [0, tol]);
            flag = 1;
        catch
            %%%% double search interval.
            tol = min(tol * 2, .1);
        end
        %%%% if all else fails:
        if tol == .1
            root = fminbnd(@(x) abs(getFunc(x)), 0, .1);
            flag = 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>  \endcond
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
