####################################################################################################################################
####################################################################################################################################
####
####   MIT License
####
####   ParaMonte: plain powerful parallel Monte Carlo library.
####
####   Copyright (C) 2012-present, The Computational Data Science Lab
####
####   This file is part of the ParaMonte library.
####
####   Permission is hereby granted, free of charge, to any person obtaining a 
####   copy of this software and associated documentation files (the "Software"), 
####   to deal in the Software without restriction, including without limitation 
####   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
####   and/or sell copies of the Software, and to permit persons to whom the 
####   Software is furnished to do so, subject to the following conditions:
####
####   The above copyright notice and this permission notice shall be 
####   included in all copies or substantial portions of the Software.
####
####   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
####   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
####   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
####   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
####   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
####   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
####   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
####
####   ACKNOWLEDGMENT
####
####   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
####   As per the ParaMonte library license agreement terms, if you use any parts of 
####   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
####   work (education/research/industry/development/...) by citing the ParaMonte 
####   library as described on this page:
####
####       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
####
####################################################################################################################################
####################################################################################################################################

import numpy as np
from scipy.fft import dctn, idctn
from scipy.optimize import brentq

####################################################################################################################################
#### kde2d
####################################################################################################################################

def kde2d(x, y, n=256, limits=None):
    """

    Return the 2d density map from discrete observations via 
    2-dimensional diffusion Kernel density estimation.

    First the input data is binned. After binning, the function 
    determines the optimal bandwidth according to the diffusion-based 
    method. It then smooths the binned data over the grid using a 
    Gaussian kernel with a standard deviation corresponding to 
    that bandwidth.

        This module is based on the KDE estimation method of 

            Z. I. Botev, J. F. Grotowski, D. P. Kroese: 
            Kernel density estimation via diffusion. 
            Annals of Statistics 38 (2010), no. 5, 2916--2957. 
            doi:10.1214/10-AOS799

        and 

            John Hennig 
            DOI: 10.5281/zenodo.3830437 

        **Parameters**

            x

                A lists of array of numbers that represent discrete 
                observations of a random variable with two coordinate 
                components. The observations are binned on a grid of 
                n*n points, where ``n`` must be a power of 2 or will 
                be coerced to the next one. If ``x`` and ``y`` are not 
                the same length, the algorithm will raise a ``ValueError``.

            y

                A lists of array of numbers that represent discrete 
                observations of a random variable with two coordinate 
                components. The observations are binned on a grid of 
                n*n points, where ``n`` must be a power of 2 or will 
                be coerced to the next one. If ``x`` and ``y`` are not 
                the same length, the algorithm will raise a ``ValueError``.

            n (optional)

                The number of grid points. It must be a power of 2. 
                Otherwise, it will be coerced to the next power of two.
                The default is 256.

            limits (optional)

                Data ``limits`` specified as a tuple of tuples denoting
                ``((xmin, xmax), (ymin, ymax))``. If any of the values 
                are ``None``, they will be inferred from the data. Each 
                tuple, or even both of them, may also be replaced by a 
                single value denoting the upper bound of a range 
                centered at zero. The default is ``None``.

        **Returns**

            A tuple whose elements are the following:

                density

                    The density map of the data.

                grid

                    The grid at which the density is computed.

                bandwidth

                    The optimal values (per axis) that the algorithm has 
                    determined. If the algorithm does not converge, it 
                    will raise a ``ValueError``.

    ---------------------------------------------------------------------------
    """

    # Convert to arrays in case lists are passed in.

    x = np.array(x)
    y = np.array(y)

    # Make sure numbers of data points are consistent.

    N = len(x)
    if len(y) != N:
        raise ValueError('x and y must have the same length.')

    # Round up number of bins to next power of two.

    n = int(2**np.ceil(np.log2(n)))

    # Determine missing data limits.

    if limits is None:

        xmin = xmax = ymin = ymax = None

    elif isinstance(limits, tuple):

        (xlimits, ylimits) = limits
        if xlimits is None:
            xmin = xmax = None
        elif isinstance(xlimits, tuple):
            (xmin, xmax) = xlimits
        else:
            xmin = -xlimits
            xmax = +xlimits
        if ylimits is None:
            ymin = ymax = None
        elif isinstance(ylimits, tuple):
            (ymin, ymax) = ylimits
        else:
            ymin = -ylimits
            ymax = +ylimits
    else:

        xmin = -limits
        xmax = +limits
        ymin = -limits
        ymax = +limits

    if None in (xmin, xmax):
        delta = x.max() - x.min()
        if xmin is None:
            xmin = x.min() - delta/4
        if xmax is None:
            xmax = x.max() + delta/4

    if None in (ymin, ymax):
        delta = y.max() - y.min()
        if ymin is None:
            ymin = y.min() - delta/4
        if ymax is None:
            ymax = y.max() + delta/4

    deltax = xmax - xmin
    deltay = ymax - ymin

    # Bin samples on regular grid.

    (binned, xedges, yedges) = np.histogram2d(x, y, bins=n, range=((xmin, xmax), (ymin, ymax)))
    grid = (xedges[:-1], yedges[:-1])

    # Compute discrete cosine transform. Adjust first component.

    transformed = dctn(binned/N)
    transformed[0, :] /= 2
    transformed[:, 0] /= 2

    # Pre-compute squared indices and transform components before solver loop.

    k  = np.arange(n, dtype="float") # float avoids integer overflow.
    k2 = k**2
    a2 = transformed**2

    # Define internal functions to be solved iteratively.

    def γ(t):
        sigma = psi(0, 2, t) + psi(2, 0, t) + 2*psi(1, 1, t)
        γ = (2*np.pi*N*sigma)**(-1/3)
        return (t - γ) / γ

    def psi(i, j, t):
        if i + j <= 4:
            sigma  = abs(psi(i+1, j, t) + psi(i, j+1, t))
            C  = (1 + 1/2**(i+j+1)) / 3
            pii = np.product(np.arange(1, 2*i, 2))
            pij = np.product(np.arange(1, 2*j, 2))
            t  = (C*pii*pij / (np.pi*N*sigma)) ** (1/(2+i+j))
        w = 0.5 * np.ones(n)
        w[0] = 1
        w = w * np.exp(-np.pi**2 * k2*t)
        wx = w * k2**i
        wy = w * k2**j
        return (-1)**(i+j) * np.pi**(2*(i+j)) * wy @ a2 @ wx

    # Solve for optimal diffusion time t*.

    try:
        ts = brentq(lambda t: t - γ(t), 0, 0.1)
    except ValueError:
        raise ValueError('Bandwidth optimization did not converge.') from None

    # Calculate diffusion times along x- and y-axis.

    psi02 = psi(0, 2, ts)
    psi20 = psi(2, 0, ts)
    psi11 = psi(1, 1, ts)
    tx1 = (psi02**(3/4) / (4*np.pi*N*psi20**(3/4) * (psi11 + np.sqrt(psi02*psi20))) )**(1/3)
    tx2 = (psi20**(3/4) / (4*np.pi*N*psi02**(3/4) * (psi11 + np.sqrt(psi02*psi20))) )**(1/3)

    # Note:
    # The above uses the nomenclature from the paper. In the Matlab
    # reference, tx1 is called t_y, while tx2 is t_x. This is a curious
    # change in notation. It may be related to the fact that image
    # coordinates are typically in (y,x) index order, whereas matrices,
    # such as the binned histogram (in Matlab as much as in Python),
    # are in (x,y) order. The Matlab code eventually does return
    # image-like index order, though it never explicitly transposes
    # the density matrix. That is implicitly handled by its custom
    # implementation of the inverse transformation (idct2d), which
    # only employs one matrix transposition, not two as its forward
    # counterpart (dct2d).

    # Apply Gaussian filter with optimized kernel.

    smoothed = transformed * np.outer(np.exp(-np.pi**2 * k2 * tx2/2), np.exp(-np.pi**2 * k2 * tx1/2))

    # Reverse transformation.

    smoothed[0, :] *= 2
    smoothed[:, 0] *= 2
    inverse = idctn(smoothed)

    # Normalize density.

    density = np.transpose(inverse) * n/deltax * n/deltay

    # Determine bandwidth from diffusion times.

    bandwidth = np.array([np.sqrt(tx2)*deltax, np.sqrt(tx1)*deltay])

    # Return results.

    return (density, grid, bandwidth)
