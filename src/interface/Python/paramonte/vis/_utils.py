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

####################################################################################################################################
#### colorbar
####################################################################################################################################

def isColorBar(ax):
    """

    Attempt to guess whether an input matplotlib
    ``Axes`` is home to a ``colorbar`` object.

        **Parameters**

            ax

                An Axes instance.

        **Returns**

            A boolean value of ``True`` if the x xor y axis
            satisfies all of the following and thus looks
            like a ``colorbar`` object:
                No ticks, no tick labels, no axis label

    """

    xcb = (len(ax.get_xticks()) == 0) and (len(ax.get_xticklabels()) == 0) and (len(ax.get_xlabel()) == 0) #and (ax.get_xlim() == (0, 1))
    ycb = (len(ax.get_yticks()) == 0) and (len(ax.get_yticklabels()) == 0) and (len(ax.get_ylabel()) == 0) #and (ax.get_ylim() == (0, 1))
    return xcb != ycb

####################################################################################################################################
#### isPresentColorBar
####################################################################################################################################

def isPresentColorBar():
    from matplotlib import pyplot as plt
    return any([isColorBar(ax) for ax in np.atleast_1d(plt.gcf().axes).flatten()])
