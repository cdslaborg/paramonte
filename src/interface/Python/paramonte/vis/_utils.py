####################################################################################################################################
####################################################################################################################################
####
####   ParaMonte: plain powerful parallel Monte Carlo library.
####
####   Copyright (C) 2012-present, The Computational Data Science Lab
####
####   This file is part of the ParaMonte library.
####
####   ParaMonte is free software: you can redistribute it and/or modify it
####   under the terms of the GNU Lesser General Public License as published
####   by the Free Software Foundation, version 3 of the License.
####
####   ParaMonte is distributed in the hope that it will be useful,
####   but WITHOUT ANY WARRANTY; without even the implied warranty of
####   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
####   GNU Lesser General Public License for more details.
####
####   You should have received a copy of the GNU Lesser General Public License
####   along with the ParaMonte library. If not, see,
####
####       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
####
####   ACKNOWLEDGMENT
####
####   As per the ParaMonte library license agreement terms,
####   if you use any parts of this library for any purposes,
####   we ask you to acknowledge the use of the ParaMonte library
####   in your work (education/research/industry/development/...)
####   by citing the ParaMonte library as described on this page:
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
