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
import typing as tp
import pandas as pd
import weakref as wref
import _pmutils as pmutils
import _paramonte as pm
from paramonte.vis.Target import Target
from paramonte.vis._BasePlot import BasePlot

Struct = pmutils.Struct
newline = pmutils.newline

####################################################################################################################################
#### DensityPlot class
####################################################################################################################################

class DensityPlot(BasePlot):
    """

    This is the DensityPlot class for generating instances
    of histogram and contour figures in two and three dimensions
    based on a wide range of plotting tools from the matplotlib,
    seaborn, and other Python libraries.

    Normally, the public users are not supposed to use this
    class directly, although they can for the purposes other
    than plotting the ParaMonte simulation files.

        **Parameters**

            plotType

                A string indicating the name of the plot to be constructed.

            dataFrame (optional)

                A pandas dataFrame whose data will be plotted.

            methodName (optional)

                The name of the ParaMonte sample requesting the BasePlot.

            reportEnabled (optional)

                A boolean whose value indicates whether guidelines should be
                printed in the standard output.

            resetPlot (optional)

                A function that resets the properties of the plot as desired
                from outside. If provided, a pointer to this function will be
                saved for future internal usage.

        **Attributes**

            xcolumns

                An attribute that determines the columns of dataFrame
                to be visualized as the X-axis. It can have three forms:

                    1.  A list of column indices in dataFrame.
                    2.  A list of column names in dataFrame.columns.
                    3.  A ``range(start,stop,step)`` of column indices.

                Examples:

                    1.  ``xcolumns = [0,1,4,3]``
                    2.  ``xcolumns = ["SampleLogFunc","SampleVariable1"]``
                    3.  ``xcolumns = range(17,7,-2)``

                The default behavior includes all columns of the dataFrame.

            ycolumns (only in kdeplot2, contour, contourf, contour3 plots)

                An attribute that determines the columns of dataFrame
                to be visualized as the Y-axis. It can have three forms:

                    1.  A list of column indices in dataFrame.
                    2.  A list of column names in dataFrame.columns.
                    3.  A ``range(start,stop,step)`` of column indices.

                Examples:

                    1.  ``ycolumns = [0,1,4,3]``
                    2.  ``ycolumns = ["SampleLogFunc","SampleVariable1"]``
                    3.  ``ycolumns = range(17,7,-2)``

                The default behavior includes all columns of the dataFrame.

            rows

                An attribute that determines the rows of dataFrame to be
                visualized. It can be either:

                    1.  A ``range(start,stop,step)``, or,
                    2.  A list of row indices in dataFrame.index.

                Examples:

                    1.  ``rows = range(17,7,-2)``
                    2.  ``rows = [i for i in range(7,17)]``

                The default behavior includes all rows of the dataFrame.

            histplot (available only in histplot objects)

                A structure with two attributes:

                    enabled

                        A boolean indicating whether a call to the ``histplot()``
                        function of the seaborn library should be made or not.

                    kws

                        A structure whose components are directly passed as
                        keyword arguments to the ``histplot()`` function.

                        Example usage:

                            .. code-block:: python
                                :linenos:

                                histplot.kws.stat = "count"

                        **NOTE**

                        If a desired property is missing among the ``kws``
                        attributes, simply add the field and its value to
                        the component.

            kdeplot (available only in ``kdeplot1`` and ``kdeplot2`` objects)

                A structure with two attributes:

                    enabled

                        A boolean indicating whether a call to the ``kdeplot()``
                        function of the seaborn library should be made or not.

                    kws

                        A structure whose components are directly passed as
                        keyword arguments to the ``kdeplot()`` function.

                        Example usage:

                            .. code-block:: python
                                :linenos:

                                kdeplot.kws.cbar = True
                                kdeplot.kws.shade = False

                        **NOTE**

                        If a desired property is missing among the ``kws``
                        attributes, simply add the field and its value to
                        the component.

            contour (available only in ``contour`` and ``contour3`` objects)

                A structure with two attributes:

                    enabled

                        A boolean indicating whether a call to the ``contour()``
                        function of the seaborn library should be made or not.

                    kws

                        A structure whose components are directly passed as
                        keyword arguments to the ``contour()`` function.

                        Example usage:

                            .. code-block:: python
                                :linenos:

                                contour.kws.cmap = "winter"
                                contour.kws.levels = 50

                        **NOTE**

                        If a desired property is missing among the ``kws``
                        attributes, simply add the field and its value to
                        the component.

            contourf (available only in ``contourf`` objects)

                A structure with two attributes:

                    enabled

                        A boolean indicating whether a call to the ``contourf()``
                        function of the seaborn library should be made or not.

                    kws

                        A structure whose components are directly passed as
                        keyword arguments to the ``contourf()`` function.

                        Example usage:

                            .. code-block:: python
                                :linenos:

                                contour.kws.cmap = "winter"
                                contour.kws.levels = 50

                        **NOTE**

                        If a desired property is missing among the ``kws``
                        attributes, simply add the field and its value to
                        the component.

            gridSize (available in ``contour``, ``contourf``, ``contour3`` objects)

                An integer indicating the grid resolution for discretization of the
                data during the kernel density estimation. It must be a power of
                two, otherwise it will be changed to the next power of two at the
                time of using it. The default value is 512.

                Example usage:

                    .. code-block:: python

                        gridSize = 512

            limits (available in ``contour``, ``contourf``, ``contour3`` objects)

                Data ``limits`` used in the kernel density estimation
                specified as a tuple of tuples ``((xmin, xmax), (ymin, ymax))``.
                If any of the values are ``None``, they will be inferred
                from the data. Each tuple, or even both of them, may
                also be replaced by a single value denoting the
                upper bound of a range centered at zero.
                The default is ``None``.

                Example usage:

                    .. code-block:: python

                        limits = ( (-10, 10),(-5,5) )

            noiseDensity (available in ``contour``, ``contourf``, ``contour3``)

                A float indicating the threshold below which the kernel density
                estimate is considered to be noise and is rounded to zero.
                The higher this value is, the less noise will be visible.

                Example usage:

                    .. code-block:: python

                        noiseDensity = 1.e-5

            set

                A structure with two attributes:

                    enabled

                        A boolean indicating whether a call to the ``set()``
                        function of the seaborn library should be made or not.

                    kws

                        A structure whose components are directly passed as
                        keyword arguments to the ``set()`` function.

                        Example usage:

                            .. code-block:: python

                                set.kws.style = "darkgrid"

                        **NOTE**

                        If a desired property is missing among the ``kws``
                        attributes, simply add the field and its value to
                        the component.

            axes (available only in 1D and 2D plots)

                A structure with one attribute:

                    kws

                        A structure whose components are directly passed as
                        keyword arguments to the ``gca()`` function of the
                        matplotlib library.

                        Example usage:

                            .. code-block:: python

                                axes.kws.facecolor = "w"

                        **NOTE**

                        If a desired property is missing among the ``kws``
                        attributes, simply add the field and its value to
                        the component.

            axes3d (available only in 3D plots)

                A structure with one attribute:

                    kws

                        A structure whose components are directly passed as
                        keyword arguments to the ``Axes3D()`` function of the
                        matplotlib library.

                        Example usage:

                            .. code-block:: python

                                axes3d.kws.facecolor = "w"

                        **NOTE**

                        If a desired property is missing among the ``kws``
                        attributes, simply add the field and its value to
                        the component.

            figure

                A structure with two attributes:

                    enabled

                        A boolean indicating whether a call to the ``figure()``
                        function of the matplotlib library should be made or not.
                        If a call is made, a new figure will be generated.
                        Otherwise, the current active figure will be used.

                    kws

                        A structure whose components are directly passed as
                        keyword arguments to the ``figure()`` function.

                        Example usage:

                            .. code-block:: python

                                figure.kws.facecolor = "w"

                        **NOTE**

                        If a desired property is missing among the ``kws``
                        attributes, simply add the field and its value to
                        the component.

            colorbar (exists only for plots that require colorbar)

                A structure with two attributes:

                    enabled

                        A boolean indicating whether a call to the ``colorbar()``
                        function of the matplotlib library should be made or not.
                        If a call is made, a new figure will be generated.
                        Otherwise, the current active figure will be used.

                    kws

                        A structure whose components are directly passed as
                        keyword arguments to the ``colorbar()`` function of
                        the matplotlib library.

                        **NOTE**

                        If a desired property is missing among the ``kws``
                        attributes, simply add the field and its value to
                        the component.

                A colorbar will be added to a plot only if a color-mappings
                is requested in the plot.

            legend (may not exist for some types of plots)

                A structure with two attributes:

                    enabled

                        A boolean indicating whether a call to the ``legend()``
                        function of the matplotlib library should be made or not.
                        If a call is made, a new figure will be generated.
                        Otherwise, the current active figure will be used.

                    kws

                        A structure whose components are directly passed as
                        keyword arguments to the ``legend()`` function.

                        Example usage:

                        .. code-block:: python

                            legend.kws.labels = ["Variable1", "Variable2"]

                        **NOTE**

                        If a desired property is missing among the ``kws``
                        attributes, simply add the field and its value to
                        the component.

                A legend will be added to a plot only if no color-mappings are
                requested in the plot.

            currentFig

                A structure whose attributes are the outputs of various plotting
                tools used to make the current figure. These include the handle
                to the current figure, the handle to the current axes in the plot,
                the handle to the colorbar (if any exists), and other Python
                plotting tools used to make to generate the figure.

            target (available only in 1D and 2D plot objects)

                A callable object of the ParaMonte library's ``Target`` class
                which can be used to add target point or lines to the current
                active plot.

        **Returns**

                An object of class ``DensityPlot``.

    ---------------------------------------------------------------------------
    """

    ################################################################################################################################
    #### __init__
    ################################################################################################################################

    def __init__( self
                , plotType      # : str
                , dataFrame     = None # : tp.Optional[ pd.DataFrame ] = None
                , methodName    = "ParaMonte" # : tp.Optional[ str ] = "ParaMonte"
                , reportEnabled = True # : tp.Optional[ bool ] = True
                , resetPlot     = None
                ):

        super().__init__( plotType = plotType
                        , dataFrame = dataFrame
                        , methodName = methodName
                        , reportEnabled = reportEnabled
                        , resetPlot = resetPlot
                        )

        self._reset()
        if resetPlot is None: self._resetPlot = self._reset
        self._progress.note()

    ################################################################################################################################
    #### _reset
    ################################################################################################################################

    def _reset(self):

        super()._reset()
        self._indexOffset = 1

        if not self._type.is3d:
            self.target = Target()
            if self._type.is1d:
                self.target.axhline.enabled = False
                self.target.scatter.enabled = False

        self.xcolumns = None

        if not (self._type.isHistplot or self._type.isKdeplot1):

            self.ycolumns = None

        if self._type.isDiffusionPlot:

            self.colorbar = Struct()
            self.colorbar.kws = Struct()
            self.colorbar.enabled = True

        if self._type.isHistplot:

            self.histplot = Struct()
            self.histplot.kws = Struct()
            self.histplot.enabled = True

        if self._type.isKdeplot1 or self._type.isKdeplot2:

            self.kdeplot = Struct()
            self.kdeplot.enabled = True
            self.kdeplot.kws = Struct()

        if self._type.isJointplot:

            self.jointplot = Struct()
            self.jointplot.enabled = True
            self.jointplot.kws = Struct()

        if self._type.isJointplot:

            self.jointplot = Struct()
            self.jointplot.enabled = True
            self.jointplot.kws = Struct()

        if self._type.isContour or self._type.isContour3:

            self.contour = Struct()
            self.contour.enabled = True
            self.contour.kws = Struct()
            self.contour.kws.cmap = "winter"

        if self._type.isContourf:

            self.contourf = Struct()
            self.contourf.enabled = True
            self.contourf.kws = Struct()

        if self._type.isDiffusionPlot:
            self.noiseDensity = 1.e-3
            self.limits = None
            self.gridSize = 512


        self._isdryrun = True
        self.make()
        self._isdryrun = False

    ################################################################################################################################
    #### __call__
    ################################################################################################################################

    def __call__( self
                , reself : tp.Optional[ bool ] = False
                , **kwargs
                ):
        """

        Call the ``make()`` method of the current instance
        of the class.

            **Parameters**

                Any arguments that can be passed to the
                ``make()`` method of the plot object.

            **Returns**

                Any return value from the ``make()`` method
                of the plot object.

        """

        return self.make(reself, **kwargs)

    ################################################################################################################################
    #### make
    ################################################################################################################################

    def make( self
            , reself : tp.Optional[ bool ] = False
            , **kwargs
            ):
        """

        Generate a plot from the selected columns
        of the object's dataFrame.

            **Parameters**

                reself

                    A logical variable. If ``True``, an instance of
                    the object will be returned  to the calling routine
                    upon exit. The default value is ``False``.

            **Returns**

                The object self if ``reself = True`` otherwise, ``None``.
                However, this method causes side-effects by manipulating
                the existing attributes of the object.

                **NOTE**

                **This method causes side-effects by manipulating
                the existing attributes of the object.**

        """

        for key in kwargs.keys():
            if hasattr(self,key):
                setattr(self, key, kwargs[key])
            elif key=="dataFrame":
                setattr( self, "_dfref", wref.ref(kwargs[key]) )
            else:
                raise Exception ( "Unrecognized input '"+key+"' class attribute detected." + newline
                                + self._getDocString()
                                )

        # set what to plot

        cEnabled = not self._type.isHistplot

        ############################################################################################################################
        #### histplot properties
        ############################################################################################################################

        if self._type.isHistplot:
            if isinstance(self.histplot.kws,Struct):
                self.histplot.kws.legend = self.legend.enabled
                #if "legend" not in vars(self.histplot.kws).keys(): self.histplot.kws.legend = False
                if "kde" not in vars(self.histplot.kws).keys(): self.histplot.kws.kde = False
                if "bins" not in vars(self.histplot.kws).keys(): self.histplot.kws.bins = "auto"
                if "stat" not in vars(self.histplot.kws).keys(): self.histplot.kws.stat = "count"
                if "kde_kws" not in vars(self.histplot.kws).keys(): self.histplot.kws.kde_kws = dict()
                if "binwidth" not in vars(self.histplot.kws).keys(): self.histplot.kws.binwidth = None
                if "binrange" not in vars(self.histplot.kws).keys(): self.histplot.kws.binrange = None
                if "multiple" not in vars(self.histplot.kws).keys(): self.histplot.kws.multiple = "stack"
                if "element" not in vars(self.histplot.kws).keys(): self.histplot.kws.element = "step"
                if "shrink" not in vars(self.histplot.kws).keys(): self.histplot.kws.shrink = 1
                if "color" not in vars(self.histplot.kws).keys(): self.histplot.kws.color = None
                if "fill" not in vars(self.histplot.kws).keys(): self.histplot.kws.fill = True
                if "common_norm" not in vars(self.histplot.kws).keys(): self.histplot.kws.common_norm = True
                if "common_bins" not in vars(self.histplot.kws).keys(): self.histplot.kws.common_bins = False
                if "line_kws" not in vars(self.histplot.kws).keys(): self.histplot.kws.line_kws = dict()
                if "linewidth" not in self.histplot.kws.line_kws.keys(): self.histplot.kws.line_kws["linewidth"] = 0
                if "linestyle" not in self.histplot.kws.line_kws.keys(): self.histplot.kws.line_kws["linestyle"] = "-"
            else:
                raise Exception ( "The histplot.kws component of the current DensityPlot object must" + newline
                                + "be an object of class Struct(), essentially a structure with components" + newline
                                + "whose names are the input arguments to the histplot() function of the" + newline
                                + "seaborn library." + newline
                                + self._getDocString()
                                )

        ############################################################################################################################
        #### jointplot properties
        ############################################################################################################################

        if self._type.isJointplot:
            if isinstance(self.jointplot.kws,Struct):
                if "height" not in vars(self.jointplot.kws).keys(): self.jointplot.kws.height = 7
                if "space" not in vars(self.jointplot.kws).keys(): self.jointplot.kws.space = 0
                if "kind" not in vars(self.jointplot.kws).keys(): self.jointplot.kws.kind = "kde"
                if "fill" not in vars(self.jointplot.kws).keys(): self.jointplot.kws.fill = True
                if self.jointplot.kws.kind=="kde":
                    if "cmap" not in vars(self.jointplot.kws).keys() or self.jointplot.kws.cmap is None: self.jointplot.kws.cmap = "Blues"
                    if "cbar" not in vars(self.jointplot.kws).keys(): self.jointplot.kws.cbar = False
                    if "shade" not in vars(self.jointplot.kws).keys(): self.jointplot.kws.shade = True
                    if "n_levels" not in vars(self.jointplot.kws).keys(): self.jointplot.kws.n_levels = 100
                    if "thresh" not in vars(self.jointplot.kws).keys(): self.jointplot.kws.thresh = 0
                    #if "shade_lowest" not in vars(self.jointplot.kws).keys(): self.jointplot.kws.shade_lowest = True # deprecated in seaborn 0.11.0
                    if "cbar_kws" not in vars(self.jointplot.kws).keys(): self.jointplot.kws.cbar_kws = dict()
                    if "label" not in self.jointplot.kws.cbar_kws.keys(): self.jointplot.kws.cbar_kws["label"] = "Density"
            else:
                raise Exception ( "The jointplot.kws component of the current DensityPlot object must" + newline
                                + "be an object of class Struct(), essentially a structure with components" + newline
                                + "whose names are the input arguments to the jointplot() function of the" + newline
                                + "seaborn library." + newline
                                + self._getDocString()
                                )

        ############################################################################################################################
        #### kdeplot1 / kdeplot2 properties
        ############################################################################################################################

        if self._type.isKdeplot1 or self._type.isKdeplot2:
            if isinstance(self.kdeplot.kws,Struct):
                if "cbar" not in vars(self.kdeplot.kws).keys(): self.kdeplot.kws.cbar = True
                if "shade" not in vars(self.kdeplot.kws).keys(): self.kdeplot.kws.shade = True
                if "kernel" not in vars(self.kdeplot.kws).keys(): self.kdeplot.kws.kernel = "gau"
                if "zorder" not in vars(self.kdeplot.kws).keys(): self.kdeplot.kws.zorder = 0
            else:
                raise Exception ( "The kdeplot.kws component of the current DensityPlot object must" + newline
                                + "be an object of class Struct(), essentially a structure with components" + newline
                                + "whose names are the input arguments to the kdeplot() function of the" + newline
                                + "seaborn library." + newline
                                + self._getDocString()
                                )
            if self._type.isKdeplot2:
                if "cmap" not in vars(self.kdeplot.kws).keys() or self.kdeplot.kws.cmap is None: self.kdeplot.kws.cmap = "Blues"
                if "n_levels" not in vars(self.kdeplot.kws).keys(): self.kdeplot.kws.n_levels = 100
                if "thresh" not in vars(self.kdeplot.kws).keys(): self.kdeplot.kws.thresh = 0
                #if "shade_lowest" not in vars(self.kdeplot.kws).keys(): self.kdeplot.kws.shade_lowest = True # deprecated in seaborn 0.11.0
                if "cbar_kws" not in vars(self.kdeplot.kws).keys(): self.kdeplot.kws.cbar_kws = dict()
                if "label" not in self.kdeplot.kws.cbar_kws.keys(): self.kdeplot.kws.cbar_kws["label"] = "Density"

        ############################################################################################################################
        #### contour / contour3 properties
        ############################################################################################################################

        if self._type.isContour or self._type.isContour3:

            # if `cmap` is specified and is not None, it takes precedence over `colors` attribute.

            if isinstance(self.contour.kws,Struct):

                if "cmap" in vars(self.contour.kws).keys():
                    if self.contour.kws.cmap is None:
                        delattr(self.contour.kws, "cmap")
                    else:
                        if hasattr(self.contour.kws, "colors"): delattr(self.contour.kws, "colors")
                else:
                    if "colors" not in vars(self.contour.kws).keys() or self.contour.kws.colors is None: self.contour.kws.cmap = "Blues"

                if "alpha" not in vars(self.contour.kws).keys(): self.contour.kws.alpha = 1
                if "levels" not in vars(self.contour.kws).keys(): self.contour.kws.levels = 50
                if "linewidths" not in vars(self.contour.kws).keys(): self.contour.kws.linewidths = 1
                if "linestyles" not in vars(self.contour.kws).keys(): self.contour.kws.linestyles = "solid"

                if self._type.isContour3:
                    if "zdir" not in vars(self.contour.kws).keys(): self.contour.kws.zdir = "z"
                    if "extend3d" not in vars(self.contour.kws).keys(): self.contour.kws.extend3d = False

            else:

                raise Exception ( "The contour.kws component of the current DensityPlot object must" + newline
                                + "be an object of class Struct(), essentially a structure with components" + newline
                                + "whose names are the input arguments to the contour() function of the" + newline
                                + "matplotlib library." + newline
                                + self._getDocString()
                                )

        ############################################################################################################################
        #### contourf properties
        ############################################################################################################################

        if self._type.isContourf:
            if isinstance(self.contourf.kws,Struct):

                if "cmap" in vars(self.contourf.kws).keys():
                    if self.contourf.kws.cmap is None:
                        delattr(self.contourf.kws, "cmap")
                    else:
                        if hasattr(self.contourf.kws, "colors"): delattr(self.contourf.kws, "colors")
                else:
                    if "colors" not in vars(self.contourf.kws).keys() or self.contourf.kws.colors is None: self.contourf.kws.cmap = "Blues"

                if "alpha" not in vars(self.contourf.kws).keys(): self.contourf.kws.alpha = 1
                if "levels" not in vars(self.contourf.kws).keys(): self.contourf.kws.levels = 50

            else:

                raise Exception ( "The contourf.kws component of the current DensityPlot object must" + newline
                                + "be an object of class Struct(), essentially a structure with components" + newline
                                + "whose names are the input arguments to the contourf() function of the" + newline
                                + "matplotlib library." + newline
                                + self._getDocString()
                                )

        ############################################################################################################################
        #### legend properties
        ############################################################################################################################

        if not self._type.isDiffusionPlot and self.legend.enabled:
            if not isinstance(self.legend.kws,Struct):
                raise Exception ( "The legend.kws component of the current DensityPlot object must" + newline
                                + "be an object of class Struct(), essentially a structure with components" + newline
                                + "whose names are the input arguments to the legend() function of the" + newline
                                + "matplotlib library." + newline
                                + self._getDocString()
                                )

        ############################################################################################################################
        #### figure properties
        ############################################################################################################################

        if self.figure.enabled:
            if isinstance(self.figure.kws, Struct):
                if "dpi" not in vars(self.figure.kws).keys(): self.figure.kws.dpi = 150
                if "facecolor" not in vars(self.figure.kws).keys(): self.figure.kws.facecolor = "w"
                if "edgecolor" not in vars(self.figure.kws).keys(): self.figure.kws.edgecolor = "w"
            else:
                raise Exception ( "The figure.kws component of the current DensityPlot object must" + newline
                                + "be an object of class Struct(), essentially a structure with components" + newline
                                + "whose names are the input arguments to the figure() function of the" + newline
                                + "matplotlib library." + newline
                                + self._getDocString()
                                )

        ############################################################################################################################
        ############################################################################################################################
        if self._isdryrun: return
        ############################################################################################################################
        ############################################################################################################################

        import seaborn as sns
        import matplotlib.pyplot as plt
        if self._type.isDiffusionPlot: from paramonte.vis.kde2d import kde2d

        ############################################################################################################################
        #### generate figure and axes if needed
        ############################################################################################################################

        self._constructBasePlot()

        ############################################################################################################################
        #### check data type
        ############################################################################################################################

        self._checkDataType()

        ############################################################################################################################
        #### check rows presence. This must be checked here, because it depends on the integrity of the in input dataFrame.
        ############################################################################################################################

        if self.rows is None: self.rows = range(len(self._dfref().index))

        ############################################################################################################################
        #### check columns presence. This must be checked here, because it depends on the integrity of the in input dataFrame.
        ############################################################################################################################

        #### assign x columns to plot

        lgxicol = 0 # legend column indicator index. By default, there is only one column, hence index 0.
        if self.xcolumns is None:
            xcolindex = []
            xcolnames = ["Count"]
        else:
            try:
                xcolnames, xcolindex = pm.dfutils.getColNamesIndex(self._dfref().columns,self.xcolumns)
            except:
                raise Exception ( "The xcolumns component of the current DensityPlot object must" + newline
                                + "point to the name of a column of the input dataFrame to the DensityPlot" + newline
                                + "class constructor." + newline
                                + self._getDocString()
                                )

        #### assign x values if it is a single column of data

        xcolindexlen = len(xcolindex)
        if xcolindexlen==0:
            self._xvalues = np.array( self._dfref().index[self.rows] + self._indexOffset ).flatten()
        elif xcolindexlen==1:
            self._xvalues = self._dfref().iloc[self.rows,xcolindex[0]].values.flatten()
        elif self._type.isKdeplot2 or self._type.isJointplot:
            raise Exception ( "The length of xcolumns cannot be larger than 1 in jointplot or kdeplot2 objects." + newline
                            + self._getDocString()
                            )

        #### assign y properties for 2D plots

        if not (self._type.isHistplot or self._type.isKdeplot1):

            # assign y columns to plot

            lgyicol = 0 # legend column indicator index. By default, there is only one column, hence index 0.
            if self.ycolumns is None:
                ycolindex = []
                ycolnames = ["Count"]
            else:
                try:
                    ycolnames, ycolindex = pm.dfutils.getColNamesIndex(self._dfref().columns,self.ycolumns)
                except:
                    raise Exception ( "The ycolumns component of the current DensityPlot object must" + newline
                                    + "point to the name of a column of the input dataFrame to the DensityPlot" + newline
                                    + "class constructor." + newline
                                    + self._getDocString()
                                    )

            # check the consistency of the lengths

            ycolindexlen = len(ycolindex)
            maxLenColumns = np.max  (   [ xcolindexlen
                                        , ycolindexlen
                                        ]
                                    )

            if xcolindexlen!=maxLenColumns and xcolindexlen>1:
                raise Exception ( "length of xcolumns must be either 1 or equal to the lengths of ycolumns." + newline
                                + self._getDocString()
                                )
            if ycolindexlen!=maxLenColumns and ycolindexlen>1:
                raise Exception ( "length of ycolumns must be either 1 or equal to the lengths of xcolumns." + newline
                                + self._getDocString()
                                )

            #if self._type.is3d and isinstance(self.zlabel, str) and len(self.zlabel)==0: self.zlabel = "Density"

            # assign y values if it is a single column of data

            ycolindexlen = len(ycolindex)
            if ycolindexlen==0:
                self._yvalues = np.array( self._dfref().index[self.rows] + self._indexOffset ).flatten()
            elif ycolindexlen==1:
                self._yvalues = self._dfref().iloc[self.rows,ycolindex[0]].values.flatten()
            else:
                raise Exception ( "The length of ycolumns cannot be larger than 1 in jointplot or kdeplot2 objects." + newline
                                + self._getDocString()
                                )

            if self._type.isDiffusionPlot:

                if self._type.isContour or self._type.isContour3: self.currentFig.contourList = []
                if self._type.isContourf: self.currentFig.contourfList = []
                self.currentFig.kde2dList = []
                kde2dTriplet = Struct()

                if xcolindexlen<2 and ycolindexlen<2:

                    ( kde2dTriplet.density
                    , kde2dTriplet.grid
                    , kde2dTriplet.bandwidth
                    ) = kde2d( self._xvalues, self._yvalues, n=self.gridSize, limits=self.limits )

                    kde2dTriplet.density[ np.where( kde2dTriplet.density < self.noiseDensity ) ] = 0.0 # remove noise
                    kde2dTriplet.xmesh, kde2dTriplet.ymesh = np.meshgrid( kde2dTriplet.grid[0], kde2dTriplet.grid[1] )

        ############################################################################################################################
        #### make plot
        ############################################################################################################################

        if not self._type.isDiffusionPlot and self.legend.enabled: self.legend._labels = []

        if self._type.isHistplot:

            ########################################################################################################################
            #### make histplot
            ########################################################################################################################

            data = self._xvalues if xcolindexlen==1 else self._dfref().loc[:, xcolnames]
            self.currentFig.axes = sns.histplot ( data = data
                                                , **vars(self.histplot.kws)
                                                )

        for i, icol in enumerate(xcolindex):

            if xcolindexlen>1:
                self._xvalues = self._dfref().iloc[self.rows,xcolindex[i]].values.flatten()

            if (not self._type.is1d) and ycolindexlen>1: # order is important in condition
                self._yvalues = self._dfref().iloc[self.rows,ycolindex[i]].values.flatten()

            if (not self._type.isDiffusionPlot) and self.legend.enabled: # order is important in condition

                if xcolindexlen>1: lgxicol = i
                self.legend._labels.append(xcolnames[lgxicol])

                if not self._type.is1d:

                    if ycolindexlen>1: lgyicol = i
                    self.legend._labels[-1] += "-" + ycolnames[lgyicol]

            ########################################################################################################################
            #### make kdeplot1 / kdeplot2
            ########################################################################################################################

            if self._type.isKdeplot1 or self._type.isKdeplot2:

                if self._type.isKdeplot2: self.kdeplot.kws.data2 = self._yvalues

                self.currentFig.kdeplot = sns.kdeplot   ( data = self._xvalues
                                                        , **vars(self.kdeplot.kws)
                                                        )

            ########################################################################################################################
            #### make jointplot
            ########################################################################################################################

            if self._type.isJointplot:

                self.currentFig.jointplot = sns.jointplot   ( x = self._xvalues
                                                            , y = self._yvalues
                                                            , **vars(self.jointplot.kws)
                                                            )

            ########################################################################################################################
            #### make contour / contourf / contour3
            ########################################################################################################################

            if self._type.isDiffusionPlot:

                if xcolindexlen>1 or ycolindexlen>1:

                    ( kde2dTriplet.density
                    , kde2dTriplet.grid
                    , kde2dTriplet.bandwidth
                    ) = kde2d( self._xvalues, self._yvalues, n=self.gridSize, limits=self.limits )

                    kde2dTriplet.density[ np.where( kde2dTriplet.density < self.noiseDensity ) ] = 0.0 # remove noise
                    kde2dTriplet.xmesh, kde2dTriplet.ymesh = np.meshgrid( kde2dTriplet.grid[0]
                                                                        , kde2dTriplet.grid[1]
                                                                        )

                self.currentFig.kde2dList.append(kde2dTriplet)

                if self._type.isContour or self._type.isContour3:

                    self.currentFig.contourList.append( plt.contour ( kde2dTriplet.xmesh
                                                                    , kde2dTriplet.ymesh
                                                                    , kde2dTriplet.density
                                                                    , **vars(self.contour.kws)
                                                                    ) )

                if self._type.isContourf:

                    self.currentFig.contourfList.append(plt.contourf( kde2dTriplet.xmesh
                                                                    , kde2dTriplet.ymesh
                                                                    , kde2dTriplet.density
                                                                    , **vars(self.contourf.kws)
                                                                    ) )

        ############################################################################################################################
        #### add colorbar
        ############################################################################################################################

        if self._type.isDiffusionPlot and self.colorbar.enabled:

            self.colorbar.kws.mappable = None
            if self._type.isContour or self._type.isContour3:
                self.colorbar.kws.mappable = self.currentFig.contourList[0]
            elif self._type.isContourf:
                self.colorbar.kws.mappable = self.currentFig.contourfList[0]

            if self.colorbar.kws.mappable is not None:
                self.colorbar.kws.ax = self.currentFig.axes
                self.currentFig.colorbar = self.currentFig.figure.colorbar( **vars(self.colorbar.kws) )
                self.currentFig.colorbar.set_label( label = "Density" )

        ############################################################################################################################
        #### set axes labels
        ############################################################################################################################

        axkws = self.axes3d.kws if self._type.is3d else self.axes.kws
        xlabelNeeded = not (hasattr(axkws,"xlabel") and isinstance(axkws.xlabel, str))
        ylabelNeeded = not (hasattr(axkws,"ylabel") and isinstance(axkws.ylabel, str))

        if xlabelNeeded:

            xlabel = xcolnames[0] if len(xcolindex)==1 else "Variable Values"

            if self._type.isJointplot:
                self.currentFig.jointplot.ax_joint.set_xlabel(xlabel)
            else:
                self.currentFig.axes.set_xlabel(xlabel)

        if ylabelNeeded:

            if self._type.isHistplot:

                if self.histplot.kws.kde:
                    self.currentFig.axes.set_ylabel("Probability Density")
                else:
                    self.currentFig.axes.set_ylabel("Count")

            elif self._type.isKdeplot1:

                self.currentFig.axes.set_ylabel("Density")

            else:

                ylabel = ycolnames[0] if len(ycolindex)==1 else "Variable Values"
                if self._type.isJointplot:
                    self.currentFig.jointplot.ax_joint.set_ylabel(ylabel)
                else:
                    self.currentFig.axes.set_ylabel(ylabel)

        if self._type.is3d:

            zlabelNeeded = not (hasattr(axkws,"zlabel") and isinstance(axkws.zlabel, str))
            if zlabelNeeded: self.currentFig.axes.set_zlabel("Density")

        ############################################################################################################################
        #### set legend and other BasePlot properties
        ############################################################################################################################

        self._finalizeBasePlot()

        if not (self._type.is3d or self._type.isJointplot): self.target.currentFig.axes = self.currentFig.axes

        ############################################################################################################################

        if reself: return self

    ################################################################################################################################
    #### _getDocString
    ################################################################################################################################

    def _getDocString(self):
        docString   = newline \
                    + "Here is the help information on the DensityPlot class:" + newline \
                    + newline \
                    + self.__doc__ \
                    + super()._getDocString()
        return docString

    ################################################################################################################################
    #### helpme
    ################################################################################################################################

    def helpme(self, topic=None):
        """

        Print the documentation for the input string topic.
        If the topic does not exist, the documentation for
        the object will be printed.

            **Parameters**

                topic (optional)

                A string containing the name of the object
                for which help is needed.

            **Returns**

                None

            **Example**

                .. code-block:: python
                    :linenos:

                    helpme()
                    helpme("make")
                    helpme("helpme")
                    helpme("getLogLinSpace")

        """

        try:
            exec("print(self."+topic+".__doc__)")
        except:
            print(self._getDocString())
        return None

    ################################################################################################################################

