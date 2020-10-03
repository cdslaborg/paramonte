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
import weakref as wref
import _paramonte as pm
import _pmutils as pmutils
from paramonte.vis.Target import Target
from paramonte.vis._BasePlot import BasePlot

Struct = pmutils.Struct
newline = pmutils.newline

####################################################################################################################################
#### LineScatterPlot class
####################################################################################################################################

class LineScatterPlot(BasePlot):
    """

    This is the LineScatterPlot class for generating instances
    of line or scatter plots or the combination of the two in
    two or three dimensions based on the visualization tools
    of the ``matplotlib`` and ``seaborn`` Python libraries.

        **Usage**

            First generate an object of this class by optionally
            passing the following parameters described below. Then call
            the ``make()`` method. The generated object is also callable 
            with the same input parameters as the object's constructor.

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

            ycolumns

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

            zcolumns (exists only in 3D plot objects)

                An attribute that determines the columns of dataFrame 
                to be visualized as the Z-axis. It can have three forms:

                    1.  A list of column indices in dataFrame.
                    2.  A list of column names in dataFrame.columns.
                    3.  A ``range(start,stop,step)`` of column indices.

                Examples:

                    1.  ``zcolumns = [0,1,4,3]``
                    2.  ``zcolumns = ["SampleLogFunc","SampleVariable1"]``
                    3.  ``zcolumns = range(17,7,-2)``

                The default behavior includes all columns of the dataFrame.

            ccolumns

                An attribute that determines the columns of dataFrame 
                to be used for color mapping. It can have three forms:

                    1.  A list of column indices in dataFrame.
                    2.  A list of column names in dataFrame.columns.
                    3.  A ``range(start,stop,step)`` of column indices.

                Examples:

                    1.  ``ccolumns = [0,1,4,3]``
                    2.  ``ccolumns = ["SampleLogFunc","SampleVariable1"]``
                    3.  ``ccolumns = range(17,7,-2)``

                If ``ccolumns`` is set to ``None``, then no color-mapping 
                will be made. If it is set to an empty list ``[]``, then 
                the values from the ``rows`` attribute will be used for 
                color-mapping.

            rows

                An attribute that determines the rows of dataFrame to be 
                visualized. It can be either:

                    1.  A ``range(start,stop,step)``, or, 
                    2.  A list of row indices in dataFrame.index.

                Examples:

                    1.  ``rows = range(17,7,-2)``
                    2.  ``rows = [i for i in range(7,17)]``

                The default behavior includes all rows of the dataFrame.

            plot (exists only for line or lineScatter plots in 2D and 3D)

                A structure with two attributes:

                    enabled

                        A boolean indicating whether a call to the ``plot()``
                        function of the matplotlib library should be made 
                        or not.

                    kws

                        A structure whose components are directly passed as 
                        keyword arguments to the ``plot()`` function.

                        Example usage:

                            .. code-block:: python

                                plot.enabled = True
                                plot.kws.linewidth = 1

                        **NOTE**

                        If a desired property is missing among the ``kws`` 
                        attributes, simply add the field and its value to 
                        the component.

            scatter (exists only for scatter / lineScatter plots in 2D and 3D)

                A structure with two attributes:

                    enabled

                        A boolean indicating whether a call to the 
                        ``scatter()`` function of the matplotlib library 
                        should be made or not.

                    kws

                        A structure whose components are directly passed as 
                        keyword arguments to the ``scatter()`` function.

                        Example usage:

                            .. code-block:: python

                                scatter.enabled = True
                                scatter.kws.s = 2

                        **NOTE**

                        If a desired property is missing among the ``kws`` 
                        attributes, simply add the field and its value to 
                        the component.

            lineCollection (exists only for 2D / 3D line / lineScatter plots)

                A structure with two attributes:

                    enabled

                        A boolean indicating whether a call to the 
                        ``LineCollection()`` class of the matplotlib 
                        library should be made or not. This will result 
                        in line plots that are color-mapped.

                    kws

                        A structure whose components are directly passed as 
                        keyword arguments to the ``LineCollection()`` class.

                        Example usage:

                            .. code-block:: python

                                lineCollection.enabled = True
                                lineCollection.kws.linewidth = 1

                        **NOTE**

                        If a desired property is missing among the ``kws`` 
                        attributes, simply add the field and its value to 
                        the component.

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

                An object of class ``LineScatterPlot``.

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
        self._indexOffset       = 1
        self.xcolumns           = None
        self.ycolumns           = None
        self.ccolumns           = [] # indicates the use of default column for colormap. To turn off

        self._xlabel = None
        self._ylabel = None
        self._zlabel = None

        self._xlimit = None
        self._ylimit = None
        self._zlimit = None

        self._xscale = None
        self._yscale = None
        self._zscale = None

        if not self._type.is3d: self.target = Target()

        if self._type.isLine:
            self.plot = Struct()
            self.plot.enabled = False
            self.plot.kws = Struct()
            self.lineCollection = Struct()
            self.lineCollection.kws = Struct()
            self.lineCollection.enabled = True

        if self._type.isScatter:
            self.scatter = Struct()
            self.scatter.enabled = True
            self.scatter.kws = Struct()

        if self._type.is3d:
            self.zcolumns = []

        self.colorbar = Struct()
        self.colorbar.kws = Struct()
        if (self._type.isScatter and self.ccolumns is not None) or (self._type.isLine and self.lineCollection.enabled):
            self.colorbar.enabled = True
        else:
            self.colorbar.enabled = False

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

        Call the ``make()`` method of the current 
        instance of the class.

            **Parameters**

                Any arguments that can be passed to the 
                ``make()`` method of the plot object.

            **Returns**

                Any return value from the ``make()`` 
                method of the plot object.

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

        Generate a line/scatter plot from the 
        selected columns of the object's dataframe.

            **Parameters**

                reself

                    A logical variable. If ``True``, an instance of 
                    the object will be returned to the calling routine 
                    upon exit. The default value is ``False``.

            **Returns**

                The object self if ``reself = True`` otherwise, ``None``.
                However, this method causes side-effects by manipulating
                the existing attributes of the object.

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

        cEnabled = self.ccolumns is not None

        from collections.abc import Iterable
        if self.ccolumns is not None and not isinstance(self.ccolumns, Iterable): self.ccolumns = [self.ccolumns]

        # if no colormap, then 

        if self._type.isLine and not cEnabled: self.plot.enabled = True

        ############################################################################################################################
        #### scatter plot properties
        ############################################################################################################################

        if self._type.isScatter:
            if isinstance(self.scatter.kws,Struct):

                if "s" not in vars(self.scatter.kws).keys(): self.scatter.kws.s = 2
                if "c" not in vars(self.scatter.kws).keys(): self.scatter.kws.c = None
                if "cmap" not in vars(self.scatter.kws).keys() or self.scatter.kws.cmap is None: self.scatter.kws.cmap = "autumn"
                if "alpha" not in vars(self.scatter.kws).keys(): self.scatter.kws.alpha = 1
                if "edgecolors" not in vars(self.scatter.kws).keys(): self.scatter.kws.edgecolors = None
                if "zorder" not in vars(self.scatter.kws).keys(): self.scatter.kws.zorder = 2

                if not cEnabled: self.scatter.kws.cmap = None

            else:

                raise Exception ( "The scatter.kws component of the current LineScatterPlot object must" + newline
                                + "be an object of class Struct(), essentially a structure with components" + newline
                                + "whose names are the input arguments to the scatter() function of the" + newline
                                + "matplotlib library." + newline
                                + self._getDocString()
                                )

        ############################################################################################################################
        #### line plot properties
        ############################################################################################################################

        if self._type.isLine:

            if isinstance(self.plot.kws,Struct):

                if "linewidth" in vars(self.plot.kws).keys():
                    if self.plot.kws.linewidth==0: self.plot.kws.linewidth = 1
                else:
                    self.plot.kws.linewidth = 1

                if "zorder" not in vars(self.plot.kws).keys(): self.plot.kws.zorder = 1

            else:

                raise Exception ( "The plot.kws component of the current LineScatterPlot object must" + newline
                                + "be an object of class Struct(), essentially a structure with components" + newline
                                + "whose names are the input arguments to the plot() function of the" + newline
                                + "matplotlib library." + newline
                                + self._getDocString()
                                )

            if isinstance(self.lineCollection.kws, Struct):

                if "cmap" not in vars(self.lineCollection.kws).keys() or self.lineCollection.kws.cmap is None: self.lineCollection.kws.cmap = "autumn"
                if "alpha" not in vars(self.lineCollection.kws).keys(): self.lineCollection.kws.alpha = 1
                if "linewidth" not in vars(self.lineCollection.kws).keys(): self.lineCollection.kws.linewidth = 1

            else:

                objectType = "LineCollection"
                if self._type.is3d: objectType = "Line3DCollection"
                raise Exception ( "The lineCollection.kws component of the current LineScatterPlot object must" + newline
                                + "be an object of class Struct(), essentially a structure with components" + newline
                                + "whose names are the input arguments to the " + objectType + "() class of the" + newline
                                + "matplotlib library." + newline
                                + self._getDocString()
                                )

        ############################################################################################################################
        #### legend properties
        ############################################################################################################################

        if self.legend.enabled:
            if not isinstance(self.legend.kws,Struct):
                raise Exception ( "The legend.kws component of the current LineScatterPlot object must" + newline
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
                raise Exception ( "The figure.kws component of the current LineScatterPlot object must" + newline
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

        from matplotlib import pyplot as plt
        from matplotlib.collections import LineCollection
        if self._type.is3d: from mpl_toolkits.mplot3d.art3d import Line3DCollection

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

        if self.xcolumns is None:
            lgxicol = 0
            xcolindex = []
            xcolnames = ["Count"]
            if self._type.isScatter : self.scatter._xvalues = np.array( self._dfref().index[self.rows] + self._indexOffset ).flatten()
            if self._type.isLine    : self.plot._xvalues    = np.array( self._dfref().index[self.rows] + self._indexOffset ).flatten()
        else:
            xcolnames, xcolindex = pm.dfutils.getColNamesIndex(self._dfref().columns,self.xcolumns)

        if self.ycolumns is None:
            lgyicol = 0
            ycolindex = []
            ycolnames = ["Count"]
            if self._type.isScatter : self.scatter._yvalues = np.array( self._dfref().index[self.rows] + self._indexOffset ).flatten()
            if self._type.isLine    : self.plot._yvalues    = np.array( self._dfref().index[self.rows] + self._indexOffset ).flatten()
        else:
            ycolnames, ycolindex = pm.dfutils.getColNamesIndex(self._dfref().columns,self.ycolumns)

        if self._type.is3d:
            if self.zcolumns is None:
                lgzicol = 0
                zcolindex = []
                zcolnames = ["Count"]
                if self._type.isScatter : self.scatter._zvalues = np.array( self._dfref().index[self.rows] + self._indexOffset ).flatten()
                if self._type.isLine    : self.plot._zvalues    = np.array( self._dfref().index[self.rows] + self._indexOffset ).flatten()
            else:
                zcolnames, zcolindex = pm.dfutils.getColNamesIndex(self._dfref().columns, self.zcolumns)

        ############################################################################################################################
        #### set colormap data
        ############################################################################################################################

        if cEnabled:
            if len(self.ccolumns)==0:
                ccolindex = []
                ccolnames = ["Count"]
                if self._type.isScatter : self.scatter.kws.c    = np.array( self._dfref().index[self.rows] + self._indexOffset ).flatten()
                if self._type.isLine    : cdata                 = np.array( self._dfref().index[self.rows] + self._indexOffset ).flatten()
            else:
                ccolnames, ccolindex = pm.dfutils.getColNamesIndex(self._dfref().columns,self.ccolumns)
        else:
            ccolindex = []
            ccolnames = []
            #self.scatter.kws.c = None

        ############################################################################################################################
        #### check the consistency of the lengths
        ############################################################################################################################

        xcolindexlen = len(xcolindex)
        ycolindexlen = len(ycolindex)
        ccolindexlen = len(ccolindex)

        maxLenColumns = np.max  (   [ xcolindexlen
                                    , ycolindexlen
                                    , ccolindexlen
                                    ]
                                )

        if xcolindexlen!=maxLenColumns and xcolindexlen>1: raise Exception("length of xcolumns must be either 1 or equal to the lengths of ycolumns or ccolumns.")
        if ycolindexlen!=maxLenColumns and ycolindexlen>1: raise Exception("length of ycolumns must be either 1 or equal to the lengths of xcolumns or ccolumns.")
        if ccolindexlen!=maxLenColumns and ccolindexlen>1: raise Exception("length of ccolumns must be either 1 or equal to the lengths of xcolumns or ycolumns.")

        if self._type.is3d:
            zcolindexlen = len(zcolindex)
            if zcolindexlen!=maxLenColumns and zcolindexlen>1: raise Exception("length of zcolumns must be either 1 or equal to the lengths of xcolumns or ycolumns.")

        ############################################################################################################################
        #### assign data in case of single column assignments
        ############################################################################################################################

        if xcolindexlen==1:
            lgxicol = 0
            if self._type.isScatter     : self.scatter._xvalues = self._dfref().iloc[self.rows,xcolindex].values.flatten()
            if self._type.isLine        : self.plot._xvalues    = self._dfref().iloc[self.rows,xcolindex].values.flatten()

        if ycolindexlen==1:
            lgyicol = 0
            if self._type.isScatter     : self.scatter._yvalues = self._dfref().iloc[self.rows,ycolindex].values.flatten()
            if self._type.isLine        : self.plot._yvalues    = self._dfref().iloc[self.rows,ycolindex].values.flatten()

        if self._type.is3d:
            if zcolindexlen==1:
                lgzicol = 0
                if self._type.isScatter : self.scatter._zvalues = self._dfref().iloc[self.rows,zcolindex].values.flatten()
                if self._type.isLine    : self.plot._zvalues    = self._dfref().iloc[self.rows,zcolindex].values.flatten()

        if cEnabled:
            if ccolindexlen==1:
                if self._type.isScatter : self.scatter.kws.c    = self._dfref().iloc[self.rows,ccolindex].values.flatten()
                if self._type.isLine    : cdata                 = self._dfref().iloc[self.rows,ccolindex].values.flatten()

        ############################################################################################################################
        #### add line/scatter plot
        ############################################################################################################################

        if self.legend.enabled: self.legend._labels = []
        if self._type.isScatter and self.scatter.enabled: self.currentFig.scatterList = []
        if self._type.isLine:
            if self.plot.enabled: self.currentFig.plotList = []
            if cEnabled and self.lineCollection.enabled:
                self.currentFig.lineCollectionList = []
                self.currentFig.plotList = []

        for i in range(maxLenColumns):

            if xcolindexlen>1:
                lgxicol = i
                if self._type.isScatter     : self.scatter._xvalues = self._dfref().iloc[self.rows,xcolindex[i]].values.flatten()
                if self._type.isLine        : self.plot._xvalues    = self._dfref().iloc[self.rows,xcolindex[i]].values.flatten()

            if ycolindexlen>1:
                lgyicol = i
                if self._type.isScatter     : self.scatter._yvalues = self._dfref().iloc[self.rows,ycolindex[i]].values.flatten()
                if self._type.isLine        : self.plot._yvalues    = self._dfref().iloc[self.rows,ycolindex[i]].values.flatten()

            if cEnabled:
                if ccolindexlen>1:
                    if self._type.isScatter : self.scatter.kws.c    = self._dfref().iloc[self.rows,ccolindex[i]].values.flatten()
                    if self._type.isLine    : cdata                 = self._dfref().iloc[self.rows,ccolindex[i]].values.flatten()

            if self.legend.enabled:
                if xcolindexlen<2 and ycolindexlen>1:
                    self.legend._labels.append(ycolnames[lgyicol])
                elif xcolindexlen>1 and ycolindexlen<2:
                    self.legend._labels.append(xcolnames[lgxicol])
                else:
                    self.legend._labels.append( xcolnames[lgxicol] + "-" + ycolnames[lgyicol] )

            if self._type.is3d:
                if zcolindexlen>1:
                    lgzicol = i
                    if self._type.isScatter : self.scatter._zvalues = self._dfref().iloc[self.rows,zcolindex[i]].values.flatten()
                    if self._type.isLine    : self.plot._zvalues    = self._dfref().iloc[self.rows,zcolindex[i]].values.flatten()
                if self.legend.enabled:
                    if zcolindexlen>1: self.legend._labels[-1] += "-" + zcolnames[lgzicol]

            ########################################################################################################################
            #### add scatter plot
            ########################################################################################################################

            if self._type.isScatter and self.scatter.enabled:

                if self._type.is3d:

                    self.currentFig.scatterList.append( self.currentFig.axes.scatter( self.scatter._xvalues
                                                                                    , self.scatter._yvalues
                                                                                    , self.scatter._zvalues
                                                                                    , **vars(self.scatter.kws)
                                                                                    ) )

                else:

                    self.currentFig.scatterList.append( self.currentFig.axes.scatter( self.scatter._xvalues
                                                                                    , self.scatter._yvalues
                                                                                    , **vars(self.scatter.kws)
                                                                                    ) )

            ########################################################################################################################
            #### add line plot
            ########################################################################################################################

            if self._type.isLine:

                if self.plot.enabled:

                    if self._type.is3d:

                        self.currentFig.plotList.append ( self.currentFig.axes.plot ( self.plot._xvalues
                                                                                    , self.plot._yvalues
                                                                                    , self.plot._zvalues
                                                                                    , **vars(self.plot.kws)
                                                                                    ) )

                    else:

                        self.currentFig.plotList.append ( self.currentFig.axes.plot ( self.plot._xvalues
                                                                                    , self.plot._yvalues
                                                                                    , **vars(self.plot.kws)
                                                                                    ) )

                if cEnabled and self.lineCollection.enabled:

                    self.lineCollection.kws.norm = norm = plt.Normalize(cdata.min(), cdata.max())

                    if self._type.is3d:

                        # properly and automatically set the axes limits via plot()

                        self.currentFig.plotList.append ( self.currentFig.axes.plot ( self.plot._xvalues
                                                                                    , self.plot._yvalues
                                                                                    , self.plot._zvalues
                                                                                    , linewidth = 0
                                                                                    ) )

                        # now add the lineCollection

                        points = np.array([self.plot._xvalues, self.plot._yvalues, self.plot._zvalues]).T.reshape(-1, 1, 3)
                        segments = np.concatenate([points[:-1], points[1:]], axis=1)
                        lineCollection = Line3DCollection( segments, **vars(self.lineCollection.kws) )

                    else:

                        # properly and automatically set the axes limits via plot()

                        self.currentFig.plotList.append ( self.currentFig.axes.plot ( self.plot._xvalues
                                                                                    , self.plot._yvalues
                                                                                    , linewidth = 0
                                                                                    ) )

                        # now add the lineCollection

                        points = np.array([self.plot._xvalues, self.plot._yvalues]).T.reshape(-1, 1, 2)
                        segments = np.concatenate([points[:-1], points[1:]], axis=1)
                        lineCollection = LineCollection( segments, **vars(self.lineCollection.kws) )

                    lineCollection.set_array(cdata)
                    #lineCollection.set_linewidth(0.5)
                    #lineCollection.set_solid_capstyle("round")
                    self.currentFig.lineCollectionList.append( self.currentFig.axes.add_collection(lineCollection) )

        ############################################################################################################################
        #### add colorbar
        ############################################################################################################################

        cbarEnabled = cEnabled and self.colorbar.enabled and (ccolindexlen<2) # and (not hasattr(self.currentFig,"colorbar"))
        if cbarEnabled:

            self.colorbar.kws.mappable = None
            if self._type.isLine and self.lineCollection.enabled:

                self.colorbar.kws.mappable = self.currentFig.lineCollectionList[0]

            elif self._type.isScatter and self.scatter.enabled:

                self.colorbar.kws.mappable = self.currentFig.scatterList[0]

            if self.colorbar.kws.mappable is not None:
                self.colorbar.kws.ax = self.currentFig.axes
                self.currentFig.colorbar = self.currentFig.figure.colorbar( **vars(self.colorbar.kws) )
                self.currentFig.colorbar.set_label( label = ", ".join(ccolnames) )

        ############################################################################################################################
        #### set axes scales
        ############################################################################################################################

        if self._xscale is not None: self.currentFig.axes.set_xscale(self._xscale)

        if self._yscale is not None: self.currentFig.axes.set_yscale(self._yscale)

        if self._zscale is not None and self._type.is3d: self.currentFig.axes.set_zscale(self._zscale)

        ############################################################################################################################
        #### set axes limits
        ############################################################################################################################

        if self._xlimit is not None:
            currentLim = list(self.currentFig.axes.get_xlim())
            if self._xlimit[0] is not None: currentLim[0] = self._xlimit[0]
            if self._xlimit[1] is not None: currentLim[1] = self._xlimit[1]
            self.currentFig.axes.set_xlim(currentLim)

        if self._ylimit is not None:
            currentLim = list(self.currentFig.axes.get_ylim())
            if self._ylimit[0] is not None: currentLim[0] = self._ylimit[0]
            if self._ylimit[1] is not None: currentLim[1] = self._ylimit[1]
            self.currentFig.axes.set_ylim(currentLim)

        if self._zlimit is not None and self._type.is3d:
            currentLim = list(self.currentFig.axes.get_zlim())
            if self._zlimit[0] is not None: currentLim[0] = self._zlimit[0]
            if self._zlimit[1] is not None: currentLim[1] = self._zlimit[1]
            self.currentFig.axes.set_zlim(currentLim)

        ############################################################################################################################
        #### add axes labels
        ############################################################################################################################

        if self._xlabel is None:
            if xcolindexlen>1:
                self.currentFig.axes.set_xlabel("Variable Values")
            else:
                self.currentFig.axes.set_xlabel(xcolnames[0])
        else:
            self.currentFig.axes.set_xlabel(self._xlabel)

        if self._ylabel is None:
            if ycolindexlen>1:
                self.currentFig.axes.set_ylabel("Variable Values")
            else:
                self.currentFig.axes.set_ylabel(ycolnames[0])
        else:
            self.currentFig.axes.set_ylabel(self._ylabel)

        if self._type.is3d:
            if self._zlabel is None:
                if zcolindexlen>1:
                    self.currentFig.axes.set_zlabel("Variable Values")
                else:
                    self.currentFig.axes.set_zlabel(zcolnames[0])
            else:
                self.currentFig.axes.set_zlabel(self._zlabel)

        ############################################################################################################################
        #### set legend and other BasePlot properties
        ############################################################################################################################

        self._finalizeBasePlot()

        if not self._type.is3d: self.target.currentFig.axes = self.currentFig.axes

        ############################################################################################################################

        if reself: return self

    ################################################################################################################################
    #### _getDocString
    ################################################################################################################################

    def _getDocString(self):
        docString   = newline \
                    + "Here is the help information on the LineScatterPlot class:" + newline \
                    + newline \
                    + self.__doc__ \
                    + newline \
                    + "Here is the help information on the parent BasePlot class:" + newline \
                    + newline \
                    + super().__doc__
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
