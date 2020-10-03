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
#import pandas as pd
import weakref as wref
from copy import deepcopy as dcopy

import _paramonte as pm
from paramonte.vis._BasePlot import BasePlot
from paramonte.vis.LineScatterPlot import LineScatterPlot
from paramonte.vis.DensityPlot import DensityPlot

Struct = pm.Struct
newline = pm.newline

####################################################################################################################################
#### GridPlot class
####################################################################################################################################

class GridPlot(BasePlot):
    """

    This is the GridPlot class for generating instances
    of histogram and contour figures in two and three dimensions
    based on a wide range of plotting tools from the matplotlib,
    seaborn, and other Python libraries.

    Normally, the public users are not supposed to use this class directly,
    although they can for the purposes other than plotting the ParaMonte
    simulation files.

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

            All of the parameters described above, are implicitly 
            stored in the object.

            columns

                An attribute that determines the columns of dataFrame to be 
                visualized. It can have three forms:

                    1.  A list of column indices in dataFrame.
                    2.  A list of column names in dataFrame.columns.
                    3.  A ``range(start,stop,step)`` of column indices.

                Examples:

                    1.  ``xcolumns = [0,1,4,3]``
                    2.  ``xcolumns = ["SampleLogFunc","SampleVariable1"]``
                    3.  ``xcolumns = range(17,7,-2)``

                The default behavior includes up to 5 columns of the dataFrame.

            ccolumn

                A string or integer that determines the column 
                of the dataFrame to serve as the color-values 
                corresponding to each point in the plot. 

                Examples:

                    1.  ``ccolumn = 7``
                    2.  ``ccolumn = "SampleLogFunc"``
                    3.  ``ccolumn = None``

                **NOTE**

                The value of ``ccolumns`` has precedence over all 
                other color properties of individual plotting tools, 
                unless ``ccolumns = None``.
                For example, any value of the ``c`` component of the ``scatter.kws``
                component of the ``LineScatterPlot`` object will be overridden by
                the specified value for ``ccolumns``.

            rows

                An attribute that determines the rows of dataFrame to be 
                visualized. It can be either:

                    1.  A ``range(start,stop,step)``, or, 
                    2.  A list of row indices in dataFrame.index.

                Examples:

                    1.  ``rows = range(17,7,-2)``
                    2.  ``rows = [i for i in range(7,17)]``

                The default behavior includes all rows of the dataFrame.

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

            figure

                A structure with two attributes:

                    enabled

                        A boolean indicating whether a call to the ``figure()``
                        function of the matplotlib library should be made or not.
                        If a call is made, a new figure will be generated.
                        Otherwise, the current active figure will be used.

                    kws

                        A structure whose components are directly passed as 
                        keyword arguments to the ``figure()`` function of 
                        the matplotlib library.

                        Example usage:

                            .. code-block:: python

                                figure.kws.facecolor = "w"

                        **NOTE**

                        If a desired property is missing among the ``kws`` 
                        attributes, simply add the field and its value to 
                        the component.

            pairgrid

                A structure with one attribute:

                    kws

                        A structure whose components are directly passed as 
                        keyword arguments to the ``PairGrid()`` class of 
                        the seaborn library to generate a pairgrid.

                        Example usage:

                            .. code-block:: python

                                pairgrid.kws.dropna = True

                        **NOTE**

                        Passing some keywords, like ``corner`` attribute,
                        to the PairGrid will cause potential conflicts with
                        other attributes of ParaMonte's GridPlot. Use this
                        attribute with caution. That said, if a desired 
                        property is missing among the ``kws`` attributes, 
                        simply add the field and its value to the component.

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

            plotType

                A structure with three attributes:

                    upper

                        A structure whose components determine what plots
                        can and should be added to the upper corner of 
                        the grid plot.

                    lower

                        A structure whose components determine what plots
                        can and should be added to the lower corner of 
                        the grid plot.

                    diag

                        A structure whose components determine what plots
                        can and should be added to the diagonal of 
                        the grid plot.

                Each of the above structures has three components:

                    enabled

                        A boolean that determines if the subplots in the
                        corresponding section of the grid plot should be 
                        added or not.

                    names

                        A list of strings that represent the types of plots
                        that can be added to the corresponding section of 
                        the grid plot. 

                    value

                        A string whose value can be one of the elements of 
                        the ``names`` attribute explained in the above. It 
                        determines the type of the subplot to be added to 
                        the corresponding section of the grid plot. 

                Example usage:

                    .. code-block:: python

                        plotType.upper.enabled = False # no upper corner plots
                        plotType.lower.value = "scatter" # add scatter to lower

            layout

                A structure whose components are directly passed to the 
                corresponding plotting tools of the ParaMonte library to 
                draw the corresponding subplots, if they are activated.

                Example usage:

                    .. code-block:: python

                        layout.contour.contour.kws.colors = "blue"

            currentFig

                A structure whose attributes are the outputs of various plotting 
                tools used to make the current figure. These include the handle 
                to the current figure, the handle to the current axes in the plot, 
                the handle to the colorbar (if any exists), and other Python 
                plotting tools used to make to generate the figure.

        **Returns**

            self

                An object of class ``GridPlot``.

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

        self.columns = None
        self.ccolumn = [] # use the default colormap data (which is the data count)

        self.pairgrid = Struct()
        self.pairgrid.kws = Struct()

        self.colorbar = Struct()
        self.colorbar.enabled = True
        self.colorbar.kws = Struct()

        ############################################################################################################################
        #### setup corner types
        ############################################################################################################################

        cornerNames =   [ "line"
                        , "scatter"
                        , "lineScatter"
                        , "contourf"
                        , "contour"
                        ]

        self.plotType = Struct()

        self.plotType.diag = Struct()
        self.plotType.diag.enabled = True
        self.plotType.diag.value = "histplot"
        self.plotType.diag.names = ["histplot"]

        self.plotType.upper = Struct()
        self.plotType.upper.enabled = True
        self.plotType.upper.value = "lineScatter"
        self.plotType.upper.names = dcopy(cornerNames)

        self.plotType.lower = Struct()
        self.plotType.lower.enabled = True
        self.plotType.lower.value = "contour"
        self.plotType.lower.names = dcopy(cornerNames)

        ############################################################################################################################
        #### setup subplot templates
        ############################################################################################################################

        template = Struct()

        # target

        template.target = Struct()
        template.target.axvline = Struct()
        template.target.axvline.kws = Struct()
        template.target.axhline = Struct()
        template.target.axhline.kws = Struct()
        template.target.scatter = Struct()
        template.target.scatter.kws = Struct()

        template.target.axvline.kws.linewidth = 0.5
        template.target.axvline.kws.linestyle = "-"
        template.target.axvline.kws.zorder = 1000
        template.target.axvline.kws.color = "orangered"

        template.target.axhline.kws.linewidth = 0.5
        template.target.axhline.kws.linestyle = "-"
        template.target.axhline.kws.zorder = 1000
        template.target.axhline.kws.color = "orangered"

        template.target.scatter.kws.s = 20
        template.target.scatter.kws.color = "orangered"
        template.target.scatter.kws.zorder = 1002

        # contour

        template.contour = Struct()
        template.contour.enabled = True
        template.contour.kws = Struct()
        template.contour.kws.cmap = "winter"
        template.contour.kws.alpha = 1
        template.contour.kws.levels = 10
        template.contour.kws.linewidths = 1
        template.contour.kws.linestyles = "solid"

        # contourf

        template.contourf = Struct()
        template.contourf.enabled = True
        template.contourf.kws = Struct()
        template.contourf.kws.cmap = "Blues"
        template.contourf.kws.alpha = 1
        template.contourf.kws.levels = 50

        # histplot

        template.histplot = Struct()
        template.histplot.enabled = True
        template.histplot.kws = Struct()
        template.histplot.kws.kde = False
        template.histplot.kws.bins = "auto"
        template.histplot.kws.stat = "count"
        template.histplot.kws.kde_kws = dict()
        template.histplot.kws.binwidth = None
        template.histplot.kws.binrange = None
        template.histplot.kws.multiple = "stack"
        template.histplot.kws.element = "step"
        template.histplot.kws.legend = False
        template.histplot.kws.shrink = 1
        template.histplot.kws.color = None
        template.histplot.kws.fill = True
        template.histplot.kws.common_norm = True
        template.histplot.kws.common_bins = False
        template.histplot.kws.line_kws = dict()
        template.histplot.kws.line_kws["linewidth"] = 0
        template.histplot.kws.line_kws["linestyle"] = "-"

        # legend template

        template.legend = Struct()
        template.legend.enabled = False
        template.legend.kws = Struct()

        # colorbar template

        template.colorbar = Struct()
        template.colorbar.enabled = False
        template.colorbar.kws = Struct()

        # line / scatter / lineScatter

        template.scatter = Struct()
        template.scatter.enabled = True
        template.scatter.kws = Struct()
        template.scatter.kws.s = 2
        template.scatter.kws.c = None
        template.scatter.kws.cmap = "winter"
        template.scatter.kws.alpha = 1
        template.scatter.kws.edgecolor = None
        template.scatter.kws.zorder = 2

        template.plot = Struct()
        template.plot.enabled = False
        template.plot.kws = Struct()
        template.plot.kws.linewidth = 1
        template.plot.kws.zorder = 1

        template.lineCollection = Struct()
        template.lineCollection.enabled = True
        template.lineCollection.kws = Struct()
        template.lineCollection.kws.cmap = "winter"
        template.lineCollection.kws.alpha = 1
        template.lineCollection.kws.linewidth = 1

        # figure

        template.figure = Struct()
        template.figure.kws = Struct()
        template.figure.enabled = False

        ############################################################################################################################
        #### setup subplots layout
        ############################################################################################################################

        self.layout = Struct()

        self.layout.contour                     = Struct()
        self.layout.contour.figure              = dcopy(template.figure)
        self.layout.contour.contour             = dcopy(template.contour)
        self.layout.contour.colorbar            = dcopy(template.colorbar)
        #self.layout.contour.target              = dcopy(template.target)
        self.layout.contour.noiseDensity        = 1.e-3
        self.layout.contour.limits              = None
        self.layout.contour.gridSize            = 512

        self.layout.contourf                    = Struct()
        self.layout.contourf.figure             = dcopy(template.figure)
        self.layout.contourf.contourf           = dcopy(template.contourf)
        self.layout.contourf.colorbar           = dcopy(template.colorbar)
        #self.layout.contourf.target             = dcopy(template.target)
        self.layout.contourf.noiseDensity       = 1.e-5
        self.layout.contourf.limits             = None
        self.layout.contourf.gridSize           = 512

        self.layout.histplot                    = Struct()
        self.layout.histplot.figure             = dcopy(template.figure)
        self.layout.histplot.histplot           = dcopy(template.histplot)
        self.layout.histplot.legend             = dcopy(template.legend)
        #self.layout.histplot.target             = dcopy(template.target)

        self.layout.line                        = Struct()
        self.layout.line.figure                 = dcopy(template.figure)
        self.layout.line.plot                   = dcopy(template.plot)
        self.layout.line.lineCollection         = dcopy(template.lineCollection)
        self.layout.line.colorbar               = dcopy(template.colorbar)
        self.layout.line.legend                 = dcopy(template.legend)
        #self.layout.line.target                 = dcopy(template.target)

        self.layout.scatter                     = Struct()
        self.layout.scatter.figure              = dcopy(template.figure)
        self.layout.scatter.scatter             = dcopy(template.scatter)
        self.layout.scatter.colorbar            = dcopy(template.colorbar)
        self.layout.scatter.legend              = dcopy(template.legend)
        #self.layout.scatter.target              = dcopy(template.target)

        self.layout.lineScatter                 = Struct()
        self.layout.lineScatter.figure          = dcopy(template.figure)
        self.layout.lineScatter.plot            = dcopy(template.plot)
        self.layout.lineScatter.scatter         = dcopy(template.scatter)
        self.layout.lineScatter.lineCollection  = dcopy(template.lineCollection)
        self.layout.lineScatter.colorbar        = dcopy(template.colorbar)
        self.layout.lineScatter.legend          = dcopy(template.legend)
        #self.layout.lineScatter.target          = dcopy(template.target)
        self.layout.lineScatter.plot.enabled    = True
        self.layout.lineScatter.plot.kws.alpha  = 0.2
        self.layout.lineScatter.plot.kws.color  = "grey"
        self.layout.lineScatter.plot.kws.linewidth = 0.75
        self.layout.lineScatter.lineCollection.enabled = False

        ############################################################################################################################

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

        Generate a grid plot from the selected 
        columns of the object's dataFrame.

            **Parameters**

                reself

                    A logical variable. If ``True``, an instance of
                    the object will be returned  to the calling routine
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

        self._cEnabled = self.ccolumn is not None

        ############################################################################################################################
        #### verify the plot types to draw
        ############################################################################################################################

        if self.plotType.upper.enabled and (not (isinstance(self.plotType.upper.value,str) and self.plotType.upper.value in self.plotType.upper.names)):
            raise Exception ( "Unrecognized input value for the \"plotType.upper.value\" of the GridPlot object." + newline
                            + "The possible plot type names are given in \"plotType.upper.value\" of the object." + newline
                            + self._getDocString()
                            )

        if self.plotType.lower.enabled and (not (isinstance(self.plotType.lower.value,str) and self.plotType.lower.value in self.plotType.lower.names)):
            raise Exception ( "Unrecognized input value for the \"plotType.lower.value\" of the GridPlot object." + newline
                            + "The possible plot type names are given in \"plotType.lower.value\" of the object." + newline
                            + self._getDocString()
                            )

        if self.plotType.diag.enabled and (not (isinstance(self.plotType.diag.value,str) and self.plotType.diag.value in self.plotType.diag.names)):
            raise Exception ( "Unrecognized input value for the \"plotType.diag.value\" of the GridPlot object." + newline
                            + "The possible plot type names are given in \"plotType.diag.value\" of the object." + newline
                            + self._getDocString()
                            )

        ############################################################################################################################
        ############################################################################################################################
        if self._isdryrun: return
        ############################################################################################################################
        ############################################################################################################################

        import seaborn as sns
        import matplotlib.pyplot as plt

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

        # assign x columns to plot

        #if isinstance(self.columns,list):
        try:
            self._colnames, self._colindex = pm.dfutils.getColNamesIndex(self._dfref().columns,self.columns)
            colindexlen = len(self._colindex)
        except:
            raise Exception ( "The columns component of the current GridPlot object must point to" + newline
                            + "the names of a set of the columns of the input dataFrame to the GridPlot" + newline
                            + "class constructor." + newline
                            + self._getDocString()
                            )

        if colindexlen==0:
            raise Exception ( "The length of the ``columns`` component of the GridPlot object cannot be zero." + newline
                            + self._getDocString()
                            )

        # set color data

        if self._cEnabled:
            if isinstance(self.ccolumn,str) or isinstance(self.ccolumn,int):
                self._ccolname, self._ccolindex = pm.dfutils.getColNamesIndex(self._dfref().columns,self.ccolumn)
                #self._cdata = self._dfref().iloc[self.rows,self._colindex[0]].values.flatten()
            elif isinstance(self.ccolumn,list) and len(self.ccolumn)==0:
                self._ccolname = "Count"
                self._ccolindex = []
            else:
                raise Exception ( "The ccolumn component of the current GridPlot object must point to" + newline
                                + "the names of a column of the input dataFrame to the GridPlot" + newline
                                + "class constructor. It represents the set of values that will " + newline
                                + "be used to colormap the subplots, where needed." + newline
                                + self._getDocString()
                                )

            self._ccolindexlen = len(self._ccolindex)
            if self._ccolindexlen>1:
                raise Exception ( "The length of the ``ccolumn`` component of the GridPlot object cannot be larger than 1." + newline
                                + self._getDocString()
                                )


        ############################################################################################################################
        #### make pairgrid plot
        ############################################################################################################################

        self.currentFig.pairgrid = sns.PairGrid( self._dfref()[self._colnames].iloc[self.rows,:], **vars(self.pairgrid.kws) )
        self.currentFig.figure = self.currentFig.pairgrid.fig

        ############################################################################################################################
        #### hide axes as requested
        ############################################################################################################################

        if not self.plotType.upper.enabled: self.hide("upper")
        if not self.plotType.lower.enabled: self.hide("lower")
        if not self.plotType.diag.enabled: self.hide("diag")

        ############################################################################################################################
        #### adjust subplot interspaces
        ############################################################################################################################

        self.currentFig.pairgrid.fig.tight_layout()
        self.currentFig.pairgrid.fig.subplots_adjust(hspace = 0.05, wspace = 0.05)

        ############################################################################################################################
        #### add subplots
        ############################################################################################################################

        self._nrow = len(self.currentFig.pairgrid.axes[:])
        nrowMinusOne = self._nrow - 1
        self._ntot = self._nrow**2
        self._ntotStr = str(self._ntot)
        self.currentFig.subplotGrid = [ [ [] for j in range(self._nrow) ] for i in range(self._nrow) ]

        counter = 0
        dataMin = np.zeros(self._nrow)
        dataMax = np.zeros(self._nrow)
        dataMargin = np.zeros(self._nrow)

        ############################################################################################################################
        #### add the diagonal subplots
        ############################################################################################################################

        avgDiagTime = 0 # must be defined here for later default use, when diag is disabled.
        if self.plotType.diag.enabled:
            #self._progress.note( msg = "generating diagonal subplots... ", end = "" )
            self._progress.timer.toc()
            self.currentFig.pairgrid.map_diag( sns.histplot, **vars(self.layout.histplot.histplot.kws) )
            self._progress.timer.toc()
            avgDiagTime = self._progress.timer.delta / self._nrow
            #if self._reportEnabled:
            #    self._progress.timer.toc()
            #    print( "done in " + str(np.round(self._progress.timer.delta,6)) + " seconds." )

        ############################################################################################################################
        #### add the off-diagonal subplots
        ############################################################################################################################

        for irow, axRow in enumerate(self.currentFig.pairgrid.axes[:]):

            # compute data range

            #dataMin[irow] = np.min( self._dfref().iloc[self.rows,self._colindex[irow]].values.flatten() )
            #dataMax[irow] = np.max( self._dfref().iloc[self.rows,self._colindex[irow]].values.flatten() )
            #dataMargin[irow] = (dataMax[irow] - dataMin[irow]) / 10

            for icol, ax in enumerate(axRow):

                counter += 1
                plotEnabled = False

                self._progress.note( msg = "generating subplot #" + str(counter) + ": (" + str(irow) + "," + str(icol) + ") out of " + self._ntotStr + "... ", end = "" )

                if icol==irow and self.plotType.diag.enabled:
                
                    plotEnabled = True
                    requestedPlotType = self.plotType.diag.value
                    if self.plotType.diag.value=="histplot": RequestedPlotClass = DensityPlot

                if icol>irow and self.plotType.upper.enabled:

                    plotEnabled = True
                    requestedPlotType = self.plotType.upper.value
                    RequestedPlotClass = DensityPlot if "contour" in self.plotType.upper.value else LineScatterPlot

                if icol<irow and self.plotType.lower.enabled:

                    plotEnabled = True
                    requestedPlotType = self.plotType.lower.value
                    RequestedPlotClass = DensityPlot if "contour" in self.plotType.lower.value else LineScatterPlot

                if plotEnabled:

                    self.currentFig.subplotGrid[irow][icol] = RequestedPlotClass( plotType = requestedPlotType
                                                                                , dataFrame = self._dfref()
                                                                                , methodName = self._methodName
                                                                                , reportEnabled = False
                                                                                , resetPlot = None
                                                                                )

                    # set properties from the template

                    propsDist = dcopy( vars( getattr(self.layout,requestedPlotType) ) )
                    for key, val in propsDist.items(): setattr(self.currentFig.subplotGrid[irow][icol], key, val)

                    self.currentFig.subplotGrid[irow][icol].rows = self.rows
                    self.currentFig.subplotGrid[irow][icol].xcolumns = self._colnames[icol]

                    if not self.currentFig.subplotGrid[irow][icol]._type.is1d:

                        self.currentFig.subplotGrid[irow][icol].ycolumns = self._colnames[irow]

                        if "contourf"==self.currentFig.subplotGrid[irow][icol]._type.name:
                            if not self._cEnabled:
                                self.currentFig.subplotGrid[irow][icol].contourf.kws.cmap = None
                                # adding colors only makes a unicolor square subplot
                                #if "colors" not in vars(self.currentFig.subplotGrid[irow][icol].contourf.kws).keys():
                                #    self.currentFig.subplotGrid[irow][icol].contourf.kws.colors = "k"

                        elif "contour" in self.currentFig.subplotGrid[irow][icol]._type.name:
                            if not self._cEnabled:
                                self.currentFig.subplotGrid[irow][icol].contour.kws.cmap = None
                                #if "colors" not in vars(self.currentFig.subplotGrid[irow][icol].contour.kws).keys():
                                #if "colors" in vars(self.layput.contour.contour.kws).keys():
                                #    self.currentFig.subplotGrid[irow][icol].contour.kws.colors = "k"

                        else:

                            self.currentFig.subplotGrid[irow][icol].ccolumns = self.ccolumn

                    # make plot

                    plt.sca(ax)
                    if icol!=irow:
                        self.currentFig.subplotGrid[irow][icol].make()

                    # set the ticks and labels

                    if self.plotType.lower.enabled:
                        if icol > 0:
                            ax.set_ylabel("")
                        if irow < nrowMinusOne:
                            ax.set_xlabel("")
                    elif self.plotType.upper.enabled:
                        if icol > irow:
                            ax.set_ylabel("")
                            ax.set_xlabel("")

                # report progress

                if self._reportEnabled:
                    self._progress.timer.toc()
                    if irow==icol:
                        delta = avgDiagTime
                    else:
                        delta = self._progress.timer.delta
                    print( "done in " + str(np.round(delta,6)) + " seconds." )

        ############################################################################################################################

        plt.subplots_adjust(hspace=0.15, wspace=0.15)

        # add colorbar

        self._addColorBar()

        ############################################################################################################################

        if reself: return self

    ################################################################################################################################
    #### _addColorBar
    ################################################################################################################################

    def _addColorBar(self):

        if self._cEnabled and self.colorbar.enabled:

            self._progress.note( msg = "generating colorbar... ", end = "" )

            mappable = None
            isDensityMappable = False

            # search for line / scatter mappable

            if mappable is None and self.plotType.upper.enabled:

                if mappable is None and "scatter" in self.plotType.upper.value.lower():
                    mappable = self.currentFig.subplotGrid[0][1].currentFig.scatterList[0]

                if mappable is None and "line"==self.plotType.upper.value and self.layout.line.lineCollection.enabled:
                    mappable = self.currentFig.subplotGrid[0][1].currentFig.lineCollectionList[0]

                if mappable is None and "lineScatter"==self.plotType.upper.value and self.layout.lineScatter.lineCollection.enabled:
                    mappable = self.currentFig.subplotGrid[0][1].currentFig.lineCollectionList[0]

            if mappable is None and self.plotType.lower.enabled:

                if mappable is None and "scatter" in self.plotType.lower.value.lower():
                    mappable = self.currentFig.subplotGrid[1][0].currentFig.scatterList[0]

                if mappable is None and "line"==self.plotType.upper.value and self.layout.line.lineCollection.enabled:
                    mappable = self.currentFig.subplotGrid[1][0].currentFig.lineCollectionList[0]

                if mappable is None and "lineScatter"==self.plotType.upper.value and self.layout.lineScatter.lineCollection.enabled:
                    mappable = self.currentFig.subplotGrid[1][0].currentFig.lineCollectionList[0]

            # search for contour mappable

            if mappable is None and self.plotType.upper.enabled:

                if mappable is None and "contour"==self.plotType.upper.value.lower():
                    mappable = self.currentFig.subplotGrid[0][1].currentFig.contourList[0]
                    isDensityMappable = True

                if mappable is None and "contourf" in self.plotType.upper.value.lower():
                    mappable = self.currentFig.subplotGrid[0][1].currentFig.contourfList[0]
                    isDensityMappable = True

            if mappable is None and self.plotType.lower.enabled:

                if mappable is None and "contour"==self.plotType.lower.value.lower():
                    mappable = self.currentFig.subplotGrid[1][0].currentFig.contourList[0]
                    isDensityMappable = True

                if mappable is None and "contourf" in self.plotType.lower.value.lower():
                    mappable = self.currentFig.subplotGrid[1][0].currentFig.contourfList[0]
                    isDensityMappable = True

            # add colorbar

            if mappable is not None:

                self.colorbar.kws.mappable = mappable
                self.colorbar.kws.ax = self.currentFig.pairgrid.axes
                self.currentFig.colorbar = self.currentFig.figure.colorbar( **vars(self.colorbar.kws) )

                if isDensityMappable:
                    colorBarLabel = "Density"
                elif self._ccolindexlen==0:
                    colorBarLabel = "Count from the Series Start"
                else:
                    colorBarLabel = self._ccolname[0]
                self.currentFig.colorbar.set_label(colorBarLabel)

                self.currentFig.figure.set_figwidth(self.currentFig.figure.get_figwidth()*1.3333)

                self._progress.timer.toc()
                delta = self._progress.timer.delta
                msg = "done in " + str(np.round(delta,6)) + " seconds."

            else:

                msg = "No mappable were found in the current GridPlot. skipping..."

            if self._reportEnabled: print( msg, end = newline )

    ################################################################################################################################
    #### hide
    ################################################################################################################################

    def hide(self,part="all"):
        """

        Hides the requested part of the grid plot.

        **Parameters**

            part

                A string with the following possible values:

                    "lower"
                        hides the lower triangle of the grid plot.

                    "upper"
                        hides the upper triangle of the grid plot.

                    "diag"
                        hides the diagonal of the grid plot.

                    "all"
                        hides all grid plots and the colorbar.

                    "colorbar"
                        hides the colorbar of the grid plot.

                    The string can also be a mix of the above keywords,
                    separated by the ``+`` sign or some other delimiter.
                    For example, ``"lower+upper+colorbar"``

        **Returns**

            None

        """

        allhidden = "all" in part
        if allhidden or ("upper" in part): self.currentFig.pairgrid.map_upper(_hide_current_axis)
        if allhidden or ("lower" in part):
            self.currentFig.pairgrid.map_lower(_hide_current_axis)
            colnames = pm.dfutils.getColNamesIndex(self._dfref().columns,self.columns)[0]
            if not allhidden:
                for i in range(len(colnames)):
                    ax=self.currentFig.pairgrid.axes[:][i][i]
                    ax.set_xlabel(colnames[i])
                    ax.set_ylabel(colnames[i])
                    ax.xaxis.set_tick_params(which="both", labelbottom=True)
                    ax.yaxis.set_tick_params(which="both", labelbottom=True)
        if allhidden or ("diag" in part): self.currentFig.pairgrid.map_diag(_hide_current_axis)
        if allhidden or ("colorbar" in part): self.currentFig.colorbar.remove()
        #self.currentFig.figure.set_figwidth(self.currentFig.figure.get_figwidth()*1.3333)

    ################################################################################################################################
    #### show
    ################################################################################################################################

    def show(self,part="all"):
        """

        Shows the requested part of the grid plot.

            **Parameters**

                part
                    a string with the following possible values:

                        "lower"
                            shows the lower triangle of the grid plot.

                        "upper"
                            shows the lower triangle of the grid plot.

                        "diag"
                            shows the diagonal of the grid plot.

                        "all"
                            shows all grid plots.

                        "colorbar"
                            shows the colorbar of the grid plot.

                        The string can also be a mix of the above keywords,
                        separated by the ``+`` sign or some other delimiter.
                        For example, ``"lower+upper+colorbar"``

            **Returns**

                None

        """

        allshown = "all" in part
        if allshown or ("upper" in part): self.currentFig.pairgrid.map_upper(_show_current_axis)
        if allshown or ("lower" in part):
            self.currentFig.pairgrid.map_lower(_show_current_axis)
            colnames = pm.dfutils.getColNamesIndex(self._dfref().columns,self.columns)[0]
            colnamesLength=len(colnames)
            for i in range(colnamesLength):
                ax=self.currentFig.pairgrid.axes[:][i][i]
                if colnamesLength-1>i>0:
                    ax.set_xlabel("")
                    ax.set_ylabel("")
                    ax.xaxis.set_tick_params(which="both", labelbottom=False)
                    ax.yaxis.set_tick_params(which="both", labelbottom=False)
                elif i==0:
                    ax.xaxis.set_tick_params(which="both", labelbottom=False)
                    ax.set_xlabel("")
                else:
                    ax.yaxis.set_tick_params(which="both", labelbottom=False)
                    ax.set_ylabel("")
        if allshown or ("diag" in part): self.currentFig.pairgrid.map_diag(_show_current_axis)
        if allshown or ("colorbar" in part): self._addColorBar()

    ################################################################################################################################
    #### addTarget
    ################################################################################################################################

    def addTarget   ( self
                    , value : tp.Union[ np.ndarray , tp.List[float] , tp.Tuple[float] ] = None
                    ):
        """

        Call the target callable associated with each element of
        currentFig.subplotGrid of the GridPlot object to add target
        points and/or lines to all plots of the ``GridPlot`` figure.
        This method is supposed to be called only after a grid plot
        has been generated.

        **Parameters**

            value (optional)

                A numpy array, or list, or tuple whose length equals the
                number of columns/rows of the grid plot, each element of
                which represents the target value associated with the
                corresponding variable on the ``x`` axis of the plot,
                from the left to the right.

                Alternatively, ``value`` can be a string with the following 
                possible value:

                    "mode"

                        The variable values corresponding to the mode of the 
                        "SampleLogFunc" column of the input dataFrame will 
                        be used. If no "SampleLogFunc" columns name exists 
                        in the input dataFrame, an exception will be raised.

                        This is the **default value** for the input variable 
                        ``value``.

                    "mean"

                        The variable values corresponding to the mode of the 
                        "SampleLogFunc" column of the input dataFrame will 
                        be used. If no "SampleLogFunc" columns name exists 
                        in the input dataFrame, an exception will be raised.

                If not provided, the default will be set to the state
                (coordinates) of the mode of the "SampleLogFunc" column 
                in the input dataFrame to the object. This would work only
                when the input dataFrame is the contents of a ParaMonte
                output chain or sample file.

        **Returns**

            None. However, this method causes side-effects by manipulating
            the existing attributes of the ``Target`` objects in 
            ``currentFig.subplotGrid`` of the GridPlot object.

        """

        if value is None: value="mode"

        if isinstance(value, str):
            if value=="mode":
                from _dfutils import getMaxLogFunc
                maxLogFunc = getMaxLogFunc(dataFrame=self._dfref().iloc[self.rows,:])
                value = maxLogFunc.dfrow[self._colindex]
            elif value=="mean":
                value = self._dfref()[self._colnames].iloc[self.rows,:].mean()
            elif value=="median":
                value = self._dfref()[self._colnames].iloc[self.rows,:].median()
        else:
            from collections.abc import Iterable
            if not isinstance(value, Iterable):
                raise Exception ( newline
                                + "The input argument ``value`` must be either a string or an array of target value." + newline
                                + "Here is the help information for addTarget():"
                                + newline
                                + self.addTarget.__doc__
                                )

        counter = 0
        for irow, plotObjectRow in enumerate(self.currentFig.subplotGrid[:]):

            for icol, plotObject in enumerate(plotObjectRow):

                counter += 0
                self._progress.note( msg = "generating target #" + str(counter) + ": (" + str(irow) + "," + str(icol) + ") out of " + self._ntotStr + "... ", end = "" )

                #import matplotlib.pyplot as plt
                #plt.sca(self.currentFig.pairgrid.axes[irow][icol])
                if icol>irow and self.plotType.upper.enabled:
                    plotObject.target( value = [ value[icol], value[irow] ], axes = self.currentFig.pairgrid.axes[irow][icol] )
                elif icol<irow and self.plotType.lower.enabled:
                    plotObject.target( value = [ value[icol], value[irow] ], axes = self.currentFig.pairgrid.axes[irow][icol] )
                elif icol==irow and self.plotType.diag.enabled:
                    plotObject.target( value = [ value[icol], value[irow] ], axes = self.currentFig.pairgrid.axes[irow][icol] )

                # report progress

                if self._reportEnabled:
                    self._progress.timer.toc()
                    delta = self._progress.timer.delta
                    print( "done in " + str(np.round(delta,6)) + " seconds.", end = newline )

    ################################################################################################################################
    #### _getDocString
    ################################################################################################################################

    def _getDocString(self):
        docString   = newline \
                    + "Here is the help information on the GridPlot class:" + newline \
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

####################################################################################################################################

def _hide_current_axis(*args, **kwds):
    from matplotlib import pyplot as plt
    plt.gca().set_visible(False)

####################################################################################################################################

def _show_current_axis(*args, **kwds):
    from matplotlib import pyplot as plt
    plt.gca().set_visible(True)

####################################################################################################################################
