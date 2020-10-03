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

import os
import numpy as np
import typing as tp
import weakref as wref
import _paramonte as pm
import _pmutils as pmutils

Struct = pmutils.Struct
newline = pmutils.newline

####################################################################################################################################
#### BasePlot class
####################################################################################################################################

class BasePlot:
    """

    This is the class for generating instances of
    basic plots with minimally one (X)-axis. It serves as the
    superclass for a wide variety of other multi-axes ParaMonte plots.

        **Parameters**

            plotType

                A string indicating the name of the plot that is to be constructed.

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

            _dfref

                A weak reference to the input dataFrame whose data is used to 
                generate plots.

            rows

                A numeric vector that determines the rows of dataFrame
                to be visualized. It can be either:

                    1.  a ``range(start,stop,step)``, or,
                    2.  a list of row indices in ``dataFrame.index``.

                Examples:

                    1.  ``rows = range(17,7,-2)``
                    2.  ``rows = [i for i in range(7,17)]``

                If not provided, the default includes all rows of the dataframe.

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
                        keyword arguments to the ``figure()`` function of 
                        the matplotlib library.

                        Example usage:

                            .. code-block:: python

                                figure.kws.facecolor = "w"

                        **NOTE**

                        If a desired property is missing among the ``kws`` 
                        attributes, simply add the field and its value to 
                        the component.

            legend (may not exist for some types of plots)

                A structure with two attributes:

                    enabled

                        A boolean indicating whether a call to the ``legend()``
                        function of the matplotlib library should be made or not.
                        If a call is made, a new figure will be generated.
                        Otherwise, the current active figure will be used.

                    kws

                        A structure whose components are directly passed as 
                        keyword arguments to the ``legend()`` function of 
                        the matplotlib library.

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

        **Returns**

            An object of ``BasePlot`` class.

    ---------------------------------------------------------------------------
    """

    ################################################################################################################################
    #### __init__
    ################################################################################################################################

    def __init__( self
                , plotType
                , dataFrame = None
                , methodName = "ParaMonte"
                , reportEnabled = True
                , resetPlot = None
                ):

            self._dfref = None if dataFrame is None else wref.ref(dataFrame)

            self._type = Struct()
            self._type.name = plotType

            plotTypeLower = plotType.lower()

            self._type.is3d         = "3"           in plotTypeLower
            self._type.isLine       = "line"        in plotTypeLower
            self._type.isScatter    = "scatter"     in plotTypeLower
            self._type.isHeatmap    = "heatmap"     in plotTypeLower
            self._type.isKdeplot1   = "kdeplot1"    in plotTypeLower
            self._type.isKdeplot2   = "kdeplot2"    in plotTypeLower
            self._type.isHistplot   = "histplot"    in plotTypeLower
            self._type.isJointplot  = "jointplot"   in plotTypeLower
            self._type.isEllipsoid  = "covmat"      in plotTypeLower
            self._type.isEllipsoid  = "cormat"      in plotTypeLower
            self._type.isGridPlot   = "grid"        in plotTypeLower
            self._type.isContour3   = "contour3"    in plotTypeLower
            self._type.isContourf   = "contourf"    in plotTypeLower
            self._type.isContour    = "contour"     in plotTypeLower and not (self._type.isContour3 or self._type.isContourf)

            self._type.is1d         = self._type.isKdeplot1 or self._type.isHistplot
            self._type.is2d         = not (self._type.isGridPlot or self._type.is1d or self._type.is3d)

            self._type.isDiffusionPlot = self._type.isContour or self._type.isContourf or self._type.isContour3

            self._methodName = methodName
            self._reportEnabled = reportEnabled

            self._progress = pm.utils.Progress  ( msg = "creating a " + self._type.name + " plot object from scratch... "
                                                , methodName = methodName
                                                , reportEnabled = reportEnabled
                                                , end = ""
                                                )

            self._reset()

            if resetPlot is None:
                self._resetPlot = self._reset
            else:
                self._resetPlot = resetPlot

    ################################################################################################################################
    #### reset
    ################################################################################################################################

    def reset   ( self
                , resetType = "soft"
                , **kwargs
                ): 
        """

        Reset the properties of the plot to the original default settings.
        Use this method when you change many attributes of the plot and
        you want to clean up and go back to the default settings.

            **Parameters**

                resetType (optional)

                    An optional string with possible value of ``"hard"``.
                    If provided, the plot object will be regenerated from scratch.
                    This includes reading the original data frame again and resetting
                    everything. If not provided, then only the plot settings will be
                    reset without reseting the dataFrame.

            **Returns**

                None

            **Example**

                .. code-block:: python
                    :linenos:

                    reset()         # reset the plot to the default settings
                    reset("soft")   # reset the plot to the default settings
                    reset("hard")   # regenerate the plot from scratch

        """

        try:
            self._resetPlot(resetType = resetType, plotNames = self._type.name) # call external reset function
        except:
            self._resetPlot(**kwargs) # call _reset(self)

    ################################################################################################################################
    #### _reset
    ################################################################################################################################

    def _reset(self):

        self.rows = None
        self._fnameOld = None # used in savefig() to generate unique random file names.

        self.set = Struct()
        self.set.enabled = True
        self.set.kws = Struct()
        if self._type.isKdeplot2 or self._type.isJointplot:
            self.set.kws.style = "ticks"
        elif self._type.is3d:
            self.set.kws.style = "white"
        else:
            self.set.kws.style = "darkgrid"

        self.figure = Struct()
        self.figure.enabled = True
        self.figure.kws = Struct()

        self.currentFig = Struct()

        if self._type.is3d:
            self.axes3d = Struct()
            self.axes3d.kws = Struct()
            self.axes3d.kws.alpha = 1
            self.axes3d.kws.visible = True
            #self.axes3d.enabled = self.figure.enabled
        else:
            self.axes = Struct()
            self.axes.kws = Struct()
            self.axes.kws.alpha = 1
            self.axes.kws.visible = True
            #self.axes.enabled = self.figure.enabled

        if not (self._type.isGridPlot or self._type.isHeatmap or self._type.isDiffusionPlot):
            self.legend = Struct()
            self.legend.enabled = False
            self.legend.kws = Struct()

    ################################################################################################################################
    #### savefig
    ################################################################################################################################

    def savefig ( self
                , reself = False
                , **savefig_kws
                ):
        """

        Export the current figure to external file.

            **Parameters**

                reself

                    A logical variable. If ``True``, the path to 
                    the output generated file will be returned.
                    The default value is ``False``.

                **savefig_kws (optional)

                The set of input arguments to the ``savefig()`` function 
                of the matplotlib library, including the ``fname``. If 
                the file name is not provided, a unique random filename 
                will be generated and used to save the figure to an 
                output file.

                Example:

                .. code-block:: python

                    savefig(fname = "gridplot.png")

        **Returns**

            None

        **Example**

            .. code-block:: python

                savefig() # use a unique random filename to output the plot.
                savefig(fname = "thisPlot.png")
                savefig(pad_inches = 0.0, bbox_inches = "tight")

        """

        if "fname" not in savefig_kws.keys() or (isinstance(self._fnameOld,str) and isinstance(savefig_kws["fname"],str) and savefig_kws["fname"]==self._fnameOld):
            savefig_kws["fname"] = os.path.join( os.getcwd(), pm.utils.getRandomFilePrefix(prefix=self._type.name+"_") + ".png" )
            self._fnameOld = savefig_kws["fname"]

        import matplotlib.pyplot as plt # xxx how is this figure activated?
        self._progress.note( msg = "saving the plot to file: \"" + savefig_kws["fname"] + "\"", end = "\n", pre = True )
        plt.savefig(**savefig_kws)
        self._progress.note( end = "\n", pre = True )

    ################################################################################################################################
    #### getLinSpace
    ################################################################################################################################

    def getLinSpace ( self
                    , skip = None
                    , lowerLim  : tp.Optional[ np.int32 ] = 1
                    , upperLim  : tp.Optional[ np.int32 ] = None
                    ):
        """

        Generate linearly-spaced **unique** integer numbers
        between the input lowerLim and upperLim. These numbers 
        can be used as the row indices in the plots.

        **Parameters**

            skip (optional)

                The linear spacing between the generated points.
                If ``skip`` is specified as input, any input value 
                for ``npoint`` will be ignored.
                The default value is ``None``.

            lowerLim (optional)

                The natural (non-logarithmic) lower limit of the
                generated linearly-spaced integer numbers.
                If not provided, the default value is ``1``.

            upperLim (optional)

                The natural (non-logarithmic) upper limit of the
                generated linearly-spaced integer numbers.
                If not provided, the default value is the maximum
                of the number of the rows of the input dataframe
                to the ``BasePlot`` constructor.

        **Returns**

            A set of unique linearly-spaced integer numbers.

        **Example**

            .. code-block:: python

                rows = getLinSpace(3, 1, 10000)

        """

        if upperLim is None and self._dfref is not None:
            upperLim = len(self._dfref().index)
        elif upperLim is None and self.rows is not None:
            upperLim = len(self.rows[-1])
        else:
            pm.abort( msg   = "Either the input argument ``upperLim`` or the dataFrame of the plot object" + newline
                            + "must be given in order to generate a linear range of values. " + newline
                            + "Here is the documentation of getLinSpace():" + newline
                            + newline
                            + self.getLinSpace.__doc__
                            + newline
                    , marginTop = 1
                    , marginBot = 1
                    , methodName = self._methodName
                    )
        if skip is None: 
            skip = (upperLim - lowerLim) // 256
            if skip<1: skip = 1

        return range(lowerLim,upperLim,skip)

    ################################################################################################################################
    #### getLogLinSpace
    ################################################################################################################################

    def getLogLinSpace  ( self
                        , base      : tp.Optional[ np.float64 ] = 1.2
                        , logskip   : tp.Optional[ np.int32 ] = 0.2
                        , lowerLim  : tp.Optional[ np.int32 ] = 1
                        , upperLim  : tp.Optional[ np.int32 ] = None
                        ):
        """

        Generate logarithmically-uniformly-spaced **unique** integer 
        numbers between the input lowerLim and upperLim. These numbers 
        are to be used as the row indices in the plots.

        **Parameters**

            base (optional)

                The base of the logarithm used for generating the
                logarithmically-uniform range of integers.
                The default value is ``1.2``.

            logskip (optional)

                The minimum logarithmic space jump between
                the generated log-linearly-spaced integer numbers.
                The default value is ``0.2``.

            lowerLim (optional)

                The natural (non-logarithmic) lower limit of the
                generated log-linearly-spaced integer numbers.
                If not provided, the default value is ``1``.

            upperLim (optional)

                The natural (non-logarithmic) upper limit of the
                generated log-linearly-spaced integer numbers.
                If not provided, the default value is the maximum
                of the number of the rows of the input dataframe
                to the BasePlot constructor.

        **Returns**

            A set of unique log-linearly-spaced integer numbers.

        **Example**

            .. code-block:: python

                rows = getLogLinSpace(1.01, 1, 1, 10000)

        """

        if upperLim is None and self._dfref is not None:
            upperLim = len(self._dfref().index)
        elif upperLim is None and self.rows is not None:
            upperLim = len(self.rows[-1])
        else:
            pm.abort( msg   = "Either the input argument ``upperLim`` or the dataFrame of the plot object" + newline
                            + "must be given in order to generate a log-linear range of values. " + newline
                            + "Here is the documentation of getLogLinSpace():" + newline
                            + newline
                            + self.getLogLinSpace.__doc__
                            + newline
                    , marginTop = 1
                    , marginBot = 1
                    , methodName = self._methodName
                    )

        return pm.utils.getLogIntSpace(base,logskip,lowerLim,upperLim)

    ################################################################################################################################
    #### _constructBasePlot
    ################################################################################################################################

    def _constructBasePlot(self):
        """

        Generate Figure and Axes instances if needed and return nothing.

            **Parameters**

                None

            **Returns**

                None

        """

        end = "\n" if "grid"==self._type.name else ""
        self._progress.note( msg = "making the " + self._type.name + " plot... ", end = end, pre = True )

        import seaborn as sns
        from matplotlib import pyplot as plt
        if self._type.is3d: from mpl_toolkits.mplot3d import Axes3D

        if self.set.enabled:
            sns.set( **vars(self.set.kws) )

        if not (self._type.isJointplot or self._type.isGridPlot): # jointplot and gridplot create their own figure and axes

            if self.figure.enabled:
                self.currentFig.figure = plt.figure( **vars(self.figure.kws) )
                if self._type.is3d:
                    self.currentFig.axes = Axes3D(self.currentFig.figure, **vars(self.axes3d.kws) )
                else:
                    self.currentFig.axes = self.currentFig.figure.gca( **vars(self.axes.kws) )
            else:
                self.currentFig.figure = plt.gcf() # self.currentFig.axes.get_figure()
                self.currentFig.axes = self.currentFig.figure.gca()

    ################################################################################################################################
    #### _finalizeBasePlot
    ################################################################################################################################

    def _finalizeBasePlot(self):
        """

        Set the base Figure properties, such as legends, ...

            **Parameters**

                None

            **Returns**

                None

        """

        ############################################################################################################################
        #### _finalizeBasePlot
        ############################################################################################################################

        if not (self._type.isGridPlot or self._type.isHeatmap or self._type.isDiffusionPlot) and self.legend.enabled:

            labelsNeeded = False
            if not hasattr(self.legend.kws, "labels"):
                labelsNeeded = True
            elif self.legend.kws.labels is None:
                labelsNeeded = True

            if labelsNeeded and hasattr(self.legend,"_labels") and self.legend._labels is not None:
                self.legend.kws.labels = self.legend._labels

            if self._type.isJointplot:
                self.currentFig.legend = self.currentFig.jointplot.ax_joint.legend( **vars(self.legend.kws) )
            else:
                self.currentFig.legend = self.currentFig.axes.legend( **vars(self.legend.kws) )

            self.currentFig.legend.set_zorder(20) # bring legend to the top

        if not self._type.is3d:
            import matplotlib.pyplot as plt
            plt.tight_layout()

        self._progress.note()

    ################################################################################################################################
    #### _checkDataType
    ################################################################################################################################

    def _checkDataType(self):
        """

        Verify the integrity of the dataFrame and return nothing.

            **Parameters**

                None

            **Returns**

                None

        """

        import pandas as pd

        fatalmsg = None
        if self._dfref is None:
            fatalmsg = "It appears that no data has been passed for plotting." + newline
        elif not isinstance(self._dfref,wref.ref):
            fatalmsg = "It appears that you have messed with the " + newline + " internal representation of data in the object." + newline
        elif not isinstance(self._dfref(), pd.DataFrame):
            fatalmsg = ""
        if fatalmsg is not None:
            raise Exception ( fatalmsg
                            + "The input data must be a pandas' dataframe." + newline
                            + "Please pass a dataFrame to the constructor or at" + newline
                            + "the time of calling the object (which is callable with" + newline
                            + "the same input arguments as the object's constructor."
                            )

    ################################################################################################################################
    #### _getDocString
    ################################################################################################################################

    def _getDocString(self):
        docString   = newline \
                    + "Here is the help information on the BasePlot class:" + newline \
                    + newline \
                    + self.__doc__
        return docString

####################################################################################################################################
