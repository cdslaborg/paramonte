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

import _message as err
import _dfutils as dfutils
import _pmutils as pmutils
from paramonte.vis.LineScatterPlot import LineScatterPlot

Struct = pmutils.Struct
newline = pmutils.newline

####################################################################################################################################
#### AutoCorr class
####################################################################################################################################

class AutoCorr:
    """

    This is the base class for generating object of type ``AutoCorr``.  
    Upon construction, it will provide methods to compute and plot the 
    autocorrelations of the selected columns of the input dataFrame.

        **Parameters**

            dataFrame

                A Pandas dataframe based upon the selected columns 
                of which the autocorrelations will be computed.

            columns (optional)

                A list that determines the columns of the input dataFrame 
                to be used in the computation of the autocorrelations. 
                It can have three forms:

                    1.  A list of column indices from the input dataFrame.
                    2.  A list of column names from dataFrame.columns.
                    3.  A ``range(start,stop,step)`` of column indices.

                Examples:

                    1.  ``columns = [0,1,4,3]``
                    2.  ``columns = ["SampleLogFunc","SampleVariable1"]``
                    3.  ``columns = range(17,7,-2)``

                The default behavior includes all columns of the dataFrame.

        **Attributes**

            All of the parameters described above, except dataFrame.

                A weak reference to the dataFrame will be implicitly 
                stored in the object. 

            df

                A pandas dataframe containing the computed autocorrelations.

            rows

                A list of the dataFrame indices that determines the rows of 
                the input dataFrame to be used in the computation of the 
                autocorrelations. It can be either:

                    1.  A ``range(start,stop,step)``, or, 
                    2.  A list of row indices from ``dataFrame.index``.

                Examples:

                    1.  ``rows = range(17,7,-2)``
                    2.  ``rows = [i for i in range(7,17)]``

                The default behavior includes all rows of the dataFrame.

            plot

                A structure containing the following plotting tools:

                    heatmap

                        A callable object of class HeatMap which will 
                        enable plotting of the correlation matrix.

        **Returns**

            self

                an object of type class CorMat.

    ---------------------------------------------------------------------------
    """

    ################################################################################################################################
    #### __init__
    ################################################################################################################################

    def __init__( self
                , dataFrame     # : tp.Optional[ pd.DataFrame ]
                , columns       = None # : tp.Optional[ tp.Union[ range , tp.List[int] , tp.List[str] ] ] = None
                , methodName    = "ParaMonte" # : tp.Optional[ str ] = "ParaMonte"
                , reportEnabled = True # : tp.Optional[ bool ] = True
                ):

        self.df = None
        self.rows = None
        self.columns = columns

        self._dfref = None if dataFrame is None else wref.ref(dataFrame)

        self._methodName = methodName
        self._reportEnabled = reportEnabled
        self._progress = pmutils.Progress   ( msg = None
                                            , methodName = methodName
                                            , reportEnabled = reportEnabled
                                            , end = "\n"
                                            )

    ################################################################################################################################
    #### __call__
    ################################################################################################################################

    def __call__( self
                , reself : tp.Optional[ bool ] = False
                , **kwargs
                ):
        """

        Call the ``get()`` method of the current instance of the class.

            **Parameters**

                Any arguments that can be passed to the ``get()`` 
                method of the object.

            **Returns**

                Any return value from the ``get()`` method of the object.

        """

        return self.get(reself, **kwargs)

    ################################################################################################################################
    #### get
    ################################################################################################################################

    def get ( self
            , reself : tp.Optional[ bool ] = False
            , **kwargs
            ):
        """

        Compute the autocorrelations of the selected columns 
        of the input dataframe to the object's constructor.

            **Parameters**

                reself

                    A logical variable. When ``True``, an instance of 
                    the object will be returned to the calling routine 
                    upon exit. The default value is ``False``.

            **Returns**

                The object if ``reself = True`` otherwise, ``None``.

                **NOTE**

                This method causes side-effects by manipulating
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

        ############################################################################################################################
        #### check columns presence
        ############################################################################################################################

        colnames, colindex = dfutils.getColNamesIndex(self._dfref().columns,self.columns)

        ############################################################################################################################
        #### check rows presence. This must be checked here, because it depends on the integrity of the in input dataFrame.
        ############################################################################################################################

        if self.rows is None: self.rows = range(len(self._dfref().index))
        #rownames = self._dfref().index[self.rows]

        ############################################################################################################################
        #### compute the autocorrelations
        ############################################################################################################################

        self._nvar = len(colnames)
        self._nlag = len(self.rows)
        acf = np.zeros((self._nvar,self._nlag))

        try:

            from scipy.signal import correlate
            for cnt, ivar in enumerate(colindex):
                xdata = self._dfref().iloc[self.rows,ivar].values.flatten() - np.mean(self._dfref().iloc[self.rows,ivar].values.flatten())
                acf[cnt] = correlate( xdata
                                    , xdata
                                    , mode = "full"
                                    )[self._nlag-1:2*self._nlag]
                acf[cnt] = acf[cnt] / acf[cnt,0]

        except:

            if self._reportEnabled:
                err.warn( msg   = "Failed to compute the Autocorrelation function of the input dataFrame." + newline
                                + "This could likely be due to an error in importing the scipy Python library." + newline
                                + "Please make sure you have the scipy library properly installed on your system." + newline
                                + "You can do so by typing the following command on your Anaconda3 or Bash command prompt:" + newline
                                + newline
                                + "    pip install --user --upgrade scipy"
                        , methodName = self._methodName
                        , marginTop = 1
                        , marginBot = 1
                        )

        self.df = pd.DataFrame(np.transpose(acf))

        ############################################################################################################################
        #### specify columns/index names
        ############################################################################################################################

        colnames = [ "ACF_"+colnames[i] for i in range(len(colnames)) ]
        self.df.columns = colnames

        ############################################################################################################################
        #### add SampleLogFunc to df for plot coloring, if exists
        ############################################################################################################################

        ccolumns = ()
        if "SampleLogFunc" in self._dfref().columns:
            ccolumns = "SampleLogFunc"
            self.df.insert  ( loc = 0
                            , column = ccolumns
                            , value = self._dfref()[[ccolumns]].values.flatten()
                            , allow_duplicates = True
                            )

        ############################################################################################################################
        #### add lags to df
        ############################################################################################################################

        self.df.insert  ( loc = 0
                        , column = "Lag"
                        , value = [ i for i in self.df.index ]
                        )

        ############################################################################################################################
        #### graphics
        ############################################################################################################################

        self._plotTypeList =    [ "line"
                                , "scatter"
                                , "lineScatter"
                                ]

        self._progress.note( msg = "adding the autocrrelation graphics tools... ", end = newline, pre = True )
        self.plot = Struct()
        self._resetPlot(resetType="hard")

        self.plot.reset = self._resetPlot

        ############################################################################################################################

        if reself: return self

    ################################################################################################################################
    #### _resetPlot
    ################################################################################################################################

    def _resetPlot  ( self
                    , resetType = "soft"
                    , plotNames = "all"
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

                plotNames (optional)

                    An optional string value or list of string values representing 
                    the names of plots to reset. If no value is provided, 
                    then all plots will be reset.

            **Returns**

                None

            **Example**

                .. code-block:: python

                    reset("hard")                    # regenerate all plots from scratch
                    reset("hard","line")             # regenerate line plot from scratch
                    reset("hard",["line","scatter"]) # regenerate line and scatter plots

        """

        if isinstance(plotNames, str):
            plotTypeLower = plotNames.lower()
            if plotTypeLower=="all":
                requestedPlotTypeList = self._plotTypeList
            elif plotNames in self._plotTypeList:
                requestedPlotTypeList = [plotNames]
            else:
                self._reportWrongPlotName(plotNames)
        elif isinstance(plotNames, list):
            for plotName in plotNames:
                if plotName not in self._plotTypeList: self._reportWrongPlotName(plotName)
        else:
            self._reportWrongPlotName("a none-string none-list object.")

        resetTypeIsHard = None
        if isinstance(resetType, str):
            resetTypeIsHard = resetType.lower()=="hard"
        else:
            err.abort   ( msg   = "The input argument resetType must be a string representing" + newline
                                + "the type of the reset to be performed on the plots." + newline
                                + "A list of possible plots includes: \"hard\", \"soft\"" + newline
                                + "Here is the help for the ``reset()`` method: " + newline
                                + newline
                                + self._resetPlot.__doc__
                        , marginTop = 1
                        , marginBot = 1
                        , methodName = self._methodName
                        )

        ############################################################################################################################
        #### reset plots
        ############################################################################################################################

        for requestedPlotType in requestedPlotTypeList:

            plotObject = None
            requestedPlotTypeLower = requestedPlotType.lower()

            isLine      = "line"        in requestedPlotTypeLower
            isScatter   = "scatter"     in requestedPlotTypeLower

            isLineScatterPlot = isLine or isScatter

            if not resetTypeIsHard:
                plotComponent = getattr(self, "plot")
                plotObject = getattr(plotComponent, requestedPlotType)
                plotObject._reset()

            ########################################################################################################################
            #### reset line / scatter
            ########################################################################################################################

            if isLineScatterPlot:

                if resetTypeIsHard:
                    plotObject = LineScatterPlot( plotType = requestedPlotType
                                                , dataFrame = self.df
                                                , methodName = self._methodName
                                                , reportEnabled = self._reportEnabled
                                                , resetPlot = self._resetPlot
                                                )

                plotObject.xcolumns = "Lag"
                plotObject.ycolumns = self.df.columns[2:]

                plotObject._xlabel = "AutoCorrelation Lag"
                plotObject._ylabel = "AutoCorrelation Function (ACF) Value"

                plotObject.ccolumns = None

                plotObject.colorbar.kws.extend = "neither"
                plotObject.colorbar.kws.orientation = "vertical"
                plotObject.colorbar.kws.spacing = "uniform"

                plotObject._xlimit = [1, None]
                plotObject._ylimit = [None, 1]
                plotObject._xscale = "log"

                plotObject.legend.enabled = True

                if isLine:
                    if isScatter:
                        plotObject.lineCollection.enabled = False
                        plotObject.plot.enabled = True
                        plotObject.plot.kws.alpha = 0.2
                        plotObject.plot.kws.color = "grey"
                        plotObject.plot.kws.linewidth = 0.75
                    else:
                        plotObject.plot.enabled = True
                        plotObject.plot.kws.linewidth = 1
                        plotObject.lineCollection.enabled = False

            ########################################################################################################################

            if plotObject is not None: setattr(self.plot, requestedPlotType, plotObject)

    ################################################################################################################################
    #### _reportWrongPlotName
    ################################################################################################################################

    def _reportWrongPlotName( self
                            , plotNames
                            ):

        err.abort   ( msg   = "The input argument plotNames must be a string representing" + newline
                            + "the name of a plot belonging to the AutoCorr class or," + newline
                            + "a list of such plot names. You have entered: " + plotNames + newline
                            + "Possible plots are: " + newline
                            + newline
                            + newline.join(self._plotTypeList) + newline
                            + newline
                            + "Here is the help for the ``reset()`` method: " + newline
                            + newline
                            + self._resetPlot.__doc__
                    , marginTop = 1
                    , marginBot = 1
                    , methodName = self._methodName
                    )

    ################################################################################################################################

