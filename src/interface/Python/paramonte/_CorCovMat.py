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

import _message as err
import _dfutils as dfutils
import _pmutils as pmutils
from paramonte.vis.HeatMapPlot import HeatMapPlot

Struct = pmutils.Struct
newline = pmutils.newline

####################################################################################################################################
#### getCorFromCov
####################################################################################################################################

def getCorFromCov(covMat):
    stdVec = np.sqrt(np.diag(covMat))
    corMat = covMat / np.outer(stdVec, stdVec)
    corMat[covMat == 0] = 0
    return corMat

####################################################################################################################################
#### CorCovMat class
####################################################################################################################################

class CorCovMat:

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
                , reself    : tp.Optional[ bool ] = False
                , **kwargs
                ):
        """

        Call the ``get()`` method of the current instance of the class.

            **Parameters**

                reself

                    A logical variable. If ``True``, an instance of the 
                    object will be returned upon exit to the calling routine.
                    The default value is False.

                Also, any attributes of the current instance of the class.

            **Returns**

                The object self if ``reself = True`` otherwise, ``None``.
                However, this method causes side-effects by manipulating 
                the existing attributes of the object.

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

        Compute the correlation / covariance matrix of the selected 
        columns of the input dataframe to the object's constructor.

            **Parameters**

                reself

                    A logical variable. If ``True``, an instance of 
                    the object will be returned  to the calling routine 
                    upon exit. The default value is ``False``.

            **Returns**

                The object self if ``reself = True`` otherwise, ``None``.

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

        if hasattr(self,"method"):
            self._isCorMat = True
            self._matrixType = "correlation"
            if self.method not in ["pearson","kendall","spearman"]:
                raise Exception ( newline
                                + "The requested correlation type must be one of the following string values,\n"
                                + "    pearson  : standard correlation coefficient\n"
                                + "    kendall  : Kendall Tau correlation coefficient\n"
                                + "    spearman : Spearman rank correlation."
                                )
        else:
            self._isCorMat = False
            self._matrixType = "covariance"

        ############################################################################################################################
        #### check columns presence
        ############################################################################################################################

        if self.columns is None:
            colnames = self._dfref().columns
            colindex = range(len(colnames))
        elif all(isinstance(element, str) for element in self.columns):
            colnames = self.columns
            colindex = dfutils.nam2num( self._dfref().columns , self.columns )
        elif all(isinstance(element,int) for element in self.columns):
            colindex = self.columns
            colnames = self._dfref().columns[colindex]
        else:
            raise Exception ( newline
                            + "The input argument 'columns' must be a list whose elements are all\n"
                            + "    1.   string-valued, each representing the name of the column from\n"
                            + "         the input dataframe to the object's constructor, to be included\n"
                            + "         in the " + self._matrixType + " matrix construction, or,\n"
                            + "    2.   integer-valued, each representing the index of the column from\n"
                            + "         the input dataframe to the object's constructor, to be included\n"
                            + "         in the " + self._matrixType + " matrix construction."
                            )

        ############################################################################################################################
        #### check rows presence. This must be checked here, because it depends on the integrity of the in input dataFrame.
        ############################################################################################################################

        if self.rows is None: self.rows = range(len(self._dfref().index))

        ############################################################################################################################
        #### construct the matrix dataframe
        ############################################################################################################################

        if  self._isCorMat:
            self.df = self._dfref().iloc[self.rows,colindex].corr(method=self.method)
        else:
            self.df = self._dfref().iloc[self.rows,colindex].cov()

        ############################################################################################################################
        #### specify columns/index names
        ############################################################################################################################

        self.df.columns = colnames
        self.df.index   = colnames

        ############################################################################################################################
        #### graphics
        ############################################################################################################################

        self._plotTypeList =    [ "heatmap"
                                ]

        self._progress.note( msg = "adding the " + self._matrixType + " graphics tools... ", end = newline, pre = True )
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
                    reset("hard","heatmap")          # regenerate heatmap plot from scratch

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

            isHeatmap  = "heatmap" in requestedPlotTypeLower

            if not resetTypeIsHard:
                plotComponent = getattr(self, "plot")
                plotObject = getattr(plotComponent, requestedPlotType)
                plotObject._reset()

            ########################################################################################################################
            #### reset heatmap
            ########################################################################################################################

            if isHeatmap:

                if resetTypeIsHard:
                    plotObject = HeatMapPlot( plotType = requestedPlotType
                                            , dataFrame = self.df
                                            , methodName = self._methodName
                                            , reportEnabled = self._reportEnabled
                                            , resetPlot = self._resetPlot
                                            )

                if self._isCorMat:

                    plotObject.annotPrecision = 2
                    plotObject.heatmap.kws.cbar_kws =   { "label": self.method.capitalize() + "'s Correlation Strength"
                                                        , "orientation": "vertical"
                                                        , "ticks": np.linspace(-1,1,9)
                                                        }

                    plotObject.heatmap.kws.vmin   = -1
                    plotObject.heatmap.kws.vmax   = +1
                    plotObject.heatmap.kws.center = 0

                else:

                    plotObject.annotPrecision = None
                    plotObject.heatmap.kws.cbar_kws =   { "label": "Covariance Strength"
                                                        , "orientation": "vertical"
                                                        }

            ########################################################################################################################

            if plotObject is not None: setattr(self.plot, requestedPlotType, plotObject)

    ################################################################################################################################
    #### _reportWrongPlotName
    ################################################################################################################################

    def _reportWrongPlotName( self
                            , plotNames
                            ):

        err.abort   ( msg   = "The input argument plotNames must be a string representing" + newline
                            + "the name of a plot belonging to the CorCovMat class or," + newline
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

####################################################################################################################################
#### CorMat class
####################################################################################################################################

class CorMat(CorCovMat):
    """

    This is the class for generating object of type ``CorMat`` which, 
    upon construction, will provide methods to compute and plot the 
    correlation matrix of the selected columns of the input dataFrame.

        **Parameters**

            dataFrame

                A Pandas dataFrame based upon the selected columns of which 
                the correlation matrix will be computed.

            columns (optional)

                optional argument that determines the columns of the input 
                dataFrame to be used in the computation of the correlation 
                matrix. It can have three forms:

                    1.  A list of column indices from the input dataFrame. 
                    2.  A list of column names from dataFrame.columns. 
                    3.  A ``range(start,stop,step)`` of column indices.

                Examples:

                    1.  ``columns = [0,1,4,3]``
                    2.  ``columns = ["SampleLogFunc","SampleVariable1"]``
                    3.  ``columns = range(17,7,-2)``

                The default behavior includes all columns of the dataFrame.

            method (optional)

                A string value representing the method to be used 
                for the computation of correlations:

                    1.  ``'pearson'``    : standard correlation coefficient, 
                    2.  ``'kendall'``    : Kendall Tau correlation coefficient, 
                    3.  ``'spearman'``   : Spearman rank correlation.

                The default value is ``'pearson'``.

        **Attributes**

            All of the parameters described above, except dataFrame.

                A reference to the dataFrame will be implicitly 
                stored in the object.

            df

                A pandas dataframe containing the computed correlation matrix.

            rows

                A list that determines the rows of the input dataFrame to be 
                used in the computation of the correlation matrix. 
                It can be either:

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

                An object of type class ``CorMat``.

    ---------------------------------------------------------------------------
    """

    def __init__( self
                , dataFrame     # : tp.Optional[ pd.DataFrame ]
                , columns       = None # : tp.Optional[ tp.Union[ range , tp.List[int] , tp.List[str] ] ] = None
                , methodName    = "ParaMonte" # : tp.Optional[ str ] = "ParaMonte"
                , reportEnabled = True # : tp.Optional[ bool ] = True
                , method        = "pearson" # : tp.Optional[str] = "pearson"
                ):

        self.method = method
        super().__init__( dataFrame     = dataFrame
                        , columns       = columns
                        , methodName    = methodName
                        , reportEnabled = reportEnabled
                        )

####################################################################################################################################
#### CovMat class
####################################################################################################################################

class CovMat(CorCovMat):
    """

    This is the class for generating object of type ``CovMat`` which, 
    upon construction, will provide methods to compute and plot the 
    covariance matrix of the selected columns of the input dataFrame.

        **Parameters**

            dataFrame

                A Pandas dataframe based upon the selected columns of 
                which the covariance matrix will be computed.

            columns (optional)

                A argument that determines the columns of the input 
                dataFrame to be used in the computation of the covariance 
                matrix. It can have three forms:

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

                A reference to the dataFrame will be implicitly 
                stored in the object.

            df

                A pandas dataframe containing the computed covariance matrix.

            rows

                A list that determines the rows of the input dataFrame to be 
                used in the computation of the covariance matrix. It can be either:

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

                An object of type class ``CovMat``.

    ---------------------------------------------------------------------------
    """

    def __init__( self
                , dataFrame     # : tp.Optional[ pd.DataFrame ]
                , columns       = None # : tp.Optional[ tp.Union[ range , tp.List[int] , tp.List[str] ] ] = None
                , methodName    = "ParaMonte" # : tp.Optional[ str ] = "ParaMonte"
                , reportEnabled = True # : tp.Optional[ bool ] = True
                ):

        super().__init__( dataFrame     = dataFrame
                        , columns       = columns
                        , methodName    = methodName
                        , reportEnabled = reportEnabled
                        )
