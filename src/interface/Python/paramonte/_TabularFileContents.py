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
import pandas as pd
import _paramonte as pm
import _CorCovMat as ccm
from _AutoCorr import AutoCorr
from _dfutils import getMaxLogFunc
from _OutputFileContents import OutputFileContents
from paramonte.vis.LineScatterPlot import LineScatterPlot
from paramonte.vis.DensityPlot import DensityPlot
from paramonte.vis.GridPlot import GridPlot

Struct = pm.Struct
newline = pm.newline

####################################################################################################################################
#### TabularFileContents class
####################################################################################################################################

class TabularFileContents(OutputFileContents):
    """

    This is the **TabularFileContents** class for generating instances
    of the ParaMonte tabular output contents. This class is NOT meant to
    be directly accessed by the ParaMonte library users. It is internally
    used by the ParaMonte library to parse the tabular contents of the 
    output files generated by the ParaMonte sampler routines. For example, 
    the ParaDRAM sampler class makes calls to this class via its 
    ``readSample()`` or ``readChain()`` or ``readMarkovChain()`` 
    or ``readProgress()`` methods to return a list of objects 
    of class ``TabularFileContents``.

        **Parameters**

            file

                The full path to the file containing the sample/chain.

            fileType

                A string containing the type of the file to be parsed.
                Current options include but are not limited to:
                ``sample``, ``chain``, ``markovChain``, ``progress``

            delimiter

                The delimiter used in the sample/chain file, which
                must be provided by the user.

            methodName

                A string representing the name of the ParaMonte sampler used
                to call the constructor of the ``TabularFileContents`` class.

            parseContents

                If set to ``True``, the contents of the file will be parsed and
                stored in a component of the object named ``contents``.
                The default value is ``True``.

            reportEnabled

                A logical input parameter indicating whether the ParaMonte
                automatic guidelines to the standard output should be provided 
                or not. The default value is ``True``.

        **Attributes**

            file

                The full path to the file containing the sample/chain.

            delimiter

                The delimiter used in the sample/chain file, which
                must be provided by the user.

            ndim

                The number of dimensions of the domain of the objective
                function from which the sample has been drawn.

            count

                The number of points (states) in the sample/chain file.
                This is essentially, the number of rows in the file
                minus one (representing the header line).

            plot

                A structure containing the graphics tools for the 
                visualization of the contents of the file.

            df

                If the input file contents is structured in a format that
                could be read as a dataframe, then the contents of the file
                will be stored in the form of a pandas-library DataFrame
                in this property (hence called ``df``).

            contents

                If the input file contents is structured in the form of columns,
                then a property named ``contents`` is also added to the object.
                Each component of contents will named via the header of the file
                and will contain data from the corresponding column of the file.

        **Returns**

            tabularFileContents

                An object of class ``TabularFileContents``.

    ----------------------------------------------------------------------
    """

    def __init__( self
                , file
                , fileType
                , delimiter
                , methodName
                , parseContents = True
                , reportEnabled = True
                ):

        super().__init__(file, methodName, reportEnabled)


        markovChainRequested = fileType=="markovChain"
        self._isProgressFile = "progress"==fileType

        self._sampleLogFuncColName = "" if self._isProgressFile else "SampleLogFunc"

        #if "sample"==fileType:
        #    fileSuffix = "sample"
        #elif fileType=="chain" or markovChainRequested:
        #    fileSuffix = "chain"
        #elif self._isProgressFile:
        #    fileSuffix = "progress"
        #else:
        #    pm.abort( msg   = "Internal error occurred. The input fileType is not recognized.\n"
        #                    + "Please report this error at:\n\n"
        #                    + "    " + pm.website.github.issues.url
        #            , methodName = self._methodName
        #            , marginTop = 1
        #            , marginBot = 1
        #            )

        if fileType!="sample" and fileType!="chain" and not (self._isProgressFile or markovChainRequested):
            pm.abort( msg   = "Internal error occurred. The input fileType is not recognized.\n"
                            + "Please report this error at:\n\n"
                            + "    " + pm.website.github.issues.url
                    , methodName = self._methodName
                    , marginTop = 1
                    , marginBot = 1
                    )

        ############################################################################################################################
        #### data
        ############################################################################################################################

        self.delimiter = delimiter

        self.df = pd.read_csv   ( self.file
                                , delimiter = self.delimiter
                                , header = 0
                                )

        if self._isProgressFile:
            self._offset = -1
        else:
            self._offset = list(self.df.columns).index(self._sampleLogFuncColName) + 1 # index of the first variable
            self.ndim = len(self.df.columns) - self._offset

        self.count  = len(self.df.iloc[:,1])
        self.ncol  = len(self.df.iloc[1,:])

        if markovChainRequested:

            CumSumWeight = np.cumsum(self.df.iloc[:,self._offset-2].values, dtype=np.int32)
            if CumSumWeight[-1] != self.count: # it is indeed a compact chain
                #dfMarkov = pd.DataFrame( columns=list(self.df.columns), index=list(range(CumSumWeight[-1])) )
                dfMarkov = np.zeros( (CumSumWeight[-1] , self.ndim+self._offset) )
                istart = 0
                for i in range(self.count):
                    iend = CumSumWeight[i]
                    #dfMarkov.iloc[istart:iend,:] = self.df.iloc[i].values
                    dfMarkov[istart:iend,:] = self.df.iloc[i].values
                    istart = iend
                columns = self.df.columns
                self.df = pd.DataFrame(dfMarkov)
                self.count = len(self.df.iloc[:,1])
                self.df.columns = columns

        self._progress.note()
        if not self._isProgressFile:
            self._progress.note( msg = "ndim = " + str(self.ndim) + ", count = " + str(self.count), end = newline, pre = True )

        # set dynamic properties

        if parseContents:
            self._progress.note( msg = "parsing file contents... ", end = newline, pre = True )
            self.contents = Struct()
            for icol, colName in enumerate(self.df.columns):
                setattr ( self.contents, colName, self.df[colName] )

        ############################################################################################################################
        #### statistics
        ############################################################################################################################

        if not self._isProgressFile:

            self.stats = Struct()

            #### add chain cormat

            self._progress.note( msg = "computing the sample correlation matrix... ", end = newline, pre = True )
            self.stats.cormat = ccm.CorMat  ( dataFrame     = self.df
                                            , columns       = range(self._offset,self._offset+self.ndim)
                                            , methodName    = self._methodName
                                            , reportEnabled = self._reportEnabled
                                            , method        = "pearson"
                                            )
            self.stats.cormat()

            #### add chain covmat

            self._progress.note( msg = "computing the sample covariance matrix... ", end = newline, pre = True )
            self.stats.covmat = ccm.CovMat  ( dataFrame     = self.df
                                            , columns       = range(self._offset,self._offset+self.ndim)
                                            , methodName    = self._methodName
                                            , reportEnabled = self._reportEnabled
                                            )
            self.stats.covmat()

            #### add chain autocorrelation

            self.stats.maxLogFunc = getMaxLogFunc(dataFrame = self.df)

            self._progress.note( msg = "computing the sample autocorrelations... ", end = newline, pre = True )
            self.stats.autocorr = AutoCorr  ( dataFrame     = self.df
                                            , columns       = range(self._offset-1,self._offset+self.ndim)
                                            , methodName    = self._methodName
                                            , reportEnabled = self._reportEnabled
                                            )
            self.stats.autocorr()

        ############################################################################################################################
        #### graphics
        ############################################################################################################################

        self._plotTypeList =    [ "line"
                                , "scatter"
                                , "lineScatter"
                                ]
        
        if not self._isProgressFile: self._plotTypeList +=  [ "line3"
                                                            , "scatter3"
                                                            , "lineScatter3"
                                                            , "jointplot"
                                                            , "distplot"
                                                            , "kdeplot1"
                                                            , "kdeplot2"
                                                            , "contour3"
                                                            , "contourf"
                                                            , "contour"
                                                            , "grid"
                                                            ]

        self._progress.note( msg = "adding the graphics tools... ", end = newline, pre = True )
        self.plot = Struct()
        self._resetPlot(resetType="hard")

        self.plot.reset = self._resetPlot
        #self.plot.helpme = self.helpme

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

                    reset("hard")                   # regenerate all plots from scratch
                    reset("hard","line3")           # regenerate line3 plot from scratch
                    reset("hard",["line","line3"])  # regenerate line and line3 plots

        """

        requestedPlotTypeList = []
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

        if isinstance(resetType, str):
            resetTypeIsHard = resetType.lower()=="hard"
        else:
            resetTypeIsHard = None
            pm.abort( msg   = "The input argument resetType must be a string representing" + newline
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

            is3d        = "3"           in requestedPlotTypeLower
            isLine      = "line"        in requestedPlotTypeLower
            isScatter   = "scatter"     in requestedPlotTypeLower

            isJointplot = "jointplot"   in requestedPlotTypeLower
            isDistplot  = "distplot"    in requestedPlotTypeLower
            isKdeplot1  = "kdeplot1"    in requestedPlotTypeLower
            isKdeplot2  = "kdeplot2"    in requestedPlotTypeLower
            isContourf  = "contourf"    in requestedPlotTypeLower
            isContour3  = "contour3"    in requestedPlotTypeLower
            isContour   = "contour"     in requestedPlotTypeLower and not (isContourf or isContour3)

            isGridPlot  = "grid"        in requestedPlotTypeLower

            isLineScatterPlot = isLine or isScatter
            isDensityPlot = isJointplot or isDistplot or isKdeplot1 or isKdeplot2 or isContourf or isContour3 or isContour

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

                plotObject.ycolumns = self.df.columns[self._offset] # :]
                plotObject.ccolumns = self._sampleLogFuncColName
                plotObject.colorbar.kws.extend = "neither"
                plotObject.colorbar.kws.orientation = "vertical"
                plotObject.colorbar.kws.spacing = "uniform"

                if is3d:
                    plotObject.zcolumns = self._sampleLogFuncColName
                    if self.ndim>1:
                        plotObject.xcolumns = self.df.columns[self._offset]
                        plotObject.ycolumns = self.df.columns[self._offset+1]


                if isLine:
                    if isScatter:
                        plotObject.lineCollection.enabled = False
                        plotObject.plot.enabled = True
                        plotObject.plot.kws.alpha = 0.2
                        plotObject.plot.kws.color = "grey"
                        plotObject.plot.kws.linewidth = 0.75
                    else:
                        plotObject.lineCollection.enabled = True
                        plotObject.plot.enabled = False

            ########################################################################################################################
            #### reset density plots: kdeplot / distplot / jointplot / contour / contourf / contour3
            ########################################################################################################################

            if isDensityPlot:

                if resetTypeIsHard:
                    plotObject = DensityPlot( plotType = requestedPlotType
                                            , dataFrame = self.df
                                            , methodName = self._methodName
                                            , reportEnabled = self._reportEnabled
                                            , resetPlot = self._resetPlot
                                            )

                plotObject.xcolumns = self.df.columns[self._offset]
                if not (isDistplot or isKdeplot1):
                    if self.ndim==1:
                        plotObject.xcolumns = self.df.columns[self._offset-1]
                        plotObject.ycolumns = self.df.columns[self._offset]
                    else:
                        plotObject.ycolumns = self.df.columns[self._offset+1]

            ########################################################################################################################
            #### reset GridPlot
            ########################################################################################################################

            if isGridPlot:

                if resetTypeIsHard:
                    plotObject = GridPlot   ( plotType = requestedPlotType
                                            , dataFrame = self.df
                                            , methodName = self._methodName
                                            , reportEnabled = self._reportEnabled
                                            , resetPlot = self._resetPlot
                                            )

                endColindex = np.min( [self._offset+3, self._offset+self.ndim] )
                plotObject.columns = self.df.columns[self._offset-1:endColindex]
                plotObject.ccolumn = self._sampleLogFuncColName

            ########################################################################################################################
            ## reset target component
            ########################################################################################################################

            if (isLineScatterPlot or isDensityPlot) and not (plotObject._type.is3d or self._isProgressFile):

                xtarget = 0 # dummy
                if isDensityPlot: xtarget = self.df[plotObject.xcolumns].values.flatten()[self.stats.maxLogFunc.idrow]

                if plotObject._type.is1d: plotObject.target.value = [ xtarget, 0 ]

                if plotObject._type.is2d:
                    ytarget = self.df[plotObject.ycolumns].values.flatten()[self.stats.maxLogFunc.idrow]
                    plotObject.target.value = [ xtarget, ytarget ]

                if isDensityPlot and plotObject._type.is1d: plotObject.target.axhline.enabled = False
                if isLine or isScatter: plotObject.target.axvline.enabled = False
                plotObject.target.label = "maxLogFunc"

            ########################################################################################################################

            if plotObject is not None: setattr(self.plot, requestedPlotType, plotObject)

    ################################################################################################################################
    #### _reportWrongPlotName
    ################################################################################################################################

    def _reportWrongPlotName( self
                            , plotNames
                            ):

        pm.abort( msg   = "The input argument plotNames must be a string representing" + newline
                        + "the name of a plot belonging to the TabularFileContents class or," + newline
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

