#**********************************************************************************************************************************
#**********************************************************************************************************************************
#
#  ParaMonte: plain powerful parallel Monte Carlo library.
#
#  Copyright (C) 2012-present, The Computational Data Science Lab
#
#  This file is part of ParaMonte library.
#
#  ParaMonte is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, version 3 of the License.
#
#  ParaMonte is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
#
#**********************************************************************************************************************************
#**********************************************************************************************************************************

import numpy as _np
import typing as _tp
import pandas as _pd

import _paramonte as _pm

_timer = _pm.pmutils.Timer(_methodName=_pm.names.paradram)

####################################################################################################################################
#### _ParaDRAMChain class
####################################################################################################################################

class _Struct:
    pass

class _ParaDRAMChain:
    """
    This is the _ParaDRAMChain class for generating instances 
    of ParaDRAM sample/chain. The ParaDRAM class's readSample() or 
    readChain() methods return an object or a list of objects of 
    class _ParaDRAMChain.

    Parameters
    ----------
        file
            full path to the file containing the sample/chain.
        delimiter
            the delimiter used in the sample/chain file, which 
            must be provided by the user.
        parseContents
            If set to True, the contents of the file will be parsed and 
            stored in a component of the object named 'contents'.
            The default value is True.
        markovChainRequested
            boolean value indicating weather the full Markov Chain
            has to be generated from the potentially-weighted sample 
            in the input file. If True, each sampled state in the resulting 
            output dataframe is guaranteed to have a weight of 1.

    Attributes
    ----------
        file
            full path to the file containing the sample/chain.
        delimiter
            the delimiter used in the sample/chain file, which 
            must be provided by the user.
        ndim
            number of dimensions of the domain of the objective 
            function from which the sample has been drawn.
        count
            number of points (states) in the sample/chain file. 
            This is essentially, the number of rows in the file 
            minus one (representing the header line).
        [df]
            if the input file contents is structured in a format that
            could be read as a dataframe, then the contents of the file
            will be stored in the form of a pandas-library DataFrame 
            in this property (hence called 'df').
        [contents]
            if the input file contents is structured in the form of columns,
            then a property named 'contents' is also added to the object.
            Each component of contents will named via the header of the file
            and will contain data from the corresponding column of the file.

    Returns
    -------
        _ParaDRAMChain

    ----------------------------------------------------------------------
    """

    def __init__( self
                , file
                , delimiter
                , parseContents = True
                , markovChainRequested = False
                , mpiDisabled = True
                ):

        ############################################################################################################################
        #### data
        ############################################################################################################################

        self.file = file
        self.delimiter = delimiter

        #with open(file, 'r') as targetFile: Line = targetFile.readlines()
        #self.colHeaderList = Line[0][:-1].split(delimiter)
        #ndimPlusOffset = len(self.colHeaderList)
        #self.dataMat = _np.zeros( (self.count,ndimPlusOffset) )
        #for isample,line in enumerate(Line[1:]): # ignore the header line
        #    self.dataMat[isample][0:ndimPlusOffset] = _np.array( line[:-1].split(delimiter) )
        #self.dataMat = self.dataMat.T

        _timer.tic( msg = "reading file contents... " )

        self.df =_pd.read_csv   ( self.file
                                , delimiter = self.delimiter
                                , header = 0
                                )
        self._offset = list(self.df.columns).index("SampleLogFunc") + 1 # index of the first variable
        self.ndim   = len(self.df.columns) - self._offset
        self.count  = len(self.df.iloc[:,1])

        if markovChainRequested:

            CumSumWeight = _np.cumsum(self.df.iloc[:,self._offset-2].values, dtype=_np.int32)
            if CumSumWeight[-1] != self.count: # it is indeed a compact chain
                #dfMarkov = _pd.DataFrame( columns=list(self.df.columns), index=list(range(CumSumWeight[-1])) )
                dfMarkov = _np.zeros( (CumSumWeight[-1] , self.ndim+self._offset) )
                istart = 0
                for i in range(self.count):
                    iend = CumSumWeight[i]
                    #dfMarkov.iloc[istart:iend,:] = self.df.iloc[i].values
                    dfMarkov[istart:iend,:] = self.df.iloc[i].values
                    istart = iend
                columns = self.df.columns
                self.df = _pd.DataFrame(dfMarkov)
                self.count = len(self.df.iloc[:,1])
                self.df.columns = columns

        _timer.toc()

        if not mpiDisabled:
            _pm.note( msg = "ndim = " + str(self.ndim) + ", count = " + str(self.count)
                    , methodName = _pm.names.paradram
                    , marginTop = 0
                    , marginBot = 1
                    , end = ""
                    )

        # set dynamic properties

        if parseContents:
            _timer.tic( msg = "parsing file contents... " )
            self.contents = _Struct()
            for icol, colName in enumerate(self.df.columns):
                setattr ( self.contents, colName, self.df[colName] )
            _timer.toc()

        ############################################################################################################################
        #### statistics
        ############################################################################################################################

        self.stats = _Struct()
        
        # add chain cormat

        self.stats.cormat = _pm.stats.CorMat( dataFrame     = self.df
                                            , columns       = range(self._offset,self._offset+self.ndim)
                                            , method        = "pearson"
                                            )
        
        _timer.tic( msg = "computing sample correlation matrix... " )
        self.stats.cormat()
        _timer.toc()

        # add chain covmat

        self.stats.covmat = _pm.stats.CovMat( dataFrame     = self.df
                                            , columns       = range(self._offset,self._offset+self.ndim)
                                            )

        _timer.tic( msg = "computing sample covariance matrix... " )
        self.stats.covmat()
        _timer.toc()

        self.stats.maxLogFunc = _pm.stats.getMaxLogFunc(dataFrame = self.df)
        #self.stats.maxLogFunc = _Struct()
        #self.stats.maxLogFunc.idrow = self.df[["SampleLogFunc"]].idxmax().values[0]
        #self.stats.maxLogFunc.value = self.df[["SampleLogFunc"]].iat[self.stats.maxLogFunc.idrow,0]
        #self.stats.maxLogFunc.dfrow = self.df.iloc[self.stats.maxLogFunc.idrow,:]
        #self.stats.maxLogFunc.state = self.df.iloc[self.stats.maxLogFunc.idrow,self._offset:]

        # add chain autocorrelation

        self.stats.acf = _pm.stats.AutoCorr ( dataFrame     = self.df
                                            , columns       = range(self._offset-1,self._offset+self.ndim)
                                            )

        _timer.tic( msg = "computing autocorrelations... " )
        self.stats.acf()
        _timer.toc()

#!DEC$ ifdef PMVIS_ENABLED

        ############################################################################################################################
        #### graphics
        ############################################################################################################################

        _timer.tic( msg = "adding graphics tools... " )

        # add HistPlot

        self.plot = _Struct()

        self.plot.hist = _pm.vis.HistPlot   ( dataFrame = self.df
                                            , columns = self.df.columns[self._offset:]
                                            )

        # add LinePlot

        self.plot.line = _pm.vis.LinePlot   ( dataFrame = self.df
                                            , ycolumns = self.df.columns[self._offset:]
                                            , ccolumns = "SampleLogFunc"
                                            , lc_kws =  {
                                                        #"linewidth":0.75,
                                                        #"cmap":"viridis",
                                                        "cmap":"autumn",
                                                        #"alpha":0.5,
                                                        }
                                            , colorbar_kws =    {
                                                                "extend":"neither",
                                                                "orientation":"vertical",
                                                                #"spacing":"uniform",
                                                                }
                                            #, legend_kws = None
                                            )

        # add ScatterPlot

        self.plot.scatter = _pm.vis.ScatterPlot ( dataFrame = self.df
                                                , ycolumns = self.df.columns[self._offset:]
                                                , ccolumns = "SampleLogFunc"
                                                #, scatter_kws = {}
                                                , colorbar_kws =    {
                                                                    "extend":"neither",
                                                                    "orientation":"vertical",
                                                                    #"spacing":"uniform",
                                                                    }
                                                #, legend_kws = None
                                                )

        # add DensityMapPlot

        xindex = self._offset
        yindex = self._offset + 1
        if self.ndim==1: xindex, yindex = yindex-1, xindex-1
            
        self.plot.density = _pm.vis.DensityMapPlot  ( dataFrame = self.df
                                                    , xcolumn = xindex
                                                    , ycolumn = yindex
                                                    )
        #print(self.stats.maxLogFunc.idrow,xindex,yindex)
        #print(self.df)
        #print(self.df.iat[self.stats.maxLogFunc.idrow,xindex])
        #print(self.df.iat[self.stats.maxLogFunc.idrow,yindex])
        self.plot.density.target.__init__   ( value = [ self.df.iat[self.stats.maxLogFunc.idrow,xindex], self.df.iat[self.stats.maxLogFunc.idrow,yindex] ]
                                            , scatter_kws = {"label":"maxLogFunc"}
                                            )

        # add GridPlot

        endColindex = _np.min( [self._offset+3, self._offset+self.ndim] )
        self.plot.grid = _pm.vis.GridPlot   ( dataFrame = self.df
                                            , columns = self.df.columns[self._offset:endColindex]
                                            , scatterplot_kws = {"ccolumns": "SampleLogFunc"}
                                            , _methodName = _pm.names.paradram
                                            )

        # add ScatterLinePlot

        self.plot._scatterline = _pm.vis.ScatterLinePlot( dataFrame = self.df
                                                        , ycolumns = self.df.columns[self._offset:]
                                                        , lccolumns = "SampleLogFunc"
                                                        , lc_kws =  {
                                                                    #"linewidth":0.75,
                                                                    #"cmap":"viridis",
                                                                    "cmap":"autumn",
                                                                    #"alpha":0.5,
                                                                    }
                                                        #, scatter_kws = {}
                                                        , colorbar_kws =    {
                                                                            "extend":"neither",
                                                                            "orientation":"vertical",
                                                                            #"spacing":"uniform",
                                                                            }
                                                        #, legend_kws = None
                                                        )

        _timer.toc()

    ################################################################################################################################

#!DEC$ endif
