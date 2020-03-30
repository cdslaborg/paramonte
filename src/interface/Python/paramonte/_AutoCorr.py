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
import weakref as _wref
#from scipy.signal import correlate as _ccor

import _message as msg
import _dfutils as dfutils
import _pmutils as pmutils
#from _LinePlot import LinePlot
#from _ScatterPlot import ScatterPlot

class _struct:
    pass

####################################################################################################################################
#### AutoCorr class
####################################################################################################################################

class AutoCorr:
    """
    This is the class for generating object of type CorMat which, 
    upon construction, will provide methods to compute and plot the 
    autocorrelations of the selected columns of the input dataFrame.

    Parameters
    ----------
        dataFrame
            a Pandas dataframe based upon the selected comlumns of which 
            the autocorrelations will be computed.
        columns
            optional argument that determines the columns of the input dataFrame to be 
            used in the computation of the autocorrelations. It can have three forms:
                1. a list of column indices from the input dataFrame.
                2. a list of column names from dataFrame.columns.
                3. a range(start,stop,step), representing the column indices in dataFrame.
            Examples:
                1. columns = [0,1,4,3]
                2. columns = ["SampleLogFunc","SampleVariable1"]
                3. columns = range(17,7,-2)
            If not provided, the default behavior includes all columns of the dataFrame.
        rows
            optional argument that determines the rows of the input dataFrame to be 
            used in the computation of the autocorrelations. It can be either:
                1. a range(start,stop,step), or, 
                2. a list of row indices from dataFrame.index.
            Examples:
                1. rows = range(17,7,-2)
                2. rows = [i for i in range(7,17)]
            If not provided, the default behavior includes all rows of the dataFrame.

    Attributes
    ----------
        All of the parameters described above, except dataFrame.
            a reference to the dataFrame will be implicitly stored in the object.
        df
            a pandas dataframe containing the computed autocorrelations.
        plot
            a structure containing the following plotting tools:
                heatmap
                    a callable object of class HeatMap which will enable 
                    plotting of the correlation matrix.

    Returns
    -------
        corMat
            an object of type class CorMat.

    ----------------------------------------------------------------------
    """

    def __init__( self
                , dataFrame     : _tp.Optional[ _pd.DataFrame ] = None
                , columns       : _tp.Optional[ _tp.Union[ range , _tp.List[int] , _tp.List[str] ] ] = None
                , rows          : _tp.Optional[ _tp.Union[ range , _tp.List[int] ] ] = None
                ):

        self._dfref = None if dataFrame is None else _wref.ref(dataFrame)
        self.columns = columns
        self.rows = rows
        self.df = None

    ################################################################################################################################

    def __call__( self
                , reself    : _tp.Optional[ bool ] = False
                , **kwargs
                ):
        """
        calls the get() method of the current instance of the class.

        Parameters
        ----------
            reself
                logical variable. If True, an instance of the object 
                will be returned upon exit to the calling routine.
                The default value is False.
            also, any attributes of the current instance of the class.

        Returns
        -------
            the object self if reself = True otherwise, None.
            However, this method causes side-effects by manipulating 
            the existing attributes of the object.

        """

        for key in kwargs.keys():
            if hasattr(self,key):
                setattr(self, key, kwargs[key])
            elif key=="dataFrame":
                setattr( self, "_dfref", _wref.ref(kwargs[key]) )
            else:
                raise Exception ( "Unrecognized input '"+key+"' class attribute detected.\n"
                                + "For allowed attributes, use help(objectname) on Python command line,\n"
                                + "where objectname should be replaced with the name of the object being called."
                                )
        self.get()
        if reself: return self

    ################################################################################################################################

    def get(self):
        """
        compute the autocorrelations of the selected columns 
        of the input dataframe to the object's constructor.

        Parameters
        ----------
            None

        Returns
        -------
            None. However, this method causes side-effects by manipulating 
            the existing attributes of the object.

        """

        # check columns presence

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
            raise Exception ( "The input argument 'columns' must be a list whose elements are all\n"
                            + "    1.   string-valued, each representing the name of the column from\n"
                            + "         the input dataframe to the object's constructor, to be included\n"
                            + "         in the computation of autocorrelations, or,\n"
                            + "    2.   integer-valued, each representing the index of the column from\n"
                            + "         the input dataframe to the object's constructor, to be included\n"
                            + "         in the computation of autocorrelations."
                            )

        # check rows presence

        if self.rows is None:
            rowindex = range(len(self._dfref().index))
        else:
            rowindex = self.rows

        # compute the autocorrelations

        nvar = len(colnames)
        nlag = len(rowindex)
        acf = _np.zeros((nvar,nlag))

        from scipy.signal import correlate as _ccor
        for icnt, ivar in enumerate(colindex):
            xdata = self._dfref().iloc[rowindex,ivar].values.flatten() - _np.mean(self._dfref().iloc[rowindex,ivar].values.flatten())
            acf[icnt] = _ccor   ( xdata
                                , xdata
                                , mode = "full"
                                )[nlag-1:2*nlag]
            acf[icnt] = acf[icnt] / acf[icnt,0]
        self.df = _pd.DataFrame(_np.transpose(acf))

        # specify columns/index names

        colnames = [ "ACF_"+colnames[i] for i in range(len(colnames)) ]
        self.df.columns = colnames

        # add SampleLogFunc to df for plot coloring, if exists

        ccolumns = ()
        if "SampleLogFunc" in self._dfref().columns:
            ccolumns = "SampleLogFunc"
            self.df.insert  ( loc = 0
                            , column = ccolumns
                            , value = self._dfref()[[ccolumns]].values.flatten()
                            , allow_duplicates = True
                            )

        # add lags to df

        self.df.insert  ( loc = 0
                        , column = "Lag"
                        , value = [ i for i in self.df.index ]
                        )

#!DEC$ ifdef PMVIS_ENABLED

        #####################
        #### add plots
        #####################

        self.plot = _struct()

        # add LinePlot

        from _LinePlot import LinePlot

        self.plot.line = LinePlot   ( dataFrame = self.df
                                    , xcolumns = "Lag"
                                    , ycolumns = colnames
                                    , ccolumns = None #ccolumns
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

        from _ScatterPlot import ScatterPlot

        self.plot.scatter = ScatterPlot ( dataFrame = self.df
                                        , xcolumns = "Lag"
                                        , ycolumns = colnames
                                        , ccolumns = None #ccolumns
                                        #, scatter_kws = {}
                                        , colorbar_kws =    {
                                                            "extend":"neither",
                                                            "orientation":"vertical",
                                                            #"spacing":"uniform",
                                                            }
                                        #, legend_kws = None
                                        )

    ################################################################################################################################

#!DEC$ endif