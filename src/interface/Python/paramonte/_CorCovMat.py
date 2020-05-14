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
####   we ask you to acknowledge the ParaMonte library's usage
####   in your work (education/research/industry/development/...)
####   by citing the ParaMonte library as described on this page:
####
####       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
####
####################################################################################################################################
####################################################################################################################################

import numpy as _np
import typing as _tp
import pandas as _pd
import weakref as _wref

import _message as msg
import _dfutils as dfutils
#from _HeatMapPlot import HeatMapPlot

class _struct:
    pass

####################################################################################################################################
#### CorCovMat class
####################################################################################################################################

class CorCovMat:

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
        compute the correlation matrix of the selected columns 
        of the input dataframe to the object's constructor.

        Parameters
        ----------
            None

        Returns
        -------
            None. However, this method causes side-effects by manipulating 
            the existing attributes of the object.

        """

        if hasattr(self,"method"):
            self._isCorMat = True
            self._matrixType = "correlation"
            if self.method not in ["pearson","kendall","spearman"]:
                raise Exception ( "The requested correlation type must be one of the following string values,\n"
                                + "    pearson  : standard correlation coefficient\n"
                                + "    kendall  : Kendall Tau correlation coefficient\n"
                                + "    spearman : Spearman rank correlation."
                                )
        else:
            self._isCorMat = False
            self._matrixType = "covariance"

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
                            + "         in the " + self._matrixType + " matrix construction, or,\n"
                            + "    2.   integer-valued, each representing the index of the column from\n"
                            + "         the input dataframe to the object's constructor, to be included\n"
                            + "         in the " + self._matrixType + " matrix construction."
                            )

        # check rows presence

        if self.rows is None:
            rowindex = range(len(self._dfref().index))
        else:
            rowindex = self.rows

        # construct the matrix dataframe

        if  self._isCorMat:
            self.df = self._dfref().iloc[rowindex,colindex].corr(method=self.method)
        else:
            self.df = self._dfref().iloc[rowindex,colindex].cov()

        # specify columns/index names

        self.df.columns = colnames
        self.df.index   = colnames

#!DEC$ ifdef PMVIS_ENABLED

        # add heatmap plot

        heatmap_kws = {}
        if self._isCorMat:

            annotPrecision = 2
            heatmap_kws["cbar_kws"] =   { "label": self.method.capitalize() + "'s Correlation Strength"
                                        , "orientation": "vertical"
                                        , "ticks": _np.linspace(-1,1,9)
                                        }

            heatmap_kws["vmin"]     = -1
            heatmap_kws["vmax"]     = +1
            heatmap_kws["center"]   = 0

        else:

            annotPrecision = None
            heatmap_kws["cbar_kws"] =   { "label": "Covariance Strength"
                                        , "orientation": "vertical"
                                        }

        from _HeatMapPlot import HeatMapPlot

        self.plot = _struct()
        self.plot.heatmap = HeatMapPlot ( dataFrame = self.df
                                        , heatmap_kws = heatmap_kws
                                        , annotPrecision = annotPrecision
                                        )

#!DEC$ endif

####################################################################################################################################
#### CorMat class
####################################################################################################################################

class CorMat(CorCovMat):
    """
    This is the class for generating object of type CorMat which, 
    upon construction, will provide methods to compute and plot the 
    correlation matrix of the selected columns of the input dataFrame.

    Parameters
    ----------
        dataFrame
            a Pandas dataframe based upon the selected comlumns of which 
            the correlation matrix will be computed.
        columns
            optional argument that determines the columns of the input dataFrame to be 
            used in the computation of the correlation matrix. It can have three forms:
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
            used in the computation of the correlation matrix. It can be either:
                1. a range(start,stop,step), or, 
                2. a list of row indices from dataFrame.index.
            Examples:
                1. rows = range(17,7,-2)
                2. rows = [i for i in range(7,17)]
            If not provided, the default behavior includes all rows of the dataFrame.
        method
            optional string value representing the method to be used 
            for the computation of correlations:
                1. 'pearson'  : standard correlation coefficient, 
                2. 'kendall'  : Kendall Tau correlation coefficient, 
                3. 'spearman' : Spearman rank correlation.
        The default value is 'pearson'.

    Attributes
    ----------
        All of the parameters described above, except dataFrame.
            a reference to the dataFrame will be implicitly stored in the object.
        df
            a pandas dataframe containing the computed correlation matrix
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
                , dataFrame     : _pd.DataFrame
                , columns       : _tp.Optional[ _tp.Union[ range , _tp.List[int] , _tp.List[str] ] ] = None
                , rows          : _tp.Optional[ _tp.Union[ range , _tp.List[int] ] ] = None
                , method        : _tp.Optional[str] = "pearson"
                ):

        self.method = method
        super().__init__( dataFrame     = dataFrame
                        , columns       = columns
                        , rows          = rows
                        )

####################################################################################################################################
#### CovMat class
####################################################################################################################################

class CovMat(CorCovMat):
    """
    This is the class for generating object of type CovMat which, 
    upon construction, will provide methods to compute and plot the 
    covariance matrix of the selected columns of the input dataFrame.

    Parameters
    ----------
        dataFrame
            a Pandas dataframe based upon the selected comlumns of which 
            the covariance matrix will be computed.
        columns
            optional argument that determines the columns of the input dataFrame to be 
            used in the computation of the covariance matrix. It can have three forms:
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
            used in the computation of the covariance matrix. It can be either:
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
            a pandas dataframe containing the computed covariance matrix
        plot
            a structure containing the following plotting tools:
                heatmap
                    a callable object of class HeatMap which will enable 
                    plotting of the correlation matrix.

    Returns
    -------
        covMat
            an object of type class CovMat.

    ----------------------------------------------------------------------
    """

    def __init__( self
                , dataFrame     : _pd.DataFrame
                , columns       : _tp.Optional[ _tp.Union[ range , _tp.List[int] , _tp.List[str] ] ] = None
                , rows          : _tp.Optional[ _tp.Union[ range , _tp.List[int] ] ] = None
                ):

        super().__init__( dataFrame     = dataFrame
                        , columns       = columns
                        , rows          = rows
                        )
