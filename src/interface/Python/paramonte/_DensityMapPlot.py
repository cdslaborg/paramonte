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

import _message as msg
import _dfutils as dfutils

#!DEC$ ifdef PMVIS_ENABLED

#import seaborn as _sns
#from matplotlib import pyplot as _plt

from _visutils import Target, ParaMonteFigure


####################################################################################################################################
#### DensityMapPlot class
####################################################################################################################################

class DensityMapPlot:
    """
    This is the DensityMapPlot class for generating instances 
    of 2D densitymap figures based on seaborns library's 
    kdeplot() and functions.

    Usage
    -----
    first generate an object of this class by optionally 
    passing the following parameters described below. Then call 
    the plot() method. The generated object is also callable with 
    the same input parameters as the object's constructor.

    Parameters
    ----------
        dataFrame
            a Pandas dataframe from which the selected data will be plotted.
        xcolumn
            optional argument that determines the columns of dataframe to serve as 
            the x-values. It can have three forms:
                1. a column index of the input dataframe.
                2. a column name of the input dataframe (as in dataFrame.columns).
            Examples:
                1. xcolumn = 0
                2. xcolumn = [0]
                3. xcolumn = ["SampleVariable1"]
        If not provided, the default will be the first column of the input dataFrame.
        ycolumn
            optional argument that determines the columns of dataframe to serve as 
            the y-values. It can have three forms:
                1. a column index of the input dataframe.
                2. a column name of the input dataframe (as in dataFrame.columns).
            Examples:
                1. ycolumn = 0
                2. ycolumn = [0]
                3. ycolumn = ["SampleVariable2"]
        If not provided, the default will be the second column of the input dataFrame.
        rows
            optional argument that determines the rows of dataframe 
            to be visualized. It can be either:
                1. a range(start,stop,step), or, 
                2. a list of row indices in dataFrame.index.
            Examples:
                1. rows = range(17,7,-2)
                2. rows = [i for i in range(7,17)]
            If not provided, the default will include all rows of the dataframe.
        set_kws
            optional dictionary of keyword arguments to be passed to seaborn's 
            set() function. For example: 
                set_kws = {"style":"darkgrid"}
            if set to None, then no call to set() function will be made.
        figure_kws
            optional dictionary of keyword arguments to be passed to matplotlib's 
            figure() function. For example: 
                figure_kws = {"facecolor":"w","dpi":150}
        kdeplot_kws
            optional dictionary of keyword arguments to be passed to seaborn's 
            kdeplot() function. For example: 
                kdeplot_kws = {"n_levels":100,"shade":True,"shade_lowest":True,"cmap":"Blues"}
            If not set, the default will be None and no kdeplot will be plotted.
            If both kdeplot_kws and jointplot_kws are set, a jointplot will be plotted.
        jointplot_kws
            optional dictionary of keyword arguments to be passed to seaborn's 
            jointplot() function. For example: 
                kdeplot_kws = {"n_levels":100,"shade":True,"shade_lowest":True,"cmap":"Blues"}
            If not set, the default will be None and no jointplot will be plotted.
            If both kdeplot_kws and jointplot_kws are set, a jointplot will be plotted.
        outputFile
            optional string representing the name of the output file in which 
            the figure will be saved. If not provided, no output file will be generated.

    Attributes
    ----------
        All of the parameters described above, except dataFrame.
            a reference to the dataFrame will be implicitly stored in the object.
        target
            a callable object of ParaMonte library's class Target which can 
            be used to add target point and lines to the current active plot.
        currentFig
            an object of class ParaMonteFigure which is initially None, but upon 
            making a plot, is populated with attributes representing all outputs 
            from matplotlib/seaborn function calls, with the following attributes:
                figure
                    the output of matplotlib's figure() function.
                kdeplot
                    the output of seaborn's kdeplot() function.
                jointplot
                    the output of seaborn's jointplot() function.

    Returns
    -------
        DensityMapPlot

    ----------------------------------------------------------------------
    """

    def __init__( self
                , dataFrame     : _tp.Optional[ _pd.DataFrame ] = None
                , xcolumn       : _tp.Optional[ _tp.Union[ str, range , _tp.List[int] , _tp.List[str] ] ] = None
                , ycolumn       : _tp.Optional[ _tp.Union[ str, range , _tp.List[int] , _tp.List[str] ] ] = None
                , rows          : _tp.Optional[ _tp.Union[ range , _tp.List[int] ] ] = None
                , set_kws       : _tp.Optional[_tp.Dict] = ()
                , figure_kws    : _tp.Optional[_tp.Dict] = ()
                , kdeplot_kws   : _tp.Optional[_tp.Dict] = ()
                , jointplot_kws : _tp.Optional[_tp.Dict] = ()
                , outputFile    : _tp.Optional[str] = None
                ):

        self._dfref = None if dataFrame is None else _wref.ref(dataFrame)
        self.xcolumn        = xcolumn
        self.ycolumn        = ycolumn
        self.rows           = rows
        self.set_kws        = set_kws
        self.figure_kws     = figure_kws
        self.kdeplot_kws    = kdeplot_kws
        self.jointplot_kws  = jointplot_kws
        self.outputFile     = outputFile
        self.currentFig     = None
        self.target         = Target()

        self._isdryrun = True
        self.plot()
        self._isdryrun = False

    ################################################################################################################################

    def __call__( self
                , reself    : _tp.Optional[ bool ] = False
                , **kwargs
                ):
        """
        calls the plot() method of the current instance of the class.

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
        self.plot()
        if reself: return self

    ################################################################################################################################

    def plot(self):
        """
        Generate a line plot, or scatter plot, or both 
        from the selected columns of the object's dataframe.

        Parameters
        ----------
            None

        Returns
        -------
            None. However, this method causes side-effects by manipulating 
            the existing attributes of the object.

        """

        # set plot type

        kdeEnabled = self.kdeplot_kws is not None

        if self.jointplot_kws is None:
            jointEnabled = False # not kdeEnabled
        else:
            kdeEnabled = False
            jointEnabled = True

        # set what to plot

        if self.kdeplot_kws==(): self.kdeplot_kws={}
        if isinstance(self.kdeplot_kws,dict):
            if "cbar"           not in self.kdeplot_kws.keys(): self.kdeplot_kws["cbar"]          = False
            if "cmap"           not in self.kdeplot_kws.keys(): self.kdeplot_kws["cmap"]          = "Blues"
            if "shade"          not in self.kdeplot_kws.keys(): self.kdeplot_kws["shade"]         = True
            if "n_levels"       not in self.kdeplot_kws.keys(): self.kdeplot_kws["n_levels"]      = 100
            if "shade_lowest"   not in self.kdeplot_kws.keys(): self.kdeplot_kws["shade_lowest"]  = True
        elif kdeEnabled:
            raise Exception ( "The input argument 'kdeplot_kws' must be None or a dictionary,\n"
                            + "each (key,value) pair of which represents an (attribute,value)\n"
                            + "of matplotlib library's plot() function."
                            )

        if self.jointplot_kws==(): self.jointplot_kws={}
        if isinstance(self.jointplot_kws,dict):
            if "height"             not in self.jointplot_kws.keys(): self.jointplot_kws["height"]        = 7
            if "space"              not in self.jointplot_kws.keys(): self.jointplot_kws["space"]         = 0
            if "kind"               not in self.jointplot_kws.keys(): self.jointplot_kws["kind"]          = "kde"
            if self.jointplot_kws["kind"] == "kde":
                if "cbar"           not in self.jointplot_kws.keys(): self.jointplot_kws["cbar"]          = False
                if "cmap"           not in self.jointplot_kws.keys(): self.jointplot_kws["cmap"]          = "Blues"
                if "shade"          not in self.jointplot_kws.keys(): self.jointplot_kws["shade"]         = True
                if "n_levels"       not in self.jointplot_kws.keys(): self.jointplot_kws["n_levels"]      = 100
                if "shade_lowest"   not in self.jointplot_kws.keys(): self.jointplot_kws["shade_lowest"]  = True
        elif jointEnabled:
            raise Exception ( "The input argument 'jointplot_kws' must be None or a dictionary,\n"
                            + "each (key,value) pair of which represents an (attribute,value)\n"
                            + "of matplotlib library's LineCollection() class."
                            )

        #########################

        if self.set_kws==(): self.set_kws={}
        if "style" not in self.set_kws.keys(): self.set_kws["style"] = "ticks"
        figEnabled = self.figure_kws is not None
        if figEnabled:
            if self.figure_kws==(): self.figure_kws={}
            if isinstance(self.figure_kws,dict):
                if "dpi" not in self.figure_kws.keys(): self.figure_kws["dpi"] = 150
                if "facecolor" not in self.figure_kws.keys(): self.figure_kws["facecolor"] = "w"
                if "edgecolor" not in self.figure_kws.keys(): self.figure_kws["edgecolor"] = "w"
            else:
                raise Exception ( "The input argument 'figure_kws' must be a dictionary,\n"
                                + "each (key,value) pair of which represents an (attribute,value)\n"
                                + "of matplotlib library's figure() function."
                                )

        #########################
        if self._isdryrun: return
        #########################

        import seaborn as _sns
        from matplotlib import pyplot as _plt

        # generate figure and axes if needed

        self.currentFig = ParaMonteFigure()
        if figEnabled:
            if self.set_kws is not None: _sns.set(**self.set_kws)
            if kdeEnabled: self.currentFig.figure = _plt.figure( **self.figure_kws )
            #self.currentFig.axes = _plt.subplot(1,1,1)
        #else:
        #    self.currentFig.axes = _plt.gca()
        #    self.currentFig.figure = self.currentFig.axes.get_figure()

        # check data type

        noDataPassed = False
        if self._dfref is None:
            noDataPassed = True
            fatalmsg = "It appears that no data has been passed for plotting.\n"
        elif not isinstance(self._dfref,_wref.ref):
            noDataPassed = True
            fatalmsg = "It appears that you have messed with the\ninternal representation of data in the object.\n"
        elif not isinstance(self._dfref(),_pd.DataFrame):
            noDataPassed = True
            fatalmsg = ""
        if noDataPassed:
            raise Exception ( fatalmsg
                            + "The input data must be a pandas' dataframe.\n"
                            + "Please pass a dataFrame to the constructor or at\n"
                            + "the time of calling the object (which is callable with\n"
                            + "the same input arguments as the object's constructor."
                            )

        #########################

        # check rows presence

        if self.rows is None:
            rowindex = range(len(self._dfref().index))
        else:
            rowindex = self.rows

        # check columns presence

        if self.xcolumn is None:
            xcolindex = [0]
            xcolnames = list(self._dfref().columns[xcolindex[0]:xcolindex[0]+1])
        else:
            xcolnames, xcolindex = dfutils.getColNamesIndex(self._dfref().columns,self.xcolumn)

        if self.ycolumn is None:
            ycolindex = [1]
            ycolnames = list(self._dfref().columns[ycolindex[0]:ycolindex[0]+1])
        else:
            ycolnames, ycolindex = dfutils.getColNamesIndex(self._dfref().columns,self.ycolumn)

        # check only one column is present

        if len(xcolnames)>1:
            raise Exception ( "Only one data column can be plotted on the x-axis.\n"
                            + "You have provided " + str(len(xcolnames)) + " columns:\n"
                            + str(xcolnames)
                            )
        if len(ycolnames)>1:
            raise Exception ( "Only one data column can be plotted on the y-axis.\n"
                            + "You have provided " + str(len(ycolnames)) + " columns:\n"
                            + str(ycolnames)
                            )

        # make plots

        if kdeEnabled:

            self.currentFig.kdeplot = _sns.kdeplot  ( data = self._dfref().iloc[rowindex,xcolindex[0]]
                                                    , data2 = self._dfref().iloc[rowindex,ycolindex[0]]
                                                    , **self.kdeplot_kws
                                                    )

        if jointEnabled:

            self.currentFig.jointplot = _sns.jointplot  ( x = self._dfref().iloc[rowindex,xcolindex[0]]
                                                        , y = self._dfref().iloc[rowindex,ycolindex[0]]
                                                        , **self.jointplot_kws
                                                        )

        _plt.tight_layout()
        if self.outputFile is not None:
            self.currentFig.figure.savefig  ( self.outputFile
                                            , bbox_inches = 'tight'
                                            , pad_inches = 0.0
                                            )

    ################################################################################################################################

#!DEC$ endif