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

from _pmutils import Timer
import _message as msg
import _dfutils as dfutils

#!DEC$ ifdef PMVIS_ENABLED

#import seaborn as _sns
#from matplotlib import pyplot as _plt

from _HistPlot import HistPlot
from _LinePlot import LinePlot
from _ScatterPlot import ScatterPlot
from _visutils import Target, ParaMonteFigure

_eol = ""

def hide_current_axis(*args, **kwds):
    from matplotlib import pyplot as _plt
    _plt.gca().set_visible(False)

def show_current_axis(*args, **kwds):
    from matplotlib import pyplot as _plt
    _plt.gca().set_visible(True)

####################################################################################################################################
#### GridPlot class
####################################################################################################################################

class GridPlot:
    """
    This is the GridPlot class for generating instances 
    of matrix-shaped figures of the following form:
        -   the upper triangle subplots can optionally represent nothing 
            or the birvariate line-and/or-scatter plots of data.
        -   the diagonal elements represent the histograms of data.
        -   the lower triangle subplots can optionally represent nothing 
            or the birvariate line-and/or-scatter plots of data.

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
        columns
            optional argument that determines the columns of dataframe to serve as 
            the x-values. It can have three forms:
                1.  a list of column indices in the input dataFrame.
                2.  a list of column names in dataFrame.columns.
                3.  a range(start,stop,step), representing the column indices 
                    in the input dataFrame.
            Examples:
                1.  columns = [0,1,4,3]
                2.  columns = ["SampleLogFunc","SampleVariable1"]
                3.  columns = range(17,7,-2)
            However, in all cases, it must have a length that is either 1 or equal 
            to the length of ycolumns. If the length is 1, then xcolumns will be 
            plotted against data corresponding to each element of ycolumns.
            If not provided, the default will be the count of the rows of the 
            input dataFrame.
        rows
            optional argument that determines the rows of dataframe 
            to be visualized. It can be either:
                1.  a range(start,stop,step), or, 
                2.  a list of row indices in dataFrame.index.
            Examples:
                1. rows = range(17,7,-2)
                2. rows = [i for i in range(7,17)]
            If not provided, the default includes all rows of the dataframe.
        set_kws
            optional dictionary of keyword arguments to be passed to seaborn's 
            set() function. For example: 
                set_kws = {"style":"darkgrid"}
            if set to None, then no call to set() function will be made.
        pairgrid_kws
            optional dictionary of keyword arguments to be passed to the PairGrid 
            class of seaborn library for making a grid of axes on which the data
            will be plotted.
            The default is {}. It is recommended that you leave this parameter 
            unchanged from the default value, otherwise, unpredicted behavior
            may occur.
        distplot_kws
            optional dictionary of keyword arguments to be passed to distplot()
            function of seaborn library for making histograms of data on the 
            diagonal elements of the grid plot. For example:
                distplot_kws =  { "kde":False
                                , hist_kws = { "histtype": "stepfilled"
                                             , "linewidth": 1
                                             }
                                }
            The default is {}. 
            If set to None, no histograms will be added to diagonal of the grid plot. 
        kdeplot_kws
            optional dictionary of keyword arguments to be passed to kdeplot()
            function of seaborn library for making histograms of data on the 
            diagonal elements of the grid plot. For example:
                kdeplot_kws =  {"n_levels": 100,
                                "shade": True,
                                "shade_lowest": True,
                                "cmap": "Blues",
                                "cbar": False,
                                "zorder": 1,
                                alpha": 1,
                                }
            The default is {}. 
            If set to None, no kdeplots will be added to the grid plot. 
        lineplot_kws
            optional dictionary of keyword arguments to be passed to an instance of 
            ParaMonte's LinePlot class for making colored line plots of the requested 
            columns of data on the upper triangle of the grid plot.
            The default is {}.
            If set to None, no line plots will be added to upper triangle. 
        scatterplot_kws
            optional dictionary of keyword arguments to be passed to an instance of 
            ParaMonte's ScatterPlot class for making scatter plots of the requested 
            columns of data on the upper triangle of the grid plot.
            The default is {}.
            If set to None, no scatter plots will be added to upper triangle. 
        colorbar_kws
            optional dictionary of keyword arguments to be passed to matplotlib's 
            figure.colorbar() function. For example: 
                colorbar_kws = {"orientation":"vertical"}
            The default is {}. Despite the potential existence of three types of colored 
            plots on the grid plot, only one colorbar will be added to the figure.
            The colorbar preference is in the following order:
                1. colored scatterplot, if requested,
                2. colored lineplot, if requested,
                3. kdeplot, if requested,
            If set to None, no colorbar will be plotted.
        outputFile
            optional string representing the name of the output file in which 
            the figure will be saved. If not provided, no output file will be generated.
        kdecorner
            optional string representing the corner (upper/lower triangle) of the grid 
            plot where the kdeplots will be added. 
            Possible values are: "upper", "lower", "auto", None. 
            The default is "auto", which is the lower or any empty corner. 
            If set to "auto", the corner will be automatically determined. 
            If set to None, no kdeplots will be added to the grid plot. 
        lscorner
            optional string representing the corner (upper/lower triangle) of the grid 
            plot where the line/scatter plots will be added. 
            Possible values are: "upper", "lower", "auto", None. 
            The default is "auto", which is the upper or any empty corner. 
            If set to "auto", the corner will be automatically determined. 
            If set to None, no scatter or line plots will be added to the grid plot. 

    Attributes
    ----------
        All of the parameters described above, except dataFrame.
            a reference to the dataFrame will be implicitly stored in the object.
        addTarget
            a method which calls the individual callable Target objects of each subplot 
            to add target point and lines to all of the plots.
        currentFig
            an object of class ParaMonteFigure which is initially None, but upon 
            making a plot, is populated with attributes representing all outputs 
            from matplotlib/seaborn function calls, with the following attributes:
                figure
                    the output of matplotlib's figure() function.
                addTarget
                    a method which calls an object of ParaMonte library's Target 
                    class to add target point and lines to all plots of the grid.
                pairgrid
                    the output of seaborn's PairGrid() class.
                scatterplotList
                    a list each element of which corresponds to one row of the PairGrid, 
                    and is itself a list of scatterplot objects.
                lineplotList
                    a list each element of which corresponds to one row of the PairGrid, 
                    and is itself a list of lineplot objects.
                colorbar
                    the output of matplotlib's Figure.colorbar() function.

    Returns
    -------
        GridPlot

    ----------------------------------------------------------------------
    """

    def __init__( self
                , dataFrame         : _tp.Optional[ _pd.DataFrame ] = None
                , columns           : _tp.Optional[ _tp.Union[ str, range , _tp.List[int] , _tp.List[str] ] ] = None
                , rows              : _tp.Optional[ _tp.Union[ range , _tp.List[int] ] ] = None
                , set_kws           : _tp.Optional[_tp.Dict] = ()
                , pairgrid_kws      : _tp.Optional[_tp.Dict] = ()
                , kdeplot_kws       : _tp.Optional[_tp.Dict] = ()
                , distplot_kws      : _tp.Optional[_tp.Dict] = ()
                , lineplot_kws      : _tp.Optional[_tp.Dict] = ()
                , scatterplot_kws   : _tp.Optional[_tp.Dict] = ()
                , colorbar_kws      : _tp.Optional[_tp.Dict] = ()
                , outputFile        : _tp.Optional[str] = None
                , kdecorner         : _tp.Optional[str] = "auto"
                , lscorner          : _tp.Optional[str] = "auto"
                , _methodName       : _tp.Optional[str] = ""
                ):

        self._dfref = None if dataFrame is None else _wref.ref(dataFrame)
        self.columns            = columns
        self.rows               = rows
        self.set_kws            = set_kws
        self.pairgrid_kws       = pairgrid_kws
        self.kdeplot_kws        = kdeplot_kws
        self.distplot_kws       = distplot_kws
        self.lineplot_kws       = lineplot_kws
        self.scatterplot_kws    = scatterplot_kws
        self.colorbar_kws       = colorbar_kws
        self.outputFile         = outputFile
        self.kdecorner          = kdecorner
        self.lscorner           = lscorner
        self.currentFig         = None
        self._methodName        = _methodName

        self._timer = Timer(_methodName = _methodName)

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
        Generate a grid plot from the selected columns of the object's dataframe.

        Parameters
        ----------
            None

        Returns
        -------
            None. However, this method causes side-effects by manipulating 
            the existing attributes of the object.

        """

        # setup figure

        scatterplotEnabled = self.scatterplot_kws is not None
        lineplotEnabled = self.lineplot_kws is not None
        distplotEnabled = self.distplot_kws is not None
        kdeplotEnabled = (self.kdeplot_kws is not None) and (self.kdecorner is not None)
        lsplotEnabled = (scatterplotEnabled or lineplotEnabled) and (self.lscorner is not None)

        if self.set_kws==(): self.set_kws={}

        kdeUpperEnabled = False
        kdeLowerEnabled = False
        if kdeplotEnabled:
            if self.kdecorner=="auto":
                kdeLowerEnabled = True
            elif isinstance(self.kdecorner,str):
                if "upper" in self.kdecorner: kdeUpperEnabled = True
                if "lower" in self.kdecorner: kdeLowerEnabled = True
            else:
                raise Exception ( "The input argument 'kdecorner' must be either 'lower', 'upper' or\n"
                                + "a concatanation of both, or otherwise None." )

        lsUpperEnabled = False
        lsLowerEnabled = False
        if lsplotEnabled:
            if self.lscorner=="auto":
                lsUpperEnabled = not kdeUpperEnabled
                lsLowerEnabled = not kdeLowerEnabled
                lsplotEnabled = lsLowerEnabled or lsUpperEnabled
            elif isinstance(self.lscorner,str):
                if "upper" in self.lscorner: lsUpperEnabled = True
                if "lower" in self.lscorner: 
                    lsLowerEnabled = True
                    if self.kdecorner=="auto":
                        kdeUpperEnabled = not lsUpperEnabled
                        kdeLowerEnabled = not lsLowerEnabled
                        kdeplotEnabled = kdeLowerEnabled or kdeUpperEnabled
            else:
                raise Exception ( "The input argument 'lscorner' must be either 'lower', 'upper' or\n"
                                + "a concatanation of both, or otherwise None." )
        elif kdeplotEnabled:
            if self.kdecorner=="auto":
                kdeUpperEnabled = True
                kdeLowerEnabled = True
            


        if self.kdeplot_kws==(): self.kdeplot_kws={}
        if isinstance(self.kdeplot_kws,dict):
            self.kdeplot_kws["cbar"] = False
            if "cmap" not in self.kdeplot_kws.keys():           self.kdeplot_kws["cmap"] = "Blues"
            if "alpha" not in self.kdeplot_kws.keys():          self.kdeplot_kws["alpha"] = 1
            if "shade" not in self.kdeplot_kws.keys():          self.kdeplot_kws["shade"] = True
            if "n_levels" not in self.kdeplot_kws.keys():       self.kdeplot_kws["n_levels"] = 100
            if "shade_lowest" not in self.kdeplot_kws.keys():   self.kdeplot_kws["shade_lowest"] = True
        elif kdeplotEnabled:
            raise Exception ( "The input argument 'kdeplot_kws' must be None or a dictionary,\n"
                            + "each (key,value) pair of which represents an (attribute,value)\n"
                            + "of seaborn library's kdeplot() function."
                            )

        if self.distplot_kws==(): self.distplot_kws={}
        if isinstance(self.distplot_kws,dict):
            if "kde" not in self.distplot_kws.keys(): self.distplot_kws["kde"] = False
            if "hist_kws" not in self.distplot_kws.keys(): self.distplot_kws["hist_kws"] = {"histtype": "stepfilled", "linewidth": 1}
        elif distplotEnabled:
            raise Exception ( "The input argument 'distplot_kws' must be None or a dictionary,\n"
                            + "each (key,value) pair of which represents an (attribute,value)\n"
                            + "of seaborn library's distplot() function."
                            )

        if self.lineplot_kws==(): self.lineplot_kws={}
        if isinstance(self.lineplot_kws,dict):
            self.lineplot_kws["figure_kws"] = None
            if "legend_kws" not in self.lineplot_kws.keys(): self.lineplot_kws["legend_kws"] = None
            if "colorbar_kws" not in self.lineplot_kws.keys(): self.lineplot_kws["colorbar_kws"] = None
            if "lc_kws" not in self.lineplot_kws.keys(): self.lineplot_kws["lc_kws"] = {"cmap":"winter"}
            if "set_kws" not in self.lineplot_kws.keys(): self.lineplot_kws["set_kws"] = self.set_kws
            if "ccolumns" not in self.lineplot_kws.keys():
                if scatterplotEnabled:
                    self.lineplot_kws["ccolumns"] = None
                else:
                    self.lineplot_kws["ccolumns"] = ()
            if self.lineplot_kws["ccolumns"] is None:
                if scatterplotEnabled:
                    if "plot_kws" in self.lineplot_kws.keys():
                        if "color" not in self.lineplot_kws["plot_kws"].keys(): self.lineplot_kws["plot_kws"]["color"] = "grey"
                        if "alpha" not in self.lineplot_kws["plot_kws"].keys(): self.lineplot_kws["plot_kws"]["alpha"] = 0.3
                    else:
                        self.lineplot_kws["plot_kws"] = {"color": "grey", "alpha": 0.3}
            else:
                if "plot_kws" not in self.lineplot_kws.keys(): self.lineplot_kws["plot_kws"] = None
        elif lineplotEnabled:
            raise Exception ( "The input argument 'lineplot_kws' must be None or a dictionary,\n"
                            + "each (key,value) pair of which represents an (attribute,value)\n"
                            + "of ParaMonte library's LinePlot() class."
                            )

        if self.scatterplot_kws==(): self.scatterplot_kws={}
        if isinstance(self.scatterplot_kws,dict):
            self.scatterplot_kws["figure_kws"] = None
            if "legend_kws" not in self.scatterplot_kws.keys(): self.scatterplot_kws["legend_kws"] = None
            if "colorbar_kws" not in self.scatterplot_kws.keys(): self.scatterplot_kws["colorbar_kws"] = None
            if "set_kws" not in self.scatterplot_kws.keys(): self.scatterplot_kws["set_kws"] = self.set_kws
            if "scatter_kws" in self.scatterplot_kws.keys():
                if "s" not in self.scatterplot_kws["scatter_kws"].keys(): self.scatterplot_kws["scatter_kws"]["s"] = 3
                if "cmap" not in self.scatterplot_kws["scatter_kws"].keys(): self.scatterplot_kws["scatter_kws"]["cmap"] = "winter"
                if "alpha" not in self.scatterplot_kws["scatter_kws"].keys(): self.scatterplot_kws["scatter_kws"]["alpha"] = 1
                if "zorder" not in self.scatterplot_kws["scatter_kws"].keys(): self.scatterplot_kws["scatter_kws"]["zorder"] = 4
            else:
                self.scatterplot_kws["scatter_kws"] = {"s":3,"cmap":"winter","alpha":1,"zorder":4}
        elif scatterplotEnabled:
            raise Exception ( "The input argument 'scatterplot_kws' must be None or a dictionary,\n"
                            + "each (key,value) pair of which represents an (attribute,value)\n"
                            + "of ParaMonte library's ScatterPlot() class."
                            )

        if self.colorbar_kws==(): self.colorbar_kws={}

        #########################
        if self._isdryrun: return
        #########################

        # check rows presence

        if self.rows is None:
            rowindex = range(len(self._dfref().index))
        else:
            rowindex = self.rows

        # check columns presence

        colnames, colindex = dfutils.getColNamesIndex(self._dfref().columns,self.columns)

        # make pairgrid plot

        import seaborn as _sns
        from matplotlib import pyplot as _plt

        self.currentFig = ParaMonteFigure()
        if self.set_kws is not None: _sns.set(**self.set_kws)
        self.currentFig.pairgrid = _sns.PairGrid( self._dfref()[colnames].iloc[rowindex,:] )
        self.currentFig.figure = self.currentFig.pairgrid.fig

        # hide axes as requested

        if (lsUpperEnabled or kdeUpperEnabled): 
            self._upperEnabled = True
        else:
            self.hide("upper")
            self._upperEnabled = False

        if (lsLowerEnabled or kdeLowerEnabled):
            self._lowerEnabled = True
        else:
            self.hide("lower")
            self._lowerEnabled = False

        if distplotEnabled: 
            self._diagEnabled = True
        else:
            self._diagEnabled = False
            self.hide("diag")

        # adjust subplot interspaces

        self.currentFig.pairgrid.fig.tight_layout()
        self.currentFig.pairgrid.fig.subplots_adjust(hspace = 0.05, wspace = 0.05)

        # set up figure layout if needed


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

        # define currentFig components

        nrow = len(self.currentFig.pairgrid.axes[:])
        nlsplotStr = str( (lsUpperEnabled+lsLowerEnabled) * nrow * (nrow-1) // 2 )
        nkdeplotStr = str( (kdeUpperEnabled+kdeLowerEnabled) * nrow * (nrow-1) // 2 )

        self.currentFig.targetList      = [ [ [] for j in range(nrow) ] for i in range(nrow) ]
       #self.currentFig.kdeplotList     = [ [ [] for j in range(nrow) ] for i in range(nrow) ]
        self.currentFig.lineplotList    = [ [ [] for j in range(nrow) ] for i in range(nrow) ]
        self.currentFig.scatterplotList = [ [ [] for j in range(nrow) ] for i in range(nrow) ]
       #self.currentFig.distplotList    = [ [] for i in range(nrow) ]

        #####################################
        #### add distplots as needed
        #####################################

        if distplotEnabled:
            self._timer.tic( msg = "adding the diagonal histograms (distplots)... ", end = _eol )
            self.currentFig.pairgrid.map_diag( _sns.distplot, **self.distplot_kws )
            self._timer.toc()

#        if distplotEnabled:
#
#            msg.note( msg = "adding the diagonal histograms (distplots)...", methodName = self._methodName )
#
#            counter = 0
#            for irow, row in enumerate(self.currentFig.pairgrid.axes[:]):
#
#                counter += 1
#                self._timer.tic( msg = "    " + str(counter) + " out of " + nlsplotStr + ": " + colnames[irow] + " - " + colnames[irow] + "... ", end = _eol )
#
#                ax = row[irow]
#                _plt.sca(ax)
#
#                currentXlim = ax.get_xlim()
#                currentYlim = ax.get_ylim()
#
#                self.currentFig.distplotList[irow] = HistPlot   ( dataFrame = self._dfref()
#                                                                , columns = colnames[irow]
#                                                                , rows = self.rows
#                                                                , figure_kws = None
#                                                                , legend_kws = None
#                                                                , set_kws = self.set_kws
#                                                                , distplot_kws = self.distplot_kws
#                                                                )
#                self.currentFig.distplotList[irow].plot()
#
#                if self._lowerEnabled:
#                    ax.set_xlabel("")
#                    ax.set_ylabel("")
#
#                ax.set_xlim(currentXlim)
#                ax.set_ylim(currentYlim)
#
#                self._timer.toc()
#
#            _plt.tight_layout()

        #####################################
        #### add line/scatter plots as needed
        #####################################

        if lsplotEnabled:

            msg.note( msg = "adding the line/scatter plots...", methodName = self._methodName )

            counter = 0
            for irow, row in enumerate(self.currentFig.pairgrid.axes[:]):

                for icol,ax in enumerate(row):

                    self.currentFig.targetList[irow][icol] = Target(axes = ax)

                    if (icol>irow and lsUpperEnabled) or (icol<irow and lsLowerEnabled):

                        counter += 1
                        self._timer.tic( msg = "    " + str(counter) + " out of " + nlsplotStr + ": " + colnames[icol] + " - " + colnames[irow] + "... ", end = _eol )

                        _plt.sca(ax)
                        #_sns.set( self.set_kws )

                        if lineplotEnabled:

                            self.currentFig.lineplotList[irow][icol] =  LinePlot( dataFrame = self._dfref()
                                                                                , xcolumns = colnames[icol]
                                                                                , ycolumns = colnames[irow]
                                                                                , rows = self.rows
                                                                                , **self.lineplot_kws )
                            self.currentFig.lineplotList[irow][icol].plot()

                        if scatterplotEnabled:

                            self.currentFig.scatterplotList[irow][icol] = ScatterPlot   ( dataFrame = self._dfref()
                                                                                        , xcolumns = colnames[icol]
                                                                                        , ycolumns = colnames[irow]
                                                                                        , rows = self.rows
                                                                                        , **self.scatterplot_kws )
                            self.currentFig.scatterplotList[irow][icol].plot()

                        #ax.get_xaxis().set_visible(False)
                        #ax.get_yaxis().set_visible(False)
                        #if irow < nrow-1: ax.set_xlabel("")
                        #if icol > 0: ax.set_ylabel("")
                        #ax.set_xlim(currentXlim)
                        #ax.set_ylim(currentYlim)

                        self._timer.toc()

            _plt.tight_layout()

        #####################################
        #### add kdeplots plots as needed
        #####################################

        if kdeplotEnabled:
            msg.note( msg = "adding kdeplots... depending on the number of plots, this may take a long while.", methodName = self._methodName )
            if kdeUpperEnabled:
                self._timer.tic( msg = "    adding the upper-triangle kdeplots...", end = _eol )
                self.currentFig.pairgrid.map_upper( _sns.kdeplot, **self.kdeplot_kws )
                self._timer.toc()
            if kdeLowerEnabled:
                self._timer.tic( msg = "    adding the lower-triangle kdeplots...", end = _eol )
                self.currentFig.pairgrid.map_lower( _sns.kdeplot, **self.kdeplot_kws )
                self._timer.toc()

        if self.outputFile is not None:
            self.currentFig.figure.savefig  ( self.outputFile
                                            , bbox_inches = 'tight'
                                            , pad_inches = 0.0
                                            )

        # set the ticks and labels

        nrowMinusOne = nrow - 1
        for irow, row in enumerate(self.currentFig.pairgrid.axes[:]):
            for icol, ax in enumerate(row):
                if self._lowerEnabled:
                    if icol > 0:
                        ax.set_ylabel("")
                    if irow < nrowMinusOne:
                        ax.set_xlabel("")
                elif self._upperEnabled:
                    if icol>irow:
                        ax.set_ylabel("")
                        ax.set_xlabel("")

        # add colorbar

        mappable = None
        for scatterplotrow, lineplotrow in zip( self.currentFig.scatterplotList[:], self.currentFig.lineplotList[:] ):
            for scatterplot, lineplot in zip(scatterplotrow, lineplotrow):
                if scatterplot != []:
                    if scatterplot.ccolumns is not None:
                        mappable = scatterplot.currentFig.scatter
                        cbarLabel = scatterplot.ccolumns
                        if len(cbarLabel)==0: cbarLabel = "Count from the Series Start"
                        break
                elif lineplot != []:
                        mappable = lineplot.currentFig.lineCollection
                        cbarLabel = lineplot.ccolumns
                        if len(cbarLabel)==0: cbarLabel = "Count from the Series Start"
                        break

#        sccolumnsEnabled = False
#        if "ccolumns" in self.scatterplot_kws.keys():
#            self.scatterplot_kws["ccolumns"] is not None: sccolumnsEnabled = True
#        if scatterplotEnabled and self.scatterplot_kws["ccolumns"] is not None:
#            if self.currentFig.scatterplotList[0][1] != []:
#                mappable = self.currentFig.scatterplotList[0][1].currentFig.scatter 
#            else:
#                mappable = self.currentFig.scatterplotList[1][0].currentFig.scatter
#            cbarLabel = self.scatterplot_kws["ccolumns"]
#        elif lineplotEnabled and self.lineplot_kws["ccolumns"] is not None:
#            if self.currentFig.lineplotList[0][1] != []:
#                mappable = self.currentFig.lineplotList[0][1].currentFig.lineCollection 
#            else:
#                mappable = self.currentFig.lineplotList[1][0].currentFig.lineCollection
#            cbarLabel = self.lineplot_kws["ccolumns"]
#        #elif kdeplotEnabled and self.kdeplot_kws is not None:

        if mappable is not None:
            self.currentFig.figure.colorbar = self.currentFig.figure.colorbar(mappable=mappable,ax=self.currentFig.pairgrid.axes)
            self.currentFig.figure.colorbar.set_label(cbarLabel)
            self.currentFig.figure.set_figwidth(self.currentFig.figure.get_figwidth()*1.3333)

    ################################################################################################################################

    def hide(self,part="all"):
        """
        hides the requested part of the grid plot.

        Parameters
        ----------
            part
                a string with the following possible values:
                    - "lower": hides the lower triangle of the grid plot.
                    - "upper": hides the lower triangle of the grid plot.
                    - "diag" : hides the diagonal of the grid plot.
                    - "all"  : hides all grid plots.

        Returns
        -------
            None.

        """

        allhidden = part=="all"
        if allhidden or part=="upper": self.currentFig.pairgrid.map_upper(hide_current_axis)
        if allhidden or part=="lower": self.currentFig.pairgrid.map_lower(hide_current_axis)
        if allhidden or part=="diag" : self.currentFig.pairgrid.map_diag(hide_current_axis)

    ################################################################################################################################

    def show(self,part="all"):
        """
        shows the requested part of the grid plot.

        Parameters
        ----------
            part
                a string with the following possible values:
                    - "lower": shows the lower triangle of the grid plot.
                    - "upper": shows the lower triangle of the grid plot.
                    - "diag" : shows the diagonal of the grid plot.
                    - "all"  : shows all grid plots.

        Returns
        -------
            None.

        """

        allshown = part=="all"
        if allshown or part=="upper": self.currentFig.pairgrid.map_upper(show_current_axis)
        if allshown or part=="lower": self.currentFig.pairgrid.map_lower(show_current_axis)
        if allshown or part=="diag" : self.currentFig.pairgrid.map_diag(show_current_axis)

    ################################################################################################################################

    def addTarget   ( self
                    , values        : _tp.Union[ _np.ndarray , _tp.List[float] , _tp.Tuple[float] ] = None
                    , **target_kws
                    ):
        """
        calls target callable associated with each element of 
        currentFig.targetList of the GridPlot object to add target 
        points and/or lines to all plots of the GridPlot figure.
        This method is supposed to called only after a grid plot 
        has been generated.

        Parameters
        ----------
            values
                a numpy array, or list, or tuple whose length equals the 
                number of columns/rows of the grid plot, each element of 
                which represents the target value associated with the 
                corresponding variable on the x axis of the plot, 
                from left to right. 
                If not provided, the default will be set to the state 
                (coordinates) of the "SampleLongFunc" column in the 
                input dataframe to the object. This would work only 
                when the input dataframe is the contents of a ParaMonte 
                output chain or sample file.
            also, any attributes of the Target class.

        Returns
        -------
            None. However, this method causes side-effects by manipulating 
            the existing attributes of the Target objects in targetList.

        """

        diag_target_kws = target_kws.copy()
        diag_target_kws["axhline_kws"] = None
        diag_target_kws["scatter_kws"] = None
        if values is None:
            from _statistics import getMaxLogFunc
            maxLogFunc = getMaxLogFunc(dataFrame=self._dfref())
            values = maxLogFunc.state
        for irow, targetRow in enumerate(self.currentFig.targetList[:]):

            for icol,target in enumerate(targetRow):

                if icol>irow and self._upperEnabled:
                    target( value = [ values[icol], values[irow] ] , **target_kws )
                elif icol<irow and self._lowerEnabled:
                    target( value = [ values[icol], values[irow] ] , **target_kws )
                elif icol==irow and self._diagEnabled:
                    target( value = [ values[icol], values[irow] ] , **diag_target_kws )

    ################################################################################################################################

#!DEC$ endif
