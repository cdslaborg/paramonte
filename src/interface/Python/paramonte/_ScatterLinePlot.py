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
from collections.abc import Iterable as _Iterable

import _dfutils as dfutils

#!DEC$ ifdef PMVIS_ENABLED

#import seaborn as _sns
#from matplotlib import pyplot as _plt
#from matplotlib.collections import LineCollection as _LineCollection

from _visutils import ParaMonteFigure

####################################################################################################################################
#### ScatterLinePlot class
####################################################################################################################################

class ScatterLinePlot:
    """
    This is the ScatterLinePlot class for generating instances 
    of line and/or scatter figures based on matplotlib library's 
    line() and scatter() functions.

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
        xcolumns
            optional argument that determines the columns of dataframe to serve as 
            the x-values. It can have three forms:
                1.  a list of column indices in the input dataFrame.
                2.  a list of column names in dataFrame.columns.
                3.  a range(start,stop,step), representing the column indices in dataFrame.
            Examples:
                1. xcolumns = [0,1,4,3]
                2. xcolumns = ["SampleLogFunc","SampleVariable1"]
                3. xcolumns = range(17,7,-2)
            However, in all cases, it must have a length that is either 1 or equal 
            to the length of ycolumns. If the length is 1, then xcolumns will be 
            plotted against data corresponding to each element of ycolumns.
            If not provided, the default will be the count of the rows of the 
            input dataFrame.
        ycolumns
            optional argument that determines the columns of dataframe to serve 
            as the y-values. It can have three forms:
                1.  a list of column indices in dataFrame.
                2.  a list of column names in dataFrame.columns.
                3.  a range(start,stop,step), representing the column indices in dataFrame.
            Examples:
                1.  ycolumns = [0,1,4,3]
                2.  ycolumns = ["SampleLogFunc","SampleVariable1"]
                3.  ycolumns = range(17,7,-2)
            However, in all cases, it must have a length that is either 1 or 
            equal to the length of xcolumns. If the length is 1, then ycolumns 
            will be plotted against data corresponding to each element of xcolumns. 
            If not provided, the default includes all columns of the dataframe. 
        lccolumns
            (line-color-columns) optional argument that determines the columns 
            of dataframe to serve as the color-values corresponding to each 
            line-segment in the plot. If provided, a matplotlib LineCollection 
            will be added to the plot. It can have three forms:
                1. a list of column indices in dataFrame.
                2. a list of column names in dataFrame.columns.
                3. a range(start,stop,step), representing the column indices in dataFrame.
            Examples:
                1. lccolumns = [0,1,4,3]
                2. lccolumns = ["SampleLogFunc","SampleVariable1"]
                3. lccolumns = range(17,7,-2)
            However, in all cases, it must have a length that is either 1 or 
            equal to the lengths of xcolumns or ycolumns, whichever is not 1. 
            If the length is 1, then the same color will be used for plotting 
            data corresponding to each element of xcolumns. If not provided, 
            the default will be the count of the rows of the input dataFrame. 
            If set to None, no colored LineCollection will be added to the plot.
        sccolumns
            (scatter-color-columns) optional argument that determines the 
            columns of dataframe to serve as the color-values corresponding to 
            each line-segment in the plot. If provided, a line collection will 
            be added to the plot. It can have three forms:
                1. a list of column indices in dataFrame.
                2. a list of column names in dataFrame.columns.
                3. a range(start,stop,step), representing the column indices in dataFrame.
            Examples:
                1. sccolumns = [0,1,4,3]
                2. sccolumns = ["SampleLogFunc","SampleVariable1"]
                3. sccolumns = range(17,7,-2)
            However, in all cases, it must have a length that is either 
            1 or equal to the lengths of xcolumns or ycolumns, whichever is not 1. 
            If the length is 1, then the same color will be used for plotting data 
            corresponding to each element of xcolumns. If not provided, the default 
            will be the count of the rows of the input dataFrame. If set to None, 
            fixed-color scatterplots will be added to the plot.
        rows
            optional argument that determines the rows of dataframe 
            to be visualized. It can be either:
                1. a range(start,stop,step), or, 
                2. a list of row indices in dataFrame.index.
            Examples:
                1. rows = range(17,7,-2)
                2. rows = [i for i in range(7,17)]
            If not provided, the default includes all rows of the dataframe.
        lc_kws
            optional dictionary of keyword arguments to be passed to matplotlib's 
            LineCollection() class. For example: 
                lc_kws = {"cmap":"autumn"}
            The default is {}. 
            If set to None, no colored LineCollection will be plotted.
        set_kws
            optional dictionary of keyword arguments to be passed to seaborn's 
            set() function. For example: 
                set_kws = {"style":"darkgrid"}
            if set to None, then no call to set() function will be made.
        figure_kws
            optional dictionary of keyword arguments to be passed to matplotlib's 
            figure() function. For example: 
                figure_kws = {"facecolor":"w","dpi":150}
        plot_kws
            optional dictionary of keyword arguments to be passed to matplotlib's 
            plot() function. For example: 
                plot_kws = {"linewidth":0.5}
            The default is {}.
            If set to None, neither lines nor colored LineCollection will be plotted. 
        legend_kws
            optional dictionary of keyword arguments to be passed to matplotlib's 
            legend() function. If it is set to None, no legend will be added
            to the plot. If it is set to {}, then the legend will be added 
            to the plot and automatically adjusted. For example: 
                legend_kws = {"labels":["Variable1","Variable2"]}
            legend will be added to plot only if simple line/scatter plots with no 
            color-mappings are requested.
        scatter_kws
            optional dictionary of keyword arguments to be passed to matplotlib's 
            scatter() function. For example: 
                scatter_kws = {"zorder":2}
            The default is {}.
            If set to None, no scatter-plots will be plotted.
        colorbar_kws
            optional dictionary of keyword arguments to be passed to matplotlib's 
            figure.colorbar() function. For example: 
                colorbar_kws = {"orientation":"vertical"}
            The default is {}. The colorbar for scatter-plot is prioritized over
            the colorbar for the colored line-plot.
            If set to None, no colorbar will be plotted.
        outputFile
            optional string representing the name of the output file in which 
            the figure will be saved. If not provided, no output file will be generated.
    
    Attributes
    ----------
        All of the parameters described above, except dataFrame.
            a reference to the dataFrame will be implicitly stored in the object.
        currentFig
            an object of class ParaMonteFigure which is initially None, but upon 
            making a plot, is populated with attributes representing all outputs 
            from matplotlib/seaborn function calls, with the following attributes:
                figure
                    the output of matplotlib's figure() function.
                axes
                    the output of matplotlib's subplot() function.
                plot
                    a list of Line2D objects which is the output 
                    of matplotlib's plot() function.
                scatter
                    an object of type PathCollection which is the 
                    output of matplotlib's scatter() function.
                lineCollection
                    an object of type LineCollection which is the 
                    output of matplotlib's LineCollection() class.
                legend
                    the output of matplotlib's legend() function.
                colorbar
                    the output of matplotlib's Figure.colorbar() function.

    Returns
    -------
        ScatterLinePlot

    ----------------------------------------------------------------------
    """

    def __init__( self
                , dataFrame     : _tp.Optional[ _pd.DataFrame ] = None
                , xcolumns      : _tp.Optional[ _tp.Union[ str, range , _tp.List[int] , _tp.List[str] ] ] = None
                , ycolumns      : _tp.Optional[ _tp.Union[ str, range , _tp.List[int] , _tp.List[str] ] ] = None
                , lccolumns     : _tp.Optional[ _tp.Union[ str, range , _tp.List[int] , _tp.List[str] ] ] = ()
                , sccolumns     : _tp.Optional[ _tp.Union[ str, range , _tp.List[int] , _tp.List[str] ] ] = ()
                , rows          : _tp.Optional[ _tp.Union[ range , _tp.List[int] ] ] = None
                , lc_kws        : _tp.Optional[_tp.Dict] = ()
                , set_kws       : _tp.Optional[_tp.Dict] = ()
                , plot_kws      : _tp.Optional[_tp.Dict] = ()
                , figure_kws    : _tp.Optional[_tp.Dict] = ()
                , legend_kws    : _tp.Optional[_tp.Dict] = ()
                , scatter_kws   : _tp.Optional[_tp.Dict] = ()
                , colorbar_kws  : _tp.Optional[_tp.Dict] = ()
                , outputFile    : _tp.Optional[str] = None
                ):

        self._dfref = None if dataFrame is None else _wref.ref(dataFrame)
        self.xcolumns       = xcolumns
        self.ycolumns       = ycolumns
        self.lccolumns      = lccolumns
        self.sccolumns      = sccolumns
        self.rows           = rows
        self.lc_kws         = lc_kws
        self.set_kws        = set_kws
        self.plot_kws       = plot_kws
        self.figure_kws     = figure_kws
        self.legend_kws     = legend_kws
        self.scatter_kws    = scatter_kws
        self.colorbar_kws   = colorbar_kws
        self.outputFile     = outputFile
        self.currentFig     = None

        self._isdryrun = True
        self.plot()
        self._isdryrun = False

    ################################################################################################################################

    def __call__(self,**kwargs):
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

        reself = False
        for key in kwargs.keys():
            if hasattr(self,key):
                setattr(self, key, kwargs[key])
            elif key=="dataFrame":
                setattr( self, "_dfref", _wref.ref(kwargs[key]) )
            elif key=="reself":
                reself = kwargs["reself"]
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

        # setup figure

        plotEnabled = self.plot_kws is not None
        figEnabled = self.figure_kws is not None

        if figEnabled:
            if self.figure_kws==(): self.figure_kws = {}
            if isinstance(self.figure_kws,dict):
                if "dpi" not in self.figure_kws.keys(): self.figure_kws["dpi"] = 150
                if "facecolor" not in self.figure_kws.keys(): self.figure_kws["facecolor"] = "w"
                if "edgecolor" not in self.figure_kws.keys(): self.figure_kws["edgecolor"] = "w"
            else:
                raise Exception ( "The input argument 'figure_kws' must be a dictionary,\n"
                                + "each (key,value) pair of which represents an (attribute,value)\n"
                                + "of matplotlib library's figure() function."
                                )

        # set what to plot

        scatterEnabled = self.scatter_kws is not None
        lcEnabled = plotEnabled and (self.lc_kws is not None) and (self.lccolumns is not None)
        scEnabled = scatterEnabled and (self.sccolumns is not None)
        lgEnabled = (self.legend_kws is not None) and (not (plotEnabled and lcEnabled)) and (not (scatterEnabled and scEnabled))

        if self.plot_kws==(): self.plot_kws={}
        if isinstance(self.plot_kws,dict):
            if "zorder" not in self.plot_kws.keys(): self.plot_kws["zorder"] = 1
            if "linewidth" not in self.plot_kws.keys(): 
                if lcEnabled:
                    self.plot_kws["linewidth"] = 0
                else:
                    self.plot_kws["linewidth"] = 1
        elif plotEnabled:
            raise Exception ( "The input argument 'plot_kws' must be None or a dictionary,\n"
                            + "each (key,value) pair of which represents an (attribute,value)\n"
                            + "of matplotlib library's plot() function."
                            )

        if self.lc_kws==(): self.lc_kws={}
        if isinstance(self.lc_kws,dict):
            if "cmap" not in self.lc_kws.keys(): self.lc_kws["cmap"] = "autumn"
            if "alpha" not in self.lc_kws.keys(): self.lc_kws["alpha"] = 1
            if "linewidth" not in self.lc_kws.keys(): self.lc_kws["linewidth"] = 1
        elif lcEnabled:
            raise Exception ( "The input argument 'lc_kws' must be None or a dictionary,\n"
                            + "each (key,value) pair of which represents an (attribute,value)\n"
                            + "of matplotlib library's LineCollection() class."
                            )

        if self.scatter_kws==(): self.scatter_kws={}
        if isinstance(self.scatter_kws,dict):
            if "s" not in self.scatter_kws.keys(): self.scatter_kws["s"] = 8
            if "cmap" not in self.scatter_kws.keys(): self.scatter_kws["cmap"] = "autumn"
            if "alpha" not in self.scatter_kws.keys(): self.scatter_kws["alpha"] = 1
            if "edgecolor" not in self.scatter_kws.keys(): self.scatter_kws["edgecolor"] = ""
            if "zorder" not in self.scatter_kws.keys(): self.scatter_kws["zorder"] = 2
            #if scEnabled: self.scatter_kws["c"] = scdata
        elif scatterEnabled:
            raise Exception ( "The input argument 'scatter_kws' must be None or a dictionary,\n"
                            + "each (key,value) pair of which represents an (attribute,value)\n"
                            + "of matplotlib library's scatter() function."
                            )

        #########################

        if self.set_kws==(): self.set_kws={}
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
        from matplotlib.collections import LineCollection as _LineCollection

        # generate figure and axes if needed

        self.currentFig = ParaMonteFigure()
        figEnabled = self.figure_kws is not None
        if figEnabled:
            if self.set_kws is not None: _sns.set(**self.set_kws)
            self.currentFig.figure = _plt.figure( **self.figure_kws )
            self.currentFig.axes = _plt.subplot(1,1,1)
        else:
            self.currentFig.axes = _plt.gca()
            self.currentFig.figure = self.currentFig.axes.get_figure()

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

        if self.xcolumns is None:
            lgxicol = 0
            xcolindex = []
            xcolnames = ["Count"]
            xdata = _np.array( self._dfref().index[rowindex] + 1 ).flatten()
        else:
            xcolnames, xcolindex = dfutils.getColNamesIndex(self._dfref().columns,self.xcolumns)

        ycolnames, ycolindex = dfutils.getColNamesIndex(self._dfref().columns,self.ycolumns)

        # set line plot color data

        if self.lccolumns is None:
            lccolindex = []
            lccolnames = []
        else:
            if isinstance(self.lccolumns, int):
                lccolumnslen = 1
            elif isinstance(self.lccolumns, _Iterable):
                lccolumnslen = len(self.lccolumns)
            else:
                lccolumnslen = 0
            if lccolumnslen==0:
                lccolindex = []
                lccolnames = ["Count"]
                lcdata = _np.array( self._dfref().index[rowindex] + 1 ).flatten()
            else:
                lccolnames, lccolindex = dfutils.getColNamesIndex(self._dfref().columns,self.lccolumns)

        # set scatter plot color data

        if self.sccolumns is None:
            sccolindex = []
            sccolnames = []
        else:
            if isinstance(self.sccolumns, int):
                sccolumnslen = 1
            elif isinstance(self.sccolumns, _Iterable):
                sccolumnslen = len(self.sccolumns)
            else:
                sccolumnslen = 0
            if sccolumnslen==0: 
                sccolindex = []
                sccolnames = ["Count"]
                scdata = _np.array( self._dfref().index[rowindex] + 1 ).flatten()
            else:
                sccolnames, sccolindex = dfutils.getColNamesIndex(self._dfref().columns,self.sccolumns)

        # check the lengths are consistent

        xcolindexlen = len(xcolindex)
        ycolindexlen = len(ycolindex)
        lccolindexlen = len(lccolindex)
        sccolindexlen = len(sccolindex)
        maxLenColumns = _np.max (   [ xcolindexlen
                                    , ycolindexlen
                                    , lccolindexlen
                                    , sccolindexlen
                                    ]
                                )

        if xcolindexlen!=maxLenColumns and xcolindexlen>1: raise Exception("length of xcolumns must be either 1 or equal to the lengths of ycolumns, lccolumns, or sccolumns.")
        if ycolindexlen!=maxLenColumns and ycolindexlen>1: raise Exception("length of ycolumns must be either 1 or equal to the lengths of xcolumns, lccolumns, or sccolumns.")
        if lccolindexlen!=maxLenColumns and lccolindexlen>1: raise Exception("length of lccolumns must be either 1 or equal to the lengths of xcolumns, ycolumns, or sccolumns.")
        if sccolindexlen!=maxLenColumns and sccolindexlen>1: raise Exception("length of sccolumns must be either 1 or equal to the lengths of xcolumns, ycolumns, or lccolumns.")

        # assign data in case of single column assignments

        if xcolindexlen==1:
            lgxicol = 0
            xdata = self._dfref().iloc[rowindex,xcolindex].values.flatten()
        if ycolindexlen==1:
            lgyicol = 0
            ydata = self._dfref().iloc[rowindex,ycolindex].values.flatten()
        if lccolindexlen==1: lcdata = self._dfref().iloc[rowindex,lccolindex].values.flatten()
        if sccolindexlen==1: scdata = self._dfref().iloc[rowindex,sccolindex].values.flatten()

        # add line and/or scatter plots

        if lgEnabled: lglabels = []
        if scEnabled: self.scatter_kws["c"] = scdata
        for i in range(maxLenColumns):

            if xcolindexlen>1:
                lgxicol = i
                xdata = self._dfref().iloc[rowindex,xcolindex[i]].values.flatten()
            if ycolindexlen>1:
                lgyicol = i
                ydata = self._dfref().iloc[rowindex,ycolindex[i]].values.flatten()
            if lccolindexlen>1: lcdata = self._dfref().iloc[rowindex,lccolindex[i]].values.flatten()
            if sccolindexlen>1:
                scdata = self._dfref().iloc[rowindex,sccolindex[i]].values.flatten()
                if scEnabled: self.scatter_kws["c"] = scdata

            if lgEnabled:
                if xcolindexlen<2 and ycolindexlen>1:
                    lglabels.append(ycolnames[lgyicol])
                elif xcolindexlen>1 and ycolindexlen<2:
                    lglabels.append(xcolnames[lgxicol])
                else:
                    lglabels.append(xcolnames[lgxicol]+"-"+ycolnames[lgyicol])

            if plotEnabled:

                self.currentFig.plot = self.currentFig.axes.plot( xdata
                                                                , ydata
                                                                , **self.plot_kws
                                                                )

                # add colored lines if requested

                if lcEnabled:

                    norm = _plt.Normalize(lcdata.min(), lcdata.max())
                    points = _np.array([xdata,ydata]).T.reshape(-1, 1, 2)
                    segments = _np.concatenate([points[:-1], points[1:]], axis=1)
                    lineCollection = _LineCollection( segments
                                                    #, cmap = "winter"
                                                    #, cmap = "viridis"
                                                    #, cmap = "autumn"
                                                    , norm = norm
                                                    , **self.lc_kws
                                                    )
                    lineCollection.set_array(lcdata)
                    #lineCollection.set_linewidth(0.5)
                    self.currentFig.lineCollection = self.currentFig.axes.add_collection(lineCollection)

            if scatterEnabled:

                self.currentFig.scatter = self.currentFig.axes.scatter  ( xdata
                                                                        , ydata
                                                                        , **self.scatter_kws
                                                                        )

        # add line colorbar

        cbarEnabled = (lcEnabled or scEnabled) and (self.colorbar_kws is not None) and (not hasattr(self.currentFig,"colorbar"))
        if cbarEnabled:
            if scEnabled:
                cbar_label = ", ".join(sccolnames)
                cbar_mappable = self.currentFig.scatter
            elif lcEnabled:
                cbar_label = ", ".join(lccolnames)
                cbar_mappable = self.currentFig.lineCollection
            if self.colorbar_kws==(): self.colorbar_kws={}
            self.currentFig.colorbar = self.currentFig.figure.colorbar(mappable = cbar_mappable, ax = self.currentFig.axes, **self.colorbar_kws) #, orientation="vertical")
            self.currentFig.colorbar.set_label(label=cbar_label)

        # add axis labels

        if xcolindexlen>1:
            self.currentFig.axes.set_xlabel("Variable Values")
        else:
            self.currentFig.axes.set_xlabel(xcolnames[0])

        if ycolindexlen>1:
            self.currentFig.axes.set_ylabel("Variable Values")
        else:
            self.currentFig.axes.set_ylabel(ycolnames[0])

        if lgEnabled:
            if self.legend_kws==(): self.legend_kws={}
            if isinstance(self.legend_kws,dict):
                if "labels" not in self.legend_kws.keys(): self.legend_kws["labels"] = lglabels
                self.currentFig.legend = self.currentFig.axes.legend(**self.legend_kws)
            else:
                raise Exception ( "The input argument 'legend_kws' must be a dictionary,\n"
                                + "each (key,value) pair of which represents an (attribute,value)\n"
                                + "of matplotlib library's legend() function."
                                )

        _plt.tight_layout()
        if self.outputFile is not None:
            self.currentFig.figure.savefig  ( self.outputFile
                                            , bbox_inches = 'tight'
                                            , pad_inches = 0.0
                                            )

    ################################################################################################################################

#!DEC$ endif
