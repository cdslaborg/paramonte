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

import _dfutils as dfutils

#!DEC$ ifdef PMVIS_ENABLED

#import seaborn as _sns
#from matplotlib import pyplot as _plt

from _visutils import Target, ParaMonteFigure

####################################################################################################################################
#### HistPlot class
####################################################################################################################################

class HistPlot():
    """
    This is the HistPlot class for generating instances 
    of histogram figures based on seaborn library's distplot().
    Usage: first generate an object of this class by optionally 
    passing the following parameters described below. Then call 
    the plot() method.

    Parameters
    ----------
        dataFrame
            a pandas dataframe containing the data to plot.
        columns
            optional argument that determines the columns of dataframe to be visualized. 
            It can have three forms:
                1. a list of column indices in dataFrame.
                2. a list of column names in dataFrame.columns.
                3. a range(start,stop,step), representing the column indices in dataFrame.
            Examples:
                1. columns = [0,1,4,3]
                2. columns = ["SampleLogFunc","SampleVariable1"]
                3. columns = range(17,7,-2)
            If not provided, the default behavior includes all columns of the dataframe.
        rows
            optional argument that determines the rows of dataframe to be visualized. 
            It can be either:
                1. a range(start,stop,step), or, 
                2. a list of row indices in dataFrame.index.
            Examples:
                1. rows = range(17,7,-2)
                2. rows = [i for i in range(7,17)]
            If not provided, the default behavior includes all rows of the dataframe.
        set_kws
            optional dictionary of keyword arguments to be passed to seaborn's 
            set() function. For example: 
                set_kws = {"style":"darkgrid"}
            if set to None, then no call to set() function will be made.
        figure_kws
            optional dictionary of keyword arguments to be passed to matplotlib's 
            figure() function. For example: 
                figure_kws = {"facecolor":"w","dpi":150}
        legend_kws
            optional dictionary of keyword arguments to be passed to matplotlib's 
            legend() function. If it is not provided, the legend will be added
            to the plot. If it is set to {}, then the legend will be added 
            to the plot and automatically adjusted. For example: 
                legend_kws = {"labels":["Variable1","Variable2"]}
        distplot_kws
            optional dictionary of keyword arguments to be passed to seaborn's 
            distplot() function. For example: 
                distplot_kws =  {"hist_kws":{ "histtype":"step"
                                            , "linewidth":0
                                            , "alpha":1
                                            , "color":"g"
                                            }
                                , "kde":False
                                }
        outputFile
            optional string representing the name of the output file in which 
            the figure will be saved. If not provided, no output file will be generated.
            and the corresponding attribute of object will be None.
        currentFig
            an object of class ParaMonteFigure containing all outputs from 
            matplotlib function calls, with the following attributes:
                figure
                    the output of matplotlib's figure() function.
                axes
                    the output of seaborns's distplot() function.
                legend
                    the output of matplotlib's legend() function.
                outputFile
                    a string value representing the name of the output file to which 
                    the figure has been written. 

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
                axes
                    the output of seaborns's heatmap() function.
                legend
                    the output of matplotlib's legend() function.

    Returns
    -------
        HistPlot

    ----------------------------------------------------------------------
    """

    def __init__( self
                , dataFrame         : _tp.Optional[ _pd.DataFrame ] = None
                , columns           : _tp.Optional[ _tp.Union[ str, range , _tp.List[int] , _tp.List[str] ] ] = None
                , rows              : _tp.Optional[ _tp.Union[ range , _tp.List[int] ] ] = None
                , set_kws           : _tp.Optional[_tp.Dict] = ()
                , figure_kws        : _tp.Optional[_tp.Dict] = ()
                , legend_kws        : _tp.Optional[_tp.Dict] = ()
                , distplot_kws      : _tp.Optional[_tp.Dict] = ()
                , outputFile        : _tp.Optional[str] = None
                ):

        self._dfref = None if dataFrame is None else _wref.ref(dataFrame)
        self.columns        = columns
        self.rows           = rows
        self.set_kws        = set_kws
        self.figure_kws     = figure_kws
        self.legend_kws     = legend_kws
        self.distplot_kws   = distplot_kws
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
        Generate a histogram from the selected columns of object's dataframe.

        Parameters
        ----------
            None

        Returns
        -------
            None. However, this method causes side-effects by manipulating 
            the existing attributes of the object.

        """

        if self.figure_kws==(): self.figure_kws={}
        if self.figure_kws is not None:
            if isinstance(self.figure_kws,dict):
                if "dpi" not in self.figure_kws.keys(): self.figure_kws["dpi"] = 150
                if "facecolor" not in self.figure_kws.keys(): self.figure_kws["facecolor"] = "w"
                if "edgecolor" not in self.figure_kws.keys(): self.figure_kws["edgecolor"] = "w"
            else:
                raise Exception ( "The input argument 'figure_kws' must be a dictionary,\n"
                                + "each (key,value) pair of which represents an (attribute,value)\n"
                                + "of matplotlib library's figure() function."
                                )

        if self.distplot_kws==(): self.distplot_kws={}
        if isinstance(self.distplot_kws,dict):
            if "kde" not in self.distplot_kws.keys(): self.distplot_kws["kde"] = False
            if "hist_kws" not in self.distplot_kws.keys(): self.distplot_kws["hist_kws"] = dict()
            if "linewidth" not in self.distplot_kws["hist_kws"].keys(): self.distplot_kws["hist_kws"]["linewidth"] = 0
            if "histtype" not in self.distplot_kws["hist_kws"].keys(): self.distplot_kws["hist_kws"]["histtype"] = "stepfilled"
            self.distplot_kws["hist_kws"]["density"] = self.distplot_kws["kde"]
        else:
            raise Exception ( "The input argument 'distplot_kws' must be a dictionary,\n"
                            + "each (key,value) pair of which represents an (attribute,value)\n"
                            + "of matplotlib library's figure() function."
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

        # check columns presence

        colnames, colindex = dfutils.getColNamesIndex(self._dfref().columns,self.columns)

        # check rows presence

        if self.rows is None:
            rowindex = range(len(self._dfref().index))
        else:
            rowindex = self.rows

        selectedRows = self._dfref().index[rowindex]
        for icol in colindex:
            self.currentFig.axes =_sns.distplot ( self._dfref().iloc[selectedRows,icol]
                                                , **self.distplot_kws
                                                #, kde = False
                                                #, hist_kws = {"histtype": "step", "linewidth": 1, "alpha": 1}#, "color": "g"}
                                                )

        if self.distplot_kws["kde"]:
            self.currentFig.axes.set_ylabel("Probability Density")
        else:
            self.currentFig.axes.set_ylabel("Count")

        if len(colindex)==1:
            self.currentFig.axes.set_xlabel(colnames[0])
        else:
            self.currentFig.axes.set_xlabel("Value")

        if self.legend_kws is not None:
            if self.legend_kws==(): self.legend_kws={}
            if isinstance(self.legend_kws,dict):
                if self.legend_kws:
                    if "labels" not in self.legend_kws.keys(): self.legend_kws["labels"] = tuple(colnames)
                    self.currentFig.legend = self.currentFig.axes.legend(**legend_kws)
                else:
                    self.currentFig.legend = self.currentFig.axes.legend(labels=tuple(colnames))
            else:
                raise Exception ( "The input argument 'legend_kws' must be a dictionary,\n"
                                + "each (key,value) pair of which represents an (attribute,value)\n"
                                + "of matplotlib library's legend() function."
                                )

        _plt.tight_layout()

        if self.outputFile is not None: 
            self.currentFig.axes.get_figure().savefig   ( self.outputFile
                                                        , bbox_inches = 'tight'
                                                        , pad_inches = 0.0
                                                        )

    ################################################################################################################################

#!DEC$ endif