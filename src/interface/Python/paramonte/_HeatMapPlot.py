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
####   we ask you to acknowledge the use of the ParaMonte library
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

from _visutils import ParaMonteFigure

####################################################################################################################################
#### HeatMapPlot class
####################################################################################################################################

class HeatMapPlot():
    """

    .. py:class:: HeatMapPlot

    This is the HeatMapPlot class for generating instances
    of histogram figures based on seaborn library's ``heatmap()``.

        **Usage**

            First generate an object of this class by optionally
            passing the following parameters described below. Then call
            the ``plot()`` method.

        **Parameters**

            dataFrame
                a numpy 2D-array or a pandas dataframe containing the data to plot.

            columns
                optional argument that determines the columns of dataframe to be visualized.
                It can have three forms:

                    1. a list of column indices in the input dataFrame.
                    2. a list of column names in dataFrame.columns.
                    3. a ``range(start,stop,step)``, representing the column indices in dataFrame.

                Examples:

                    1. ``columns = [0,1,4,3]``
                    2. ``columns = ["SampleLogFunc","SampleVariable1"]``
                    3. ``columns = range(17,7,-2)``

                If not provided, the default behavior includes all columns of the dataframe.

            rows
                optional argument that determines the rows of dataframe to be visualized.
                It can be either:
                    1. a ``range(start,stop,step)``, or,
                    2. a list of row indices in dataFrame.index.

                Examples:

                    1. ``rows = range(17,7,-2)``
                    2. ``rows = [i for i in range(7,17)]``

                If not provided, the default behavior includes all rows of the dataframe.

            set_kws
                optional dictionary of keyword arguments to be passed to seaborn's
                ``set()`` function. For example:

                .. code-block:: python

                    set_kws = {"style":"darkgrid"}

                if set to ``None``, then no call to ``set()`` function will be made.

            figure_kws
                optional dictionary of keyword arguments to be passed to matplotlib's
                ``figure()`` function. For example:

                .. code-block:: python

                    figure_kws = {"facecolor":"w","dpi":150}

            legend_kws
                optional dictionary of keyword arguments to be passed to matplotlib's
                ``legend()`` function. If it is not provided, the legend will be added
                to the plot. If it is set to ``{}``, then the legend will be added
                to the plot and automatically adjusted. For example:

                .. code-block:: python

                    legend_kws = {"labels":["Variable1","Variable2"]}

            heatmap_kws
                optional dictionary of keyword arguments to be passed to seaborn's
                ``heatmap()`` function. For example:

                .. code-block:: python

                    heatmap_kws = {"square": True}

            xticklabels_kws
                optional dictionary of keyword arguments to be passed to matplotlib's
                ``set_xticklabels()`` method of the Axes object of the plot. For example:

                .. code-block:: python

                    xticklabels_kws = {"verticalalignment": "baseline"}

            yticklabels_kws
                optional dictionary of keyword arguments to be passed to matplotlib's
                ``set_yticklabels()`` method of the Axes object of the plot. For example:

                .. code-block:: python

                    yticklabels_kws = {"verticalalignment": "baseline"}

            annotPrecision
                optional integer which determines the precision of the output numbers
                when heatmap annotation is requested via ``'annot'`` key in heatmap_kws.

            outputFile
                optional string representing the name of the output file in which
                the figure will be saved. If not provided, no output file will be generated.
                and the corresponding attribute of object will be ``None``.

        **Attributes**

            All of the parameters described above, except dataFrame.
                a reference to the dataFrame will be implicitly stored in the object.

            currentFig
                an object of class ParaMonteFigure which is initially None, but upon
                making a plot, is populated with attributes representing all outputs
                from matplotlib/seaborn function calls, with the following attributes:

                    figure
                        the output of matplotlib's ``figure()`` function.

                    axes
                        the output of seaborns's ``heatmap()`` function.

        **Returns**

            self
                an object of class ``HeatMapPlot``.

    ---------------------------------------------------------------------------
    """

    def __init__( self
               #, dataFrame         : _tp.Optional[ _tp.Union[_np.array,_pd.DataFrame] ] = None
                , dataFrame         : _tp.Optional[_pd.DataFrame] = None
                , columns           : _tp.Optional[ _tp.Union[ str, range , _tp.List[int] , _tp.List[str] ] ] = None
                , rows              : _tp.Optional[ _tp.Union[ range , _tp.List[int] ] ] = None
                , set_kws           : _tp.Optional[_tp.Dict] = ()
                , figure_kws        : _tp.Optional[_tp.Dict] = ()
                , legend_kws        : _tp.Optional[_tp.Dict] = ()
                , heatmap_kws       : _tp.Optional[_tp.Dict] = ()
                , xticklabels_kws   : _tp.Optional[_tp.Dict] = ()
                , yticklabels_kws   : _tp.Optional[_tp.Dict] = ()
                , annotPrecision    : _tp.Optional[int] = None
                , outputFile        : _tp.Optional[str] = None
                ):

        self._dfref             = None if dataFrame is None else _wref.ref(dataFrame)
        self.columns            = columns
        self.rows               = rows
        self.set_kws            = set_kws
        self.figure_kws         = figure_kws
        self.legend_kws         = legend_kws
        self.heatmap_kws        = heatmap_kws
        self.xticklabels_kws    = xticklabels_kws
        self.yticklabels_kws    = yticklabels_kws
        self.annotPrecision     = annotPrecision
        self.outputFile         = outputFile
        self.currentFig         = None

        self._colorStart        = 20
        self._colorCount        = 200
        self._colorEnd          = 220

        self._isdryrun = True
        self.plot()
        self._isdryrun = False

    ################################################################################################################################

    def __call__( self
                , reself    : _tp.Optional[ bool ] = False
                , **kwargs
                ):
        """

        .. py:method:: __call__(self, reself = False, **kwargs)

        calls the ``plot()`` method of the current instance of the class.

            **Parameters**

                reself
                    logical variable. If True, an instance of the object
                    will be returned upon exit to the calling routine.
                    The default value is ``False``.
                also, any attributes of the current instance of the class.

            **Returns**

                the object self if ``reself = True`` otherwise, None.
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

        .. py:method:: plot(self)

        Generate a heatmap of the requested columns from the selected columns of
        object's dataframe.

            **Parameters**

                None

            **Returns**

                None. However, this method causes side-effects by manipulating
                the existing attributes of the object.

        """

        if self.xticklabels_kws==(): self.xticklabels_kws={}
        if isinstance(self.xticklabels_kws,dict):
            if "horizontalalignment" not in self.xticklabels_kws.keys(): self.xticklabels_kws["horizontalalignment"] = "right"
            if "rotation" not in self.xticklabels_kws.keys(): self.xticklabels_kws["rotation"] = 45
        else:
            raise Exception ( "The input argument 'xticklabels_kws' must be a dictionary,\n"
                            + "each (key,value) pair of which represents an (attribute,value)\n"
                            + "of matplotlib library's Axes.set_xticklabels() method."
                            )

        if self.yticklabels_kws==(): self.yticklabels_kws={}
        if isinstance(self.yticklabels_kws,dict):
            if "horizontalalignment" not in self.yticklabels_kws.keys(): self.yticklabels_kws["horizontalalignment"] = "right"
            if "rotation" not in self.yticklabels_kws.keys(): self.yticklabels_kws["rotation"] = 45
        else:
            raise Exception ( "The input argument 'yticklabels_kws' must be a dictionary,\n"
                            + "each (key,value) pair of which represents an (attribute,value)\n"
                            + "of matplotlib library's Axes.set_yticklabels() method."
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

        # moved here from the top because of dependency on seaborn _sns

        if self.heatmap_kws==(): self.heatmap_kws={}
        if isinstance(self.heatmap_kws,dict):
            if "square" not in self.heatmap_kws.keys(): self.heatmap_kws["square"] = True
            if "cmap" not in self.heatmap_kws.keys():
                self.heatmap_kws["cmap"] = _sns.diverging_palette   ( h_neg = self._colorStart
                                                                    , h_pos = self._colorEnd
                                                                    , n = self._colorCount
                                                                    )
        else:
            raise Exception ( "The input argument 'heatmap_kws' must be a dictionary,\n"
                            + "each (key,value) pair of which represents an (attribute,value)\n"
                            + "of Seaborn library's heatmap() function."
                            )

        # generate figure and axes if needed

        self.currentFig = ParaMonteFigure()
        if self.set_kws is not None: _sns.set(**self.set_kws)
        figEnabled = self.figure_kws is not None
        if figEnabled:
            self.currentFig.figure = _plt.figure( **self.figure_kws )
            self.currentFig.axes = _plt.subplot(1,1,1)
        else:
            self.currentFig.axes = _plt.gca()
            self.currentFig.figure = self.currentFig.axes.get_figure()

        # check data type

        noDataPassed = ""
        if self._dfref is None:
            noDataPassed = "It appears that no data has been passed for plotting.\n"
        elif not isinstance(self._dfref,_wref.ref):
            noDataPassed = "It appears that you have messed with the\ninternal representation of data in the object.\n"
        elif not isinstance(self._dfref(),_pd.DataFrame):
            noDataPassed = "The input data is not a pandas' dataframe.\n"
        if noDataPassed:
            raise Exception ( noDataPassed
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
        rownames = self._dfref().index[rowindex]

        # set up tick labels

        xtickExists = True
        if "xticklabels" in self.heatmap_kws.keys():
            if not any(self.heatmap_kws["xticklabels"]): xtickExists = False
        else:
            self.heatmap_kws["xticklabels"] = colnames

        ytickExists = True
        if "yticklabels" in self.heatmap_kws.keys():
            if not any(self.heatmap_kws["yticklabels"]): ytickExists = False
        else:
            self.heatmap_kws["yticklabels"] = rownames

        # plot data

        if self.annotPrecision is None:
            self.currentFig.axes = _sns.heatmap ( data = self._dfref().iloc[rowindex,colindex]
                                                , **self.heatmap_kws
                                                )
        else:
            self.currentFig.axes =_sns.heatmap  ( data = self._dfref().iloc[rowindex,colindex].round( decimals = self.annotPrecision )
                                                , **self.heatmap_kws
                                                )

        # configure the tick labels (orientation, ...)

        self.currentFig.axes.set_xticklabels( self.currentFig.axes.get_xticklabels()
                                            , **self.xticklabels_kws
                                            );

        self.currentFig.axes.set_yticklabels( self.currentFig.axes.get_yticklabels()
                                            , **self.yticklabels_kws
                                            );

        _plt.tight_layout()

        if figEnabled:

            # default figure size

            figWidth = 6.4  # inches
            figHeight = 4.8 # inches
            figWidthScale = 1
            figHeightScale = 1
            threshDimension = 10

            # scale only if ticklabels are present

            if xtickExists:
                figWidthScale = max( 1 , self._dfref().shape[1]/threshDimension )
                figWidth *= figWidthScale

            if ytickExists:
                figHeightScale = max( 1 , self._dfref().shape[0]/threshDimension )
                figHeight *= figHeightScale

            self.currentFig.figure.set_size_inches(figWidth,figHeight)

        if self.outputFile is not None:
            self.currentFig.figure.savefig  ( self.outputFile
                                            , bbox_inches = 'tight'
                                            , pad_inches = 0.0
                                            )

    ################################################################################################################################

#!DEC$ endif