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

import typing as tp
import weakref as wref

import _paramonte as pm
import _pmutils as pmutils
from paramonte.vis._BasePlot import BasePlot

Struct = pmutils.Struct
newline = pmutils.newline

####################################################################################################################################
#### HeatMapPlot class
####################################################################################################################################

class HeatMapPlot(BasePlot):
    """

    This is the HeatMapPlot class for generating instances of 
    heatmap figures based on the seaborn library's ``heatmap()``.

        **Usage**

            First generate an object of this class by optionally
            passing the following parameters described below. 
            Then call the ``make()`` method.

        **Parameters**

            plotType

                A string indicating the name of the plot to be constructed.

            dataFrame (optional)

                A Pandas dataFrame whose data will be plotted.

            methodName (optional)

                The name of the ParaMonte sample requesting the BasePlot.

            reportEnabled (optional)

                A boolean whose value indicates whether guidelines 
                should be printed in the standard output.

            resetPlot (optional)

                A function that resets the properties of the plot as desired 
                from outside. If provided, a pointer to this function will be
                saved for future internal usage.

        **Attributes**

            columns

                An attribute that determines the columns of dataFrame 
                to be visualized. It can have three forms:

                    1.  A list of column indices in dataFrame.
                    2.  A list of column names in dataFrame.columns.
                    3.  A ``range(start,stop,step)`` of column indices.

                Examples:

                    1.  ``xcolumns = [0,1,4,3]``
                    2.  ``xcolumns = ["SampleLogFunc","SampleVariable1"]``
                    3.  ``xcolumns = range(17,7,-2)``

                The default behavior includes all columns of the dataFrame.

            rows

                An attribute that determines the rows of dataFrame to be 
                visualized. It can be either:

                    1.  A ``range(start,stop,step)``, or, 
                    2.  A list of row indices in dataFrame.index.

                Examples:

                    1.  ``rows = range(17,7,-2)``
                    2.  ``rows = [i for i in range(7,17)]``

                The default behavior includes all rows of the dataFrame.

            heatmap

                A structure with one attribute:

                    kws

                        A structure whose components are directly passed as 
                        keyword arguments to the ``heatmap()`` function of 
                        the seaborn library.

                        Example usage:

                            .. code-block:: python

                                heatmap.square = True

                        **NOTE**

                        If a desired property is missing among the ``kws`` 
                        attributes, simply add the field and its value to 
                        the component.

            annotPrecision

                optional integer which determines the precision with which 
                 the output numbers corresponding to each pixel are written 
                 to the heatmap when annotation is requested via 
                 ``heatmap.kws.annot = True``. The default is 2, meaning 
                 that only two digits after the decimal are considered.

            xticklabels

                A structure with one attribute:

                    kws

                        A structure whose components are directly passed as 
                        keyword arguments to the ``set_xticklabels()`` 
                        function of the matplotlib library.

                        Example usage:

                            .. code-block:: python

                                xticklabels.verticalalignment = "right"

                        **NOTE**

                        If a desired property is missing among the ``kws`` 
                        attributes, simply add the field and its value to 
                        the component.

            yticklabels

                A structure with one attribute:

                    kws

                        A structure whose components are directly passed as 
                        keyword arguments to the ``set_yticklabels()`` 
                        function of the matplotlib library.

                        Example usage:

                            .. code-block:: python

                                yticklabels.verticalalignment = "right"

                        **NOTE**

                        If a desired property is missing among the ``kws`` 
                        attributes, simply add the field and its value to 
                        the component.

            set

                A structure with two attributes:

                    enabled

                        A boolean indicating whether a call to the ``set()``
                        function of the seaborn library should be made or not.

                    kws

                        A structure whose components are directly passed as 
                        keyword arguments to the ``set()`` function.

                        Example usage:

                            .. code-block:: python

                                set.kws.style = "darkgrid"

                        **NOTE**

                        If a desired property is missing among the ``kws`` 
                        attributes, simply add the field and its value to 
                        the component.

            axes (available only in 1D and 2D plots)

                A structure with one attribute:

                    kws

                        A structure whose components are directly passed as 
                        keyword arguments to the ``gca()`` function of the 
                        matplotlib library.

                        Example usage:

                            .. code-block:: python

                                axes.kws.facecolor = "w"

                        **NOTE**

                        If a desired property is missing among the ``kws`` 
                        attributes, simply add the field and its value to 
                        the component.

            figure

                A structure with two attributes:

                    enabled

                        A boolean indicating whether a call to the ``figure()``
                        function of the matplotlib library should be made or not.
                        If a call is made, a new figure will be generated.
                        Otherwise, the current active figure will be used.

                    kws

                        A structure whose components are directly passed as 
                        keyword arguments to the ``figure()`` function.

                        Example usage:

                            .. code-block:: python

                                figure.kws.facecolor = "w"

                        **NOTE**

                        If a desired property is missing among the ``kws`` 
                        attributes, simply add the field and its value to 
                        the component.

            currentFig

                A structure whose attributes are the outputs of various plotting 
                tools used to make the current figure. These include the handle 
                to the current figure, the handle to the current axes in the plot, 
                the handle to the colorbar (if any exists), and other Python 
                plotting tools used to make to generate the figure.

        **Returns**

                An object of class ``HeatMapPlot``.

    ---------------------------------------------------------------------------
    """

    ################################################################################################################################
    #### __init__
    ################################################################################################################################

    def __init__( self
                , plotType      # : str
                , dataFrame     = None # : tp.Optional[ pd.DataFrame ] = None
                , methodName    = "ParaMonte" # : tp.Optional[ str ] = "ParaMonte"
                , reportEnabled = True # : tp.Optional[ bool ] = True
                , resetPlot     = None
                ):

        super().__init__( plotType = plotType
                        , dataFrame = dataFrame
                        , methodName = methodName
                        , reportEnabled = reportEnabled
                        , resetPlot = resetPlot
                        )

        self._reset()
        if resetPlot is None: self._resetPlot = self._reset
        self._progress.note()

    ################################################################################################################################
    #### _reset
    ################################################################################################################################

    def _reset(self):

        super()._reset()

        self.heatmap = Struct()
        self.heatmap.enabled = True
        self.heatmap.kws = Struct()

        self.xticklabels = Struct()
        self.xticklabels.kws = Struct()

        self.yticklabels = Struct()
        self.yticklabels.kws = Struct()

        self.columns = ""

        self._colorStart    = 20
        self._colorCount    = 200
        self._colorEnd      = 220

        self.annotPrecision = 2

        self._isdryrun = True
        self.make()
        self._isdryrun = False

    ################################################################################################################################
    #### __call__
    ################################################################################################################################

    def __call__( self
                , reself : tp.Optional[ bool ] = False
                , **kwargs
                ):
        """

        Call the ``make()`` method of the 
        current instance of the class.

            **Parameters**

                Any arguments that can be passed to the 
                ``make()`` method of the plot object.

            **Returns**

                Any return value from the ``make()`` method of the plot object.

        """

        return self.make(reself, **kwargs)

    ################################################################################################################################
    #### make
    ################################################################################################################################

    def make( self
            , reself : tp.Optional[ bool ] = False
            , **kwargs
            ):
        """

        Generate a heatmap plot from the selected 
        columns of the object's dataframe.

            **Parameters**

                reself

                    A logical variable. If ``True``, an instance of the 
                    object will be returned upon exit to the calling 
                    routine. The default value is ``False``.

            **Returns**

                the object self if ``reself = True`` otherwise, ``None``.
                However, this method causes side-effects by manipulating
                the existing attributes of the object.

        """

        for key in kwargs.keys():
            if hasattr(self,key):
                setattr(self, key, kwargs[key])
            elif key=="dataFrame":
                setattr( self, "_dfref", wref.ref(kwargs[key]) )
            else:
                raise Exception ( newline
                                + "Unrecognized input '"+key+"' class attribute detected." + newline
                                + self._getDocString()
                                )

        # set what to plot

        ############################################################################################################################
        #### xticklabels / yticklabels properties
        ############################################################################################################################

        if isinstance(self.xticklabels.kws,Struct):
            if "horizontalalignment" not in vars(self.xticklabels.kws).keys(): self.xticklabels.kws.horizontalalignment = "right"
            if "rotation" not in vars(self.xticklabels.kws).keys(): self.xticklabels.kws.rotation = 45
        else:
            raise Exception ( newline
                            + "The xticklabels.kws component of the current HeatMapPlot object must" + newline
                            + "be an object of class Struct(), essentially a structure with components" + newline
                            + "whose names are the input arguments to the set_xticklabels() method of the" + newline
                            + "Axes class of the matplotlib library." + newline
                            + self._getDocString()
                            )

        if isinstance(self.yticklabels.kws,Struct):
            if "horizontalalignment" not in vars(self.yticklabels.kws).keys(): self.yticklabels.kws.horizontalalignment = "right"
            if "rotation" not in vars(self.yticklabels.kws).keys(): self.yticklabels.kws.rotation = 45
        else:
            raise Exception ( newline
                            + "The yticklabels.kws component of the current HeatMapPlot object must" + newline
                            + "be an object of class Struct(), essentially a structure with components" + newline
                            + "whose names are the input arguments to the set_yticklabels() method of the" + newline
                            + "Axes class of the matplotlib library." + newline
                            + self._getDocString()
                            )

        ############################################################################################################################
        #### heatmap properties
        ############################################################################################################################

        if isinstance(self.heatmap.kws,Struct):
            if "square" not in vars(self.heatmap.kws).keys(): self.heatmap.kws.square = True
            if "cmap" not in vars(self.heatmap.kws).keys():
                try:
                    import seaborn as sns
                    self.heatmap.kws.cmap = sns.diverging_palette   ( h_neg = self._colorStart
                                                                    , h_pos = self._colorEnd
                                                                    , n = self._colorCount
                                                                    )
                except:
                    if self._isdryrun:
                        self.heatmap.kws.cmap = None
                    else:
                        raise Exception ( newline
                                        + "Failed to set the heatmap.kws.cmap component of the current HeatMapPlot object." + newline
                                        + "This component depends on the external seaborn Python library. Therefore, it is " + newline
                                        + "likely that the seaborn library or one of the required components of it, such as " + newline
                                        + "the matplotlib Python library is not properly installed on your system. Please " + newline
                                        + "fix this issue, otherwise, the visualization tools of the ParaMonte library " + newline
                                        + "will not work as expected. You can install the seaborn library by typing " + newline
                                        + "the following commands in your Anaconda3 or Bash command prompt: " + newline
                                        + newline
                                        + "    pip install --user --upgrade matplotlib"
                                        + "    pip install --user --upgrade seaborn"
                                        + newline
                                        + self._getDocString()
                                        )
        else:
            raise Exception ( newline
                            + "The heatmap.kws component of the current HeatMapPlot object must" + newline
                            + "be an object of class Struct(), essentially a structure with components" + newline
                            + "whose names are the input arguments to the heatmap() function of the" + newline
                            + "seaborn library." + newline
                            + self._getDocString()
                            )

        ############################################################################################################################
        #### figure properties
        ############################################################################################################################

        if self.figure.enabled:
            if isinstance(self.figure.kws, Struct):
                if "dpi" not in vars(self.figure.kws).keys(): self.figure.kws.dpi = 150
                if "facecolor" not in vars(self.figure.kws).keys(): self.figure.kws.facecolor = "w"
                if "edgecolor" not in vars(self.figure.kws).keys(): self.figure.kws.edgecolor = "w"
            else:
                raise Exception ( newline
                                + "The figure.kws component of the current DensityPlot object must" + newline
                                + "be an object of class Struct(), essentially a structure with components" + newline
                                + "whose names are the input arguments to the figure() function of the" + newline
                                + "matplotlib library." + newline
                                + self._getDocString()
                                )

        ############################################################################################################################
        ############################################################################################################################
        if self._isdryrun: return
        ############################################################################################################################
        ############################################################################################################################

        import seaborn as sns
        import matplotlib.pyplot as plt

        ############################################################################################################################
        #### generate figure and axes if needed
        ############################################################################################################################

        self._constructBasePlot()

        ############################################################################################################################
        #### check data type
        ############################################################################################################################

        self._checkDataType()

        ############################################################################################################################
        #### check rows presence. This must be checked here, because it depends on the integrity of the in input dataFrame.
        ############################################################################################################################

        if self.rows is None: self.rows = range(len(self._dfref().index))
        rownames = self._dfref().index[self.rows]

        ############################################################################################################################
        #### check columns presence. This must be checked here, because it depends on the integrity of the in input dataFrame.
        ############################################################################################################################

        colnames, colindex = pm.dfutils.getColNamesIndex(self._dfref().columns,self.columns)

        ############################################################################################################################
        #### set up tick labels
        ############################################################################################################################

        xtickExists = True
        if "xticklabels" in vars(self.heatmap.kws).keys():
            if not any(self.heatmap.kws.xticklabels): xtickExists = False
        else:
            self.heatmap.kws.xticklabels = colnames

        ytickExists = True
        if "yticklabels" in vars(self.heatmap.kws).keys():
            if not any(self.heatmap.kws.yticklabels): ytickEyists = False
        else:
            self.heatmap.kws.yticklabels = colnames

        ############################################################################################################################
        #### plot data
        ############################################################################################################################

        if self.annotPrecision is None:
            data = self._dfref().iloc[self.rows,colindex]
        else:
            data = self._dfref().iloc[self.rows,colindex].round( decimals = self.annotPrecision )

        self.currentFig.axes = sns.heatmap  ( data = data
                                            , **vars(self.heatmap.kws)
                                            )

        ############################################################################################################################
        #### configure the tick labels (orientation, ...)
        ############################################################################################################################

        self.currentFig.axes.set_xticklabels( self.currentFig.axes.get_xticklabels()
                                            , **vars(self.xticklabels.kws)
                                            );

        self.currentFig.axes.set_yticklabels( self.currentFig.axes.get_yticklabels()
                                            , **vars(self.yticklabels.kws)
                                            );

        plt.tight_layout()

        if self.figure.enabled:

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

        ############################################################################################################################

        if reself: return self

    ################################################################################################################################
    #### _getDocString
    ################################################################################################################################

    def _getDocString(self):
        docString   = newline \
                    + "Here is the help information on the HeatMapPlot class:" + newline \
                    + newline \
                    + self.__doc__ \
                    + newline \
                    + "Here is the help information on the parent BasePlot class:" + newline \
                    + newline \
                    + super().__doc__
        return docString

    ################################################################################################################################
    #### helpme
    ################################################################################################################################

    def helpme(self, topic=None):
        """

        Print the documentation for the input string topic. 
        If the topic does not exist, the documentation for
        the object will be printed.
     
            **Parameters**
         
                topic (optional)

                A string containing the name of the object 
                for which help is needed.

            **Returns**
         
                None
         
            **Example**

                .. code-block:: python
                    :linenos:

                    helpme()
                    helpme("make")
                    helpme("helpme")
                    helpme("getLogLinSpace")

        """

        try:
            exec("print(self."+topic+".__doc__)")
        except:
            print(self._getDocString())
        return None

    ################################################################################################################################
