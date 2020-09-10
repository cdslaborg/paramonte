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
import _pmutils as pmutils

Struct = pmutils.Struct
newline = pmutils.newline

####################################################################################################################################
#### Target class
####################################################################################################################################

class Target:
    """

    This is the Target class for generating instances
    of a target on the current active figure and axis.

        **Usage**

        First generate an object of this class by optionally
        passing the following parameters described below. Then call
        the ``make()`` method. The generated object is also callable 
        with the same input parameters as the object's constructor.

        **Parameters**

            values (optional)

                A pair (list, array, or tuple) of floats, 
                representing the (x,y) coordinates of the target.

            axes

                The axes object to which the target must be added.
                The default is ``None`` in which case the output of the 
                matplotlib's ``gca()`` function will be used to get the 
                current active axes.

        **Attributes**

            All of the parameters described above, except ``axes``.

                The input axes object (whether user-provided or automatically fetched)
                will be stored as a component of the object's attribute ``currentFig``.

            axvline (optional)

                A structure with two components:

                    enabled

                        A logical value, indicating whether an ``axvline`` 
                        will be added to the plot or not. If set to ``False``
                        no vertical target line will be added.

                    kws

                        A structure whose components will be passed to the 
                        matplotlib's ``axvline()`` function. For example:

                        .. code-block:: python
                            :linenos:

                            axvline.kws.linewidth = 0.5
                            axvline.kws.color = "orangered"
                            axvline.kws.linestyle = ":"
                            axvline.kws.ymin = 0.1
                            axvline.kws.ymax = 1.9

                        **NOTE**

                        If a desired property is missing in the 
                        structure, simply add the property and 
                        its value to the structure.

            axhline (optional)

                A structure with two components:

                    enabled

                        A logical value, indicating whether an ``axhline`` 
                        will be added to the plot or not. If set to ``False``
                        no horizontal target line will be added.

                    kws

                        A structure whose components will be passed to the 
                        matplotlib's ``axhline()`` function. For example:

                        .. code-block:: python
                            :linenos:

                            axhline.kws.linewidth = 0.5
                            axhline.kws.color = "orangered"
                            axhline.kws.linestyle = ":"
                            axhline.kws.xmin = 0.1
                            axhline.kws.xmax = 1.9

                        **NOTE**

                        If a desired property is missing in the 
                        structure, simply add the property and 
                        its value to the structure.

            scatter (optional)

                A structure with two components:

                    enabled

                        A logical value, indicating whether an ``scatter`` 
                        will be added to the plot or not. If set to ``False``
                        no target scatter circle will be added.

                    kws

                        A structure whose components will be passed to the 
                        matplotlib's ``scatter()`` function. For example:

                        .. code-block:: python

                            scatter.kws.s = 20
                            scatter.kws.color = "orangered"

                        NOTE: If a desired property is missing in the 
                        NOTE: structure, simply add the property and 
                        NOTE: its value to the structure.

            currentFig

                An structure which is initially ``None``, but upon making a plot, 
                is populated with attributes representing all outputs from 
                matplotlib function calls, with the following attributes:

                    axes

                        The output of matplotlib's ``gca()`` function or,
                        the input axes object by the user.

                    axhline

                        The output of matplotlib's ``axhline()`` function.

                    axvline

                        The output of matplotlib's ``axvline()`` function.

                    scatter

                        The output of matplotlib's ``scatter()`` function.

        **Returns**

            self

                An object of class ``Target``.

    ---------------------------------------------------------------------------
    """

    ################################################################################################################################
    #### __init__
    ################################################################################################################################

    def __init__( self
                , values : tp.Optional[ tp.Union[ tp.List[int] , tp.List[float] ] ] = [0,0]
                , axes = None #: tp.Optional[ plt.axes ] = None
                ):

        self.label = None
        self.value = values

        self.axvline = Struct()
        self.axvline.kws = Struct()
        self.axvline.enabled = True

        self.axhline = Struct()
        self.axhline.kws = Struct()
        self.axhline.enabled = True

        self.scatter = Struct()
        self.scatter.kws = Struct()
        self.scatter.enabled = True

        self.currentFig = Struct()
        self.currentFig.axes = axes
        if self.currentFig.axes is not None: self.currentFig.figure = self.currentFig.axes.get_figure()

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
                ``make()`` method of the target object.

            **Returns**

                Any return value from the ``make()``
                method of the target object.

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

        Add a target on an existing plot (the current active axes object)
        based on the ``values`` attribute of the target object.

            **Parameters**

                reself

                    A logical variable. If ``True``, an instance of 
                    the object will be returned upon exit to the calling 
                    routine. The default value is ``False``.

            **Returns**

                the object self if ``reself = True`` otherwise, ``None``.
                However, this method causes side-effects by manipulating
                the existing attributes of the object.

        """

        for key, val in kwargs.items():
            if hasattr(self,key):
                setattr(self, key, val)
            elif key=="axes":
                self.currentFig.axes = val
                self.currentFig.figure = self.currentFig.axes.get_figure()
            else:
                raise Exception ( "Unrecognized input '"+key+"' class attribute detected." + newline
                                + self._getDocString()
                                )

        # set what to plot

        if not isinstance(self.value, list): self.value = [self.value]

        ############################################################################################################################
        #### setup the vertical line
        ############################################################################################################################

        if self.axvline.enabled:
            if isinstance(self.axvline.kws,Struct):
                #if targetExists: self.axvline.kws.x = self.value[0]
                if "linewidth" not in vars(self.axvline.kws).keys(): self.axvline.kws.linewidth = 0.5
                if "linestyle" not in vars(self.axvline.kws).keys(): self.axvline.kws.linestyle = "-"
                if "zorder" not in vars(self.axvline.kws).keys(): self.axvline.kws.zorder = 1000
                if "color" not in vars(self.axvline.kws).keys(): self.axvline.kws.color = "orangered"
            else:
                raise Exception ( "The axvline.kws component of the current Target object must" + newline
                                + "be an object of class Struct(), essentially a structure with components" + newline
                                + "whose names are the input arguments to the axvline() function of the" + newline
                                + "matplotlib library." + newline
                                + self._getDocString()
                                )

        ############################################################################################################################
        #### setup the horizontal line
        ############################################################################################################################

        if self.axhline.enabled:
            if isinstance(self.axhline.kws,Struct):
                #if targetExists: self.axhline.kws.y = self.value[1]
                if "linewidth" not in vars(self.axhline.kws).keys(): self.axhline.kws.linewidth = 0.5
                if "linestyle" not in vars(self.axhline.kws).keys(): self.axhline.kws.linestyle = "-"
                if "color" not in vars(self.axhline.kws).keys(): self.axhline.kws.color = "orangered"
                if "zorder" not in vars(self.axhline.kws).keys(): self.axhline.kws.zorder = 1001
            else:
                raise Exception ( "The axhline.kws component of the current Target object must" + newline
                                + "be an object of class Struct(), essentially a structure with components" + newline
                                + "whose names are the input arguments to the axhline() function of the" + newline
                                + "matplotlib library." + newline
                                + self._getDocString()
                                )

        ############################################################################################################################
        #### setup the scatter line
        ############################################################################################################################

        if self.scatter.enabled:
            if isinstance(self.scatter.kws,Struct):
                #if targetExists: self.scatter.kws.x = self.value[0]
                #if targetExists: self.scatter.kws.y = self.value[1]
                if "s" not in vars(self.scatter.kws).keys(): self.scatter.kws.s = 20
                if "color" not in vars(self.scatter.kws).keys(): self.scatter.kws.color = "orangered"
                if "zorder" not in vars(self.scatter.kws).keys(): self.scatter.kws.zorder = 1002
            else:
                raise Exception ( "The scatter.kws component of the current Target object must" + newline
                                + "be an object of class Struct(), essentially a structure with components" + newline
                                + "whose names are the input arguments to the scatter() function of the" + newline
                                + "matplotlib library." + newline
                                + self._getDocString()
                                )

        ############################################################################################################################
        ############################################################################################################################
        if self._isdryrun: return
        ############################################################################################################################
        ############################################################################################################################

        from matplotlib import pyplot as plt

        # generate figure and axes if needed

        if self.currentFig.axes is None: self.currentFig.axes = plt.gca()
        self.currentFig.figure = self.currentFig.axes.get_figure()

        try:
            plt.sca(self.currentFig.axes)
        except:
            self.currentFig.axes = plt.gca()
            self.currentFig.figure = plt.gcf()
        #plt.figure(self.currentFig.figure.number)

        ############################################################################################################################
        #### add target
        ############################################################################################################################

        xlimCurrent = self.currentFig.axes.get_xlim()
        ylimCurrent = self.currentFig.axes.get_ylim()

        if self.axvline.enabled:

                #if "ymin" not in vars(self.axvline.kws).keys(): self.axvline.kws.ymin = ylimCurrent[0]
                #if "ymax" not in vars(self.axvline.kws).keys(): self.axvline.kws.ymax = ylimCurrent[1]
                self.currentFig.axvline = self.currentFig.axes.axvline  ( self.value[0]
                                                                        , **vars(self.axvline.kws)
                                                                        )

        if self.axhline.enabled:

                #if "xmin" not in vars(self.axhline.kws).keys(): self.axhline.kws.xmin = xlimCurrent[0]
                #if "xmax" not in vars(self.axhline.kws).keys(): self.axhline.kws.xmax = xlimCurrent[1]
                self.currentFig.axhline = self.currentFig.axes.axhline  ( self.value[1]
                                                                        , **vars(self.axhline.kws)
                                                                        )

        if self.scatter.enabled:

                self.currentFig.scatter = self.currentFig.axes.scatter  ( self.value[0]
                                                                        , self.value[1]
                                                                        , **vars(self.scatter.kws)
                                                                        )

        self.currentFig.axes.set_xlim(xlimCurrent)
        self.currentFig.axes.set_ylim(ylimCurrent)

        ############################################################################################################################

        if reself: return self

    ################################################################################################################################
    #### _getDocString
    ################################################################################################################################

    def _getDocString(self):
        docString   = newline \
                    + "Here is the help information on the Target class:" + newline \
                    + newline \
                    + self.__doc__
        return docString

    ################################################################################################################################
