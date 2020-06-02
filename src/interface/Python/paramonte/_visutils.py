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
#from matplotlib import pyplot as _plt

####################################################################################################################################
#### colorbar
####################################################################################################################################

def isColorBar(ax):
    """
    Guesses whether a set of Axes is home to a colorbar

    Parameters
    ----------
        ax
            Axes instance

    Returns
    -------
        True if the x xor y axis satisfies all of the 
        following and thus looks like a colorbar:
            No ticks, no tick labels, no axis label

    """

    xcb = (len(ax.get_xticks()) == 0) and (len(ax.get_xticklabels()) == 0) and (len(ax.get_xlabel()) == 0) #and (ax.get_xlim() == (0, 1))
    ycb = (len(ax.get_yticks()) == 0) and (len(ax.get_yticklabels()) == 0) and (len(ax.get_ylabel()) == 0) #and (ax.get_ylim() == (0, 1))
    return xcb != ycb

def isPresentColorBar():
    from matplotlib import pyplot as _plt
    return any([isColorBar(ax) for ax in _np.atleast_1d(_plt.gcf().axes).flatten()])

####################################################################################################################################
#### ParaMonteFigure class
####################################################################################################################################

class ParaMonteFigure:
    pass

####################################################################################################################################
#### Target class
####################################################################################################################################

class Target:
    """
    This is the Target class for generating instances 
    of a target on the current active figure and axis.

    Usage
    -----
    first generate an object of this class by optionally 
    passing the following parameters described below. Then call 
    the plot() method. The generated object is also callable with 
    the same input parameters as the object's constructor.

    Parameters
    ----------
        value
            a pair (list, array, or tuple) of floats, representing the (x,y) cordinates of the target.
        axes
            the axes object to which the target must be added.
            The default is None in which case the output of matplotlib's gca()
            function wil be used to get the current active axes.
        axhline_kws
            optional dictionary of keyword arguments to be passed to matplotlib's 
            axhline() function. For example: 
                axhline = { "linewidth": 0.5,
                            "color": "orangered",
                            "linestyle": ":",
                            "xmin": 0.1
                            "xmax": 1.9
                            }
            The default is {}, which will be appropriately populated.
            If set to None, no horizontal target line will be added.
        axvline_kws
            optional dictionary of keyword arguments to be passed to matplotlib's 
            axvline() function. For example: 
                axvline = { "linewidth": 0.5,
                            "color": "orangered",
                            "linestyle": ":",
                            "ymin": 0.1
                            "ymax": 1.9
                            }
            The default is {}, which will be appropriately populated.
            If set to None, no vertical target line will be added.
        scatter_kws
            optional dictionary of keyword arguments to be passed to matplotlib's 
            scatter() function. For example: 
                scatter_kws = {"s":10,"color":"green"}
            The default is {}.
            If set to None, no target point will be plotted.
        outputFile
            optional string representing the name of the output file in which 
            the figure will be saved. If not provided, no output file will be generated.

    Attributes
    ----------
        All of the parameters described above, except axes.
            The input axes object (whether provided by the user or fetched by the class)
            will be stored as a component of the object's attribute currentFig.
        currentFig
            an object of class ParaMonteFigure which is initially None, but upon 
            making a plot, is populated with attributes representing all outputs 
            from matplotlib/seaborn function calls, with the following attributes:
                axes
                    the output of matplotlib's gca() function or,
                    the input axes object by the user.
                axhline
                    the output of matplotlib's axhline() function.
                axvline
                    the output of matplotlib's axvline() function.
                scatter
                    the output of matplotlib's scatter() function.

    Returns
    -------
        Target

    ----------------------------------------------------------------------
    """

    def __init__( self
                , value         : _tp.Optional[ _tp.Union[ _tp.List[int] , _tp.List[float] ] ] = None
                , axes          = None #: _tp.Optional[ _plt.axes ] = None
                , axhline_kws   : _tp.Optional[_tp.Dict] = ()
                , axvline_kws   : _tp.Optional[_tp.Dict] = ()
                , scatter_kws   : _tp.Optional[_tp.Dict] = ()
                , outputFile    : _tp.Optional[str] = None
                ):

        self.value          = value
        self.axhline_kws    = axhline_kws
        self.axvline_kws    = axvline_kws
        self.scatter_kws    = scatter_kws
        self.outputFile     = outputFile
        if axes is None:
            self.currentFig = None
        else:
            self.currentFig = ParaMonteFigure()
            self.currentFig.axes = axes
            self.currentFig.figure = axes.get_figure()

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
            elif key=="axes": 
                self.currentFig = ParaMonteFigure()
                self.currentFig.axes = kwargs[key]
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
        Add a target on an existing plot (the current active axes object)
        based on the 'value' attribute of the target object.

        Parameters
        ----------
            None

        Returns
        -------
            None. However, this method causes side-effects by manipulating 
            the existing attributes of the object.

        """

        # set what to plot

        targetExists = self.value is not None
        axhlineEnabled = self.axhline_kws is not None
        axvlineEnabled = self.axvline_kws is not None
        scatterEnabled = self.scatter_kws is not None

        # setup the horizontal line

        if self.axhline_kws==(): self.axhline_kws={}
        if isinstance(self.axhline_kws,dict):
            if targetExists: self.axhline_kws["y"] = self.value[1]
            if "linewidth" not in self.axhline_kws.keys():  self.axhline_kws["linewidth"] = 0.5
            if "linestyle" not in self.axhline_kws.keys():  self.axhline_kws["linestyle"] = "-"
            if "color" not in self.axhline_kws.keys():      self.axhline_kws["color"] = "orangered"
        elif axhlineEnabled:
            raise Exception ( "The input argument 'axhline_kws' must be None or a dictionary,\n"
                            + "each (key,value) pair of which represents an (attribute,value)\n"
                            + "of matplotlib library's axhline() function."
                            )

        # setup the vertical line

        if self.axvline_kws==(): self.axvline_kws={}
        if isinstance(self.axvline_kws,dict):
            if targetExists: self.axvline_kws["x"] = self.value[0]
            if "linewidth" not in self.axvline_kws.keys():  self.axvline_kws["linewidth"] = 0.5
            if "linestyle" not in self.axvline_kws.keys():  self.axvline_kws["linestyle"] = "-"
            if "color" not in self.axvline_kws.keys():      self.axvline_kws["color"] = "orangered"
        elif axvlineEnabled:
            raise Exception ( "The input argument 'axvline_kws' must be None or a dictionary,\n"
                            + "each (key,value) pair of which represents an (attribute,value)\n"
                            + "of matplotlib library's axvline() function."
                            )

        # setup the scatter plot

        if self.scatter_kws==(): self.scatter_kws={}
        if isinstance(self.scatter_kws,dict):
            if targetExists: self.scatter_kws["x"] = self.value[0]
            if targetExists: self.scatter_kws["y"] = self.value[1]
            if "s" not in self.scatter_kws.keys():          self.scatter_kws["s"] = 20
            if "color" not in self.scatter_kws.keys():      self.scatter_kws["color"] = "orangered"
        elif scatterEnabled:
            raise Exception ( "The input argument 'scatter_kws' must be None or a dictionary,\n"
                            + "each (key,value) pair of which represents an (attribute,value)\n"
                            + "of matplotlib library's scatter() function."
                            )

        #########################
        if self._isdryrun: return
        #########################

        from matplotlib import pyplot as _plt

        # generate figure and axes if needed

        if self.currentFig is None:
            self.currentFig = ParaMonteFigure()
            self.currentFig.axes = _plt.gca()
        else:
            _plt.figure(self.currentFig.figure.number)
            _plt.sca(self.currentFig.axes)

        #########################

        # add target

        xlimCurrent = self.currentFig.axes.get_xlim()
        ylimCurrent = self.currentFig.axes.get_ylim()

        legendEnabled = False

        if axhlineEnabled:

                if "xmin" not in self.axhline_kws.keys(): self.axhline_kws["xmin"] = xlimCurrent[0]
                if "xmax" not in self.axhline_kws.keys(): self.axhline_kws["xmax"] = xlimCurrent[1]
                if "zorder" not in self.axhline_kws.keys(): self.axhline_kws["zorder"] = 21
                self.currentFig.axhline = self.currentFig.axes.axhline ( **self.axhline_kws )
                legendEnabled = "label" in self.axhline_kws.keys()

        if axvlineEnabled:

                if "xmin" not in self.axvline_kws.keys(): self.axvline_kws["ymin"] = ylimCurrent[0]
                if "xmax" not in self.axvline_kws.keys(): self.axvline_kws["ymax"] = ylimCurrent[1]
                if "zorder" not in self.axvline_kws.keys(): self.axvline_kws["zorder"] = 21
                self.currentFig.axvline = self.currentFig.axes.axvline ( **self.axvline_kws )
                legendEnabled = "label" in self.axvline_kws.keys()

        if scatterEnabled:

                if "zorder" not in self.scatter_kws.keys(): self.scatter_kws["zorder"] = 21
                self.currentFig.scatter = self.currentFig.axes.scatter( **self.scatter_kws )
                legendEnabled = "label" in self.scatter_kws.keys()

        self.currentFig.axes.set_xlim(xlimCurrent)
        self.currentFig.axes.set_ylim(ylimCurrent)

        # add legend if requested

        if legendEnabled: self.currentFig.axes.legend()

        if self.outputFile is not None:
            self.currentFig.figure.savefig  ( self.outputFile
                                            , bbox_inches = 'tight'
                                            , pad_inches = 0.0
                                            )
