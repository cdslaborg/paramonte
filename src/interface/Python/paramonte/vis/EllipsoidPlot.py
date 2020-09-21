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

import numpy as np
import typing as tp
import weakref as wref
import _paramonte as pm
from paramonte.vis.Target import Target
from paramonte.vis._BasePlot import BasePlot

Struct = pm.Struct
newline = pm.newline

####################################################################################################################################
#### EllipsoidPlot class
####################################################################################################################################

class EllipsoidPlot(BasePlot):
    """

    This is the EllipsoidPlot class for generating instances of
    ellipsoid-evolution plots via matplotlib's builtin function 
    `plot()`. It generates a plot of 2D ellipsoids corresponding 
    to the input list of covariance/correlation matrices.

    **NOTE**

    This is a low-level ParaMonte class and is NOT meant
    to be directly instantiated by the users,
    although it is possible.

        **Parameters**

            matrix

                A 3D numpy array representing the covariance/correlation 
                matrices of the ellipsoids to plot. The size of the array 
                is ``(count, ndim, ndim)`` where ``count`` represents the 
                number of ellipsoids and ``ndim`` is in the number of 
                dimensions of the (hyper-)ellipsoids.

            plotType

                A string indicating the type of the plot to be made. 
                Possible values include: 

                    ``"covmat2"``
                    ``"covmat3"``
                    ``"cormat2"``
                    ``"cormat3"``

            methodName (optional)

                A string indicating the ParaMonte sampler that 
                is instantiating this EllipsoidPlot object. This 
                variable is solely used for progress reporting 
                and has no relevance to the plotting.
                The default value is ``"ParaMonte"``.

            reportEnabled (optional)

                A logical input parameter indicating whether the 
                ParaMonte automatic guidelines to the standard output 
                should be provided or not. The default value is ``True``.

            resetPlot (optional)

                An user-specified function that can be later called to 
                reset the settings of the instantiated Ellipsoid plot.
                If ``None`` is provided, then the `reset()`` of the 
                Ellipsoid object will be set to the default ``_reset()``
                function of the Ellipsoid object.
                The default value is ``None``.

        **Attributes**

            ndim

                An integer, representing the number of dimensions of the 
                domain of the objective function as inferred from the 
                contents of the dataFrame. It is the rank of the 
                input matrix.

            matrix

                The same as the corresponding input parameter ``matrix``
                to the ``EllipsoidPlot`` class constructor in the above.

            center

                A 2D numpy array representing the centers of the ellipsoids
                to plot. The size of the array is ``(count, ndim)`` where
                ``count`` represents the number of ellipsoids and ``ndim``
                is in the number of dimensions of the (hyper-)ellipsoids.

            zdata (optional, exists only in the case of 3D Ellipsoid plots)

                Optional property that, if provided, results in the creation 
                of 3D ellipsoid plots and its data represents the z-axis 
                values of the plot. It must be a numeric vector of the 
                same length as ``center`` and ``matrix``.

                The default value is the index of the covariance/correlation 
                matrices in the input data frame.

            dimensionPair

                A pair of indices (vector of length 2) whose value determine 
                the rows and columns from the covariance/correlation matrix 
                which will be plotted. 
                The default value is ``dimensionPair = [1,2]`` 
                (corresponding to Python indices of [0,1]).

                Example usage:

                    1.  ``dimensionPair = [1,2]``
                    2.  ``dimensionPair = [3,1]``

                **WARNING**

                In all cases, the indices must be distinct from 
                each other and less than ``ndim`` where ``ndim`` 
                is the rank of the covariance/correlation matrices.

            npoint

                The number of points used to represent the ellipsoids. 
                The higher the value of ``npoint``, the higher-resolution 
                the resulting ellipsoids will look like.

                The default value is ``100``.

            colorbar (exists only for plots that require colorbar)

                A structure with two attributes:

                    enabled

                        A boolean indicating whether a call to the ``colorbar()``
                        function of the matplotlib library should be made or not.
                        If a call is made, a new figure will be generated.
                        Otherwise, the current active figure will be used.

                    kws

                        A structure whose components are directly passed as 
                        keyword arguments to the ``colorbar()`` function of 
                        the matplotlib library.

                        **NOTE**

                        If a desired property is missing among the ``kws`` 
                        attributes, simply add the field and its value to 
                        the component.

                A colorbar will be added to a plot only if a color-mappings 
                is requested in the plot.

            colormap

                A structure with the following components:
                
                    enabled

                        A logical value indicating whether the ellipsoids 
                        must be color-mapped or not.

                    name

                        A string representing the name of the color mapping 
                        from the matplotlib library's ``cm`` class, to be 
                        applied to the plot. 

                    values

                        A numeric vector of the same length as the number 
                        of ellipsoids in ``matrix``, each element of which 
                        will determine the color of the corresponding 
                        ellipsoid in the plot.

                        If the values are set to ``None`` or an empty list, 
                        then a range of values corresponding to the count of 
                        ellipsoids will be used.

                    _rgba

                        The RGB vector corresponding to the values. This is 
                        not meant to be directly manipulated by the user.

                Example usage:

                    1.  ``colormap.enabled = True``
                    2.  ``colormap.name = "winter"``

            legend (may not exist for some types of plots)

                A structure with two attributes:

                    enabled

                        A boolean indicating whether a call to the ``legend()``
                        function of the matplotlib library should be made or not.
                        If a call is made, a new figure will be generated.
                        Otherwise, the current active figure will be used.

                    kws

                        A structure whose components are directly passed as 
                        keyword arguments to the ``legend()`` function.

                        Example usage:

                        .. code-block:: python

                            legend.kws.labels = ["Variable1", "Variable2"]

                        **NOTE**

                        If a desired property is missing among the ``kws`` 
                        attributes, simply add the field and its value to 
                        the component.

                A legend will be added to a plot only if no color-mappings 
                are requested in the plot.

            plot

                A structure with two components:

                    enabled

                        A logical value. If ``true``, 
                        the ellipsoids will be plotted.

                    kws
                        A structure whose components will be 
                        passed to the ``plot()`` function of 
                        the matplotlib library.

                        Example usage:

                            ``plot.kws.linewidth = 2``

                        **NOTE**

                        If a desired property is missing among the ``kws`` 
                        attributes, simply add the field and its value to 
                        the component.

            title

                A structure with the following components:

                    enabled

                        A logical value. If ``true``, 
                        a title will be added to the plot.

                    label

                        A string representing the content of 
                        the title to be added to the plot.

                    kws

                        A structure whose components will be 
                        passed to the ``title()`` function 
                        of the matplotlib library.

                        **NOTE**

                        If a desired property is missing among the ``kws`` 
                        attributes, simply add the field and its value to 
                        the component.

            target

                An object of class ``Target`` for adding target 
                values to the plots. For more information, see the 
                documentation of the ``Target()`` class.

        **Superclass Attributes**

            See the documentation for the ``BasePlot`` class.

        **Returns**

            An object of ``EllipsoidPlot`` class.

    ---------------------------------------------------------------------------
    """

    ################################################################################################################################
    #### __init__
    ################################################################################################################################

    def __init__( self
                , matrix        : tp.Union[ np.ndarray, tp.List[np.ndarray] ]
                , plotType      : str
                , methodName    : tp.Optional[ str ] = "ParaMonte"
                , reportEnabled : tp.Optional[ bool ] = True
                , resetPlot     = None
                ):

        self.matrix = matrix
        if isinstance(matrix, np.ndarray):
            self.numEllipsoid = np.shape(matrix)[0]
        #elif instance(matrix, list):
        #    self.numEllipsoid = len(matrix)
        else:
            pm.abort( msg   = "The input argument ``matrix`` must be a numpy 3D" + newline
                            + "array of matrices representing the covariance/correlation " + newline
                            + "matrices of the ellipsoids or. The first dimension of array" + newline
                            + "must count over the ellipsoids." + newline
                    , marginTop = 1
                    , marginBot = 1
                    , methodName = methodName
                    )
        import pandas as pd
        self._dfIndex = pd.DataFrame.from_dict({ "indices" : list(range( self.numEllipsoid )) })

        super().__init__( plotType = plotType
                        , dataFrame = self._dfIndex
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

        self.npoint = None
        self.center = None
        self.dimensionPair = None
        if self._type.is3d: self.zdata = None

        super()._reset()
        if not self._type.is3d: self.target = Target()

        self.plot = Struct()
        self.plot.enabled = True
        self.plot.kws = Struct()

        self.colormap = Struct()
        self.colormap.enabled = True
        self.colormap.values = None
        self.colormap.name = None

        self.colorbar = Struct()
        self.colorbar.enabled = True
        self.colorbar.kws = Struct()

        self.title = Struct()
        self.title.enabled = False
        self.title.label = ""
        self.title.kws = Struct()
        self.title.kws.loc = "center"

        self._isdryrun = True
        self.make()
        self._isdryrun = False

        self.legend.enabled = False
        self.legend._labels = None

    ################################################################################################################################
    #### __call__
    ################################################################################################################################

    def __call__( self
                , reself : tp.Optional[ bool ] = False
                , **kwargs
                ):
        """

        Call the ``make()`` method of the current 
        instance of the class.

            **Parameters**

                Any arguments that can be passed to the 
                ``make()`` method of the plot object.

            **Returns**

                Any return value from the ``make()`` 
                method of the plot object.

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

        Generate an EllipsoidPlot from the input data the 
        EllipsoidPlot object.

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
                raise Exception ( "Unrecognized input '"+key+"' class attribute detected." + newline
                                + self._getDocString()
                                )

        if self.dimensionPair is None or (isinstance(self.dimensionPair,list) and len(self.dimensionPair)!=2): self.dimensionPair = [0,1]
        self.dimensionPair = np.int32(self.dimensionPair)
        if not isinstance(self.npoint, int): self.npoint = 100
        matrixShape = np.shape(self.matrix)
        self.numEllipsoid = matrixShape[0]
        self.ndim = matrixShape[1]

        # set what to plot

        ############################################################################################################################
        #### line plot properties
        ############################################################################################################################

        if isinstance(self.plot.kws,Struct):

            if "linewidth" in vars(self.plot.kws).keys():
                if self.plot.kws.linewidth==0: self.plot.kws.linewidth = 1
            else:
                self.plot.kws.linewidth = 1

            if "zorder" not in vars(self.plot.kws).keys(): self.plot.kws.zorder = 1

        else:

            raise Exception ( "The plot.kws component of the current EllipsoidPlot object must" + newline
                            + "be an object of class Struct(), essentially a structure with components" + newline
                            + "whose names are the input arguments to the plot() function of the" + newline
                            + "matplotlib library."
                            )

        ############################################################################################################################
        #### legend properties
        ############################################################################################################################

        if self.legend.enabled:
            if not isinstance(self.legend.kws,Struct):
                raise Exception ( "The legend.kws component of the current EllipsoidPlot object must" + newline
                                + "be an object of class Struct(), essentially a structure with components" + newline
                                + "whose names are the input arguments to the legend() function of the" + newline
                                + "matplotlib library."
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
                raise Exception ( "The figure.kws component of the current EllipsoidPlot object must" + newline
                                + "be an object of class Struct(), essentially a structure with components" + newline
                                + "whose names are the input arguments to the figure() function of the" + newline
                                + "matplotlib library."
                                )

        ############################################################################################################################
        ############################################################################################################################
        if self._isdryrun: return
        ############################################################################################################################
        ############################################################################################################################

        from matplotlib import pyplot as plt

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

        ############################################################################################################################
        #### check the consistency of zdata
        ############################################################################################################################

        if self._type.is3d:
            if self.zdata is None:
                self.zdata = self._dfIndex.values
                self._zlabel = "Count"
            elif self.numEllipsoid != len(self.zdata):
                raise Exception ( "The lengths of the first dimension of the matrix and the zdata arrays must agree." + newline
                                + "This length is representative of the number of ellipsoids to plot."
                                )
            else:
                self._zlabel = None

        ############################################################################################################################
        #### set colormap data
        ############################################################################################################################

        if self.colormap.enabled:

            if self.colormap.values is None:
                self.colormap.values = self._dfIndex.values
                self.colorbar.label = "Count"
            else:
                self.colorbar.label = None

            self.colormap.values = np.double( self.colormap.values ) # ensure vectorized-indexing is possible

            if self.colormap.name is None or not isinstance(self.colormap.name, str):
                self.colormap.name = "winter"

            import matplotlib.colors as colors
            import matplotlib.cm as cmx

            colormapObject = plt.get_cmap(self.colormap.name)
            cNorm = colors.Normalize( vmin = np.min(self.colormap.values), vmax = np.max(self.colormap.values), clip = True )
            self.colorbar.kws.mappable = cmx.ScalarMappable( norm = cNorm, cmap = colormapObject )
            #self.colormap._rgba = colormapObject( self.colormap.values )
            self.colormap._rgba = np.zeros( ( len(self.rows) , 4) )
            for i, row in enumerate(self.rows):
                self.colormap._rgba[i,:] = self.colorbar.kws.mappable.to_rgba(self.colormap.values[row])

        ############################################################################################################################
        #### check the consistency of the lengths
        ############################################################################################################################

        if self.center is None: self.center = np.zeros( (self.numEllipsoid, self.ndim) )

        if self.numEllipsoid != np.shape(self.center)[0]: 
            raise Exception ( "The lengths of the first dimensions of the matrix and center vectors must agree." + newline
                            + "This length is representative of the number of ellipsoids to plot."
                            )

        ############################################################################################################################
        #### make ellipsoids
        ############################################################################################################################

        if self.plot.enabled: self.currentFig.plotList = []
        for i, row in enumerate(self.rows):

            if self._type.is3d: zdata = np.ones(self.npoint) * self.zdata[row]

            if self.colormap.enabled: self.plot.kws.color = self.colormap._rgba[i]

            fulMatrix = self.matrix[row,:,:]
            subMatrix = fulMatrix[:,self.dimensionPair][self.dimensionPair]

            # get ellipsoid boundary

            bcrd = self.getEllipsoidBoundary( covMat = subMatrix
                                            , meanVec = self.center[row,:]
                                            , npoint = self.npoint
                                            )

            if self.plot.enabled:

                if self._type.is3d:

                    self.currentFig.plotList.append( self.currentFig.axes.plot  ( bcrd[0,:]
                                                                                , bcrd[1,:]
                                                                                , zdata
                                                                                , **vars(self.plot.kws)
                                                                                ) )

                else:

                    self.currentFig.plotList.append( self.currentFig.axes.plot  ( bcrd[0,:]
                                                                                , bcrd[1,:]
                                                                                , **vars(self.plot.kws)
                                                                                ) )

        ############################################################################################################################
        #### add axes and title labels
        ############################################################################################################################

        # add axis labels

        self.currentFig.axes.set_xlabel("Dimension " + str(self.dimensionPair[0]+1))
        self.currentFig.axes.set_ylabel("Dimension " + str(self.dimensionPair[1]+1))
        if self._type.is3d and isinstance(self._zlabel, str): self.currentFig.axes.set_zlabel(self._zlabel)

        # set title properties

        if self.title.enabled:
            if not isinstance(self.title.label, str):
                raise Exception ("The title component of an EllipsoidPlot object must be a string or character vector.")
            else:
                plt.title( self.title.label, **vars(self.title.kws) )

        ############################################################################################################################
        #### add colorbar
        ############################################################################################################################

        if self.colormap.enabled and self.colorbar.enabled:
            #self.colorbar.kws.cax = self.currentFig.axes
            self.currentFig.colorbar = self.currentFig.figure.colorbar( **vars(self.colorbar.kws) )
            if isinstance(self.colorbar.label, str): self.currentFig.colorbar.set_label( label = self.colorbar.label )

        ############################################################################################################################
        #### finalize plot
        ############################################################################################################################

        self._finalizeBasePlot()

        ############################################################################################################################

        if reself: return self

    ################################################################################################################################
    #### getEllipsoidBoundary
    ################################################################################################################################

    def getEllipsoidBoundary( self
                            , covMat
                            , meanVec
                            , npoint = 50
                            ):
        """
        
        Return the coordinates of the boundary of an 
        ellipsoid represented by the input ``covMat``
        covariance matrix and the ``meanVec`` center.

            **Parameters**

                covMat

                    The representative covariance matrix 
                    of the ellipsoid.

                meanVec

                    A vector representing the center 
                    of the ellipsoid.

                npoint

                    The number of points with which the 
                    ellipsoid is represented. The more,
                    the higher the resolution will be.
                    The default value is 50.

            **Returns**

                A matrix of size ``(2, npoint)`` that contains 
                the boundary points of the generated ellipsoid 
                corresponding to the input matrix. 

        """

        from scipy import linalg as la

        independentVariable = np.linspace(0, 2*np.pi, npoint, endpoint = True)
        ap = np.zeros( (2, npoint) )
        ap[0,:] = np.cos(independentVariable)
        ap[1,:] = np.sin(independentVariable)
        eigenValues, eigenVecMat = la.eig(covMat)
        eigenValMat = np.eye(2)
        eigenValMat[0,0] = np.sqrt(eigenValues[0].real)
        eigenValMat[1,1] = np.sqrt(eigenValues[1].real)
        meanVecRepMat = np.ones( (2,npoint) )
        meanVecRepMat[0,:] = meanVec[0]
        meanVecRepMat[1,:] = meanVec[1]
        return np.matmul( np.matmul(eigenVecMat.real , eigenValMat) , ap ) + meanVecRepMat 

    ################################################################################################################################
    #### _getDocString
    ################################################################################################################################

    def _getDocString(self):
        docString   = newline \
                    + "Here is the help information on the Ellipsoid class:" + newline \
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
