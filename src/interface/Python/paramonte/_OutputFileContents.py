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

import numpy as np
import _paramonte as pm

class OutputFileContents:
    """

    This is the **OutputFileContents** base class for the ParaMonte 
    sampler output file contents classes. **This class is NOT meant to be 
    directly accessed or called by the user of the ParaMonte library.**
    However, its children, such the ParaDRAM sampler class will be 
    indirectly accessible to the public. 

        **Parameters**

            file

                The full path to the file.

            methodName

                A string representing the name of the ParaMonte sampler used
                to call the constructor of the ``OutputFileContents`` class.

            reportEnabled

                A logical input parameter indicating whether the ParaMonte
                automatic guidelines to the standard output should be provided 
                or not. The default value is ``True``.

        **Methods**

            See below for information on the methods.  

        **Returns**

            Object of class OutputFileContents

    """

    def __init__( self
                , file          : str
                , methodName    : str
                , reportEnabled : bool
                ):

        self.file = file
        self._methodName = methodName
        self._reportEnabled = reportEnabled
        self._progress = pm.utils.Progress  ( msg = "reading the file contents... "
                                            , methodName = methodName
                                            , reportEnabled = reportEnabled
                                            , end = ""
                                            )

    ################################################################################################################################
