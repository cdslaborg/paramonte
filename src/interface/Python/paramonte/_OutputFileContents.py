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

import _paramonte as pm

####################################################################################################################################
#### OutputFileContents
####################################################################################################################################

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

####################################################################################################################################
