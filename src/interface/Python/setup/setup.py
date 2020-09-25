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

import setuptools

with open("README.md", "r") as fh: readmeFileContents = fh.read()

with open("paramonte/auxil/.VERSION_INTERFACE") as versionFile: version = versionFile.readline()

setuptools.setup( name                          = "paramonte"
                , version                       = version
                , author                        = "Amir Shahmoradi, Fatemeh Bagheri, Joshua Alexander Osborne"
                , author_email                  = "shahmoradi@utexas.edu"
                , description                   = "Plain Powerful Parallel Monte Carlo and adaptive MCMC Library"
                , long_description              = readmeFileContents
                , long_description_content_type = "text/markdown"
                , url                           = "https://github.com/cdslaborg/paramonte"
                , packages                      = setuptools.find_packages()
                , python_requires               = '>=3.0'
                , license                       = "License :: OSI Approved :: MIT License"
                , include_package_data          = True
                , install_requires              =   [ "numpy"
                                                    , "pandas"
                                                    #, "scipy"
                                                    #, "seaborn"
                                                    #, "matplotlib"
                                                    #, OTHERS:
                                                    #, "platform"
                                                    #, "weakref"
                                                    #, "typing"
                                                    #, "time"
                                                    #, "sys"
                                                    #, "os"
                                                    ]
                , classifiers                   =   [ "Development Status :: 5 - Production/Stable"
                                                    , "Natural Language :: English"
                                                    , "Programming Language :: Python :: 3"
                                                    , "Programming Language :: Fortran"
                                                    , "Programming Language :: C"
                                                    , "License :: OSI Approved :: MIT License"
                                                    , "Operating System :: OS Independent"
                                                    , "Topic :: Scientific/Engineering"
                                                    , "Topic :: Scientific/Engineering :: Physics"
                                                    , "Topic :: Scientific/Engineering :: Astronomy"
                                                    , "Topic :: Scientific/Engineering :: Mathematics"
                                                    , "Topic :: Scientific/Engineering :: Visualization"
                                                    , "Topic :: Scientific/Engineering :: Bio-Informatics"
                                                    , "Topic :: Scientific/Engineering :: Atmospheric Science"
                                                    , "Operating System :: OS Independent"
                                                    , "Operating System :: Microsoft :: Windows"
                                                    , "Operating System :: POSIX :: Linux"
                                                    , "Operating System :: MacOS"
                                                    , "Operating System :: Unix"
                                                    ]
                )
