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

import setuptools

with open("README.md", "r") as fh:
    readmeFileContents = fh.read()

with open("paramonte/.VERSION") as versionFile:
    version = versionFile.readline()
versionFile.close()

setuptools.setup( name                          = "paramonte"
                , version                       = version
                , author                        = "Amir Shahmoradi, Fatemeh Bagheri"
                , author_email                  = "shahmoradi@utexas.edu, Fatemeh.Bagheri@uta.edu"
                , description                   = "Plain Powerful Parallel Monte Carlo Library"
                , long_description              = readmeFileContents
                , long_description_content_type = "text/markdown"
                , url                           = "https://github.com/cdslaborg/paramonte"
                , packages                      = setuptools.find_packages()
                , python_requires               = '>=3.0'
                , license                       = "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)"
                , include_package_data          = True
                , install_requires              =   [ "numpy"
                                                    , "scipy"
                                                    , "pandas"
                                                    , "seaborn"
                                                    , "matplotlib"
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
                                                    , "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)"
                                                    , "Operating System :: OS Independent"
                                                    , "Topic :: Scientific/Engineering"
                                                    , "Topic :: Scientific/Engineering :: Physics"
                                                    , "Topic :: Scientific/Engineering :: Mathematics"
                                                    , "Topic :: Scientific/Engineering :: Visualization"
                                                    , "Operating System :: OS Independent"
                                                    , "Operating System :: Microsoft :: Windows"
                                                    , "Operating System :: POSIX :: Linux"
                                                    , "Operating System :: MacOS"
                                                    , "Operating System :: Unix"
                                                    ]
                )
