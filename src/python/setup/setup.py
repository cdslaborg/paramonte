####################################################################################################################################
####################################################################################################################################
####                                                                                                                            ####
####    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           ####
####                                                                                                                            ####
####    Copyright (C) 2012-present, The Computational Data Science Lab                                                          ####
####                                                                                                                            ####
####    This file is part of the ParaMonte library.                                                                             ####
####                                                                                                                            ####
####    LICENSE                                                                                                                 ####
####                                                                                                                            ####
####       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          ####
####                                                                                                                            ####
####################################################################################################################################
####################################################################################################################################

import setuptools

with open("README.md", "r") as fh: readmeFileContents = fh.read()

with open("paramonte/auxil/VERSION.md") as versionFile: version = versionFile.readline()

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
