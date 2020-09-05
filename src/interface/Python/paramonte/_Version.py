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

import os

####################################################################################################################################
#### Version
####################################################################################################################################

class Version:
    """

    This is the Version class for generating objects
    that contain the methods for getting and dumping
    the python-interface or kernel versions of the
    ParaMonte library installation on the system.

        **Parameters**

            versionPath

                A string containing the path to either the
                ParaMonte kernel or interface version file.

            versionType

                A string containing the type of the version
                file. It can be one of the following values:

                    "interface"

                        implying the Python-interface version
                        number of the ParaMonte library.

                    "kernel"

                        implying the kernel-routines version
                        number of the ParaMonte library.

    """

    ################################################################################################################################
    #### __init__
    ################################################################################################################################

    def __init__(self,versionPath,versionType):
        self._versionList = ["interface","kernel"]
        self._versionPath = versionPath
        self._versionType = versionType
        self._versionSave = None
        self._checkVersionType()

    ################################################################################################################################
    #### get
    ################################################################################################################################

    def get(self):
        """

        Get the Python-interface or kernel version of the
        ParaMonte library, in verbose format.

            **Parameters**

                None

            **Returns**

                None

        """
        return "ParaMonte Python " + self._versionType.capitalize() + " Version " + self.dump()

    ################################################################################################################################
    #### dump
    ################################################################################################################################

    def dump(self):
        """

        Dump **only the version number** of either
        the Python-interface or kernel of the
        ParaMonte library.

            **Parameters**

                None

            **Returns**

                None

        """
        for versionType in self._versionList:
            if versionType==self._versionType:
                if self._versionSave is None:
                    versionFileName = ".VERSION_" + versionType.upper()
                    versionFilePath = os.path.join(self._versionPath, versionFileName)
                    try:
                        with open(versionFilePath,"r") as versionFile:
                            self._versionSave = versionFile.readline().strip("\n")
                    except:
                        self._versionSave = "UNKNOWN"
                    return self._versionSave
                else:
                    return self._versionSave

    ################################################################################################################################
    #### _checkVersionType
    ################################################################################################################################

    def _checkVersionType(self):
        versionTypeNotFound = True
        for versionType in self._versionList:
            if versionType==self._versionType:
                versionTypeNotFound = False
                break
        if versionTypeNotFound:
            raise Exception ( "The input versionType is not a valid recognized version type. Possible values:\n"
                            + "\n".join(self._versionList)
                            )
