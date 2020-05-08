#**********************************************************************************************************************************
#**********************************************************************************************************************************
#
#  ParaMonte: plain powerful parallel Monte Carlo library.
#
#  Copyright (C) 2012-present, The Computational Data Science Lab
#
#  This file is part of ParaMonte library.
#
#  ParaMonte is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, version 3 of the License.
#
#  ParaMonte is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
#
#**********************************************************************************************************************************
#**********************************************************************************************************************************

import os as _os
import sys as _sys
import numpy as _np
import typing as _tp
import pandas as _pd
import platform as _platform

fileAbsDir = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.append(fileAbsDir)

from _message import note, warn, abort
import _visualization as vis
import _statistics as stats
import _dfutils as dfutils
import _pmutils as pmutils

####################################################################################################################################

class _Struct:
    pass

####################################################################################################################################

arch = "x86" if "32" in _platform.architecture()[0] else "x64"

####################################################################################################################################

names = _Struct()
names.paramonte = "ParaMonte"
names.paradram = "ParaDRAM"
names.paranest = "ParaNest"
names.paratemp = "ParaTemp"

####################################################################################################################################

# get version

class Version:

    def __init__(self):
        self._savedStruct = _Struct()
        from collections import OrderedDict
        self._versionFile = OrderedDict()
        self._versionFile["interface"] = ".VERSION"
        self._versionFile["kernel"] = ".VERSION_KERNEL"
        for key in self._versionFile.keys():
            setattr(self._savedStruct,key,None)
            setattr(self,key,self.get(key))

    def get(self, which : _tp.Optional[str] = None):
        this = self._checkVersionType("get",which)
        return "ParaMonte Python " + this.capitalize() + " Version " + self.dump(this)

    def dump(self, which : _tp.Optional[str] = None):
        this = self._checkVersionType("dump",which)
        versionFilePath = _os.path.join( fileAbsDir, self._versionFile[this] )
        with open(versionFilePath) as versionFile:
            version = versionFile.readline()
            version = "".join(version.splitlines())
        setattr(self._savedStruct,this,version)
        return version

    def _checkVersionType(self,versionClassMethodName,which):
        keys = list(self._versionFile.keys())
        if which is None or len(which)==0: return keys[0]
        if isinstance(which, str):
            for key in keys:
                if which.lower()==key: return key
        abort   ( msg   = "The " + versionClassMethodName + "() method of the Version() class only takes one input string argument with two possible values:\n\n"
                        + "    " + versionClassMethodName + "(\"" + keys[0] + "\")\n"
                        + "    " + versionClassMethodName + "(\"" + keys[1] + "\")"
                , methodName = names.paradram
                , marginTop = 1
                , marginBot = 1
                )

version = Version()

####################################################################################################################################
