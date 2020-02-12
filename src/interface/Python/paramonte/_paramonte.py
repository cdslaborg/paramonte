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
_sys.path.append(_os.path.dirname(__file__))

from _message import note, warn, abort
import _visualization as vis
import _statistics as stats
import _dfutils as dfutils
import _pmutils as pmutils

versionFileName = _os.path.join( _os.path.dirname(_os.path.abspath(__file__)), ".VERSION" )
with open(versionFileName) as versionFile:
    version = versionFile.readline()
    version = "".join(version.splitlines())
versionFile.close()

arch = "x86" if "32" in _platform.architecture()[0] else "x64"

class _Struct:
    pass

names = _Struct()
names.paramonte = "ParaMonte"
names.paradram = "ParaDRAM"
names.paranest = "ParaNest"
names.paratemp = "ParaTemp"
