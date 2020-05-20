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

import os as _os
import sys as _sys
import numpy as _np
import typing as _tp
import pandas as _pd
import platform as _platform

from _message import note, warn, abort
import _visualization as vis
import _statistics as stats
import _dfutils as dfutils
import _pmutils as pmutils

####################################################################################################################################

class _Struct:
    pass

####################################################################################################################################

from pathlib import Path as _Path
path = _Struct()
path.root = _os.path.dirname(_os.path.abspath(__file__))
path.auxil = _os.path.join(path.root,"auxil")
#path.home = _os.path.expanduser("~")
path.home = str(_Path.home())
path.lib = path.root

_sys.path.append(path.root)

####################################################################################################################################

arch = "x86" if "32" in _platform.architecture()[0] else "x64"

####################################################################################################################################

names = _Struct()
names.paramonte = "ParaMonte"
names.paradram = "ParaDRAM"
names.paranest = "ParaNest"
names.paratemp = "ParaTemp"

####################################################################################################################################

from _Version import Version
version = _Struct()
for versionType in ["interface","kernel"]: setattr(version,versionType,Version(path.auxil,versionType))
