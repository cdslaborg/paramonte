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

class _Struct: pass

####################################################################################################################################

from pathlib import Path as _Path

path = _Struct()
path.root = _os.path.dirname(_os.path.abspath(__file__))
path.auxil = _os.path.join(path.root,"auxil")
path.home = str(_Path.home()) # path.home = _os.path.expanduser("~")
path.lib = path.root

_sys.path.append(path.root)

####################################################################################################################################

platform = _Struct()
platform.arch = "x86" if "32" in _platform.architecture()[0] else "x64"
platform.name = _sys.platform.lower()
platform.isWin32 = True if platform.name=="win32" else False
platform.isLinux = True if platform.name=="linux" else False
platform.isMacOS = True if platform.name=="darwin" else False
platform.osname = "windows" if platform.isWin32 else platform.name

####################################################################################################################################

names = _Struct()
names.paramonte = "ParaMonte"
names.paradram = "ParaDRAM"
names.paradise = "ParaDISE"
names.paranest = "ParaNest"
names.paratemp = "ParaTemp"

####################################################################################################################################

website = _Struct()

website.home = _Struct()
website.home.url = "https://www.cdslab.org/paramonte/"
website.home.install = _Struct()
website.home.install.url = website.home.url + "notes/installation/"
website.home.install.macos = _Struct()
website.home.install.macos.url = website.home.install.url + "macos/"
website.home.install.macos.prereqs = _Struct()
website.home.install.macos.prereqs.url = website.home.install.macos.url + "#the-compile-time-and-runtime-prerequisites"
website.home.install.macos.prereqs.cmd = _Struct()
website.home.install.macos.prereqs.cmd.url = website.home.install.macos.url + "#prereqs-install"

website.github = _Struct()
website.github.issues = _Struct()
website.github.issues.url = "https://github.com/cdslaborg/paramonte/issues"

website.intel = _Struct()
website.intel.mpi = _Struct()
website.intel.mpi.home = _Struct()
website.intel.mpi.home.url = "https://software.intel.com/en-us/mpi-library"

website.intel.mpi.windows = _Struct()
website.intel.mpi.windows.url = "https://software.intel.com/en-us/get-started-with-mpi-for-windows"

website.openmpi = _Struct()
website.openmpi.home = _Struct()
website.openmpi.home.url = "https://www.open-mpi.org/"

####################################################################################################################################

from _Version import Version
version = _Struct()
for versionType in ["interface","kernel"]: setattr(version,versionType,Version(path.auxil,versionType))
