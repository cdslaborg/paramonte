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

import os
import sys
import platform as _platform

import _dfutils as dfutils
import _pmutils as utils
from _message import note, warn, abort
from _pmutils import Struct, newline, creturn

####################################################################################################################################

from pathlib import Path as _Path

path = Struct()
path.root = os.path.dirname(os.path.abspath(__file__))
path.auxil = os.path.join(path.root,"auxil")
path.home = str(_Path.home()) # path.home = os.path.expanduser("~")
path.lib = path.root

sys.path.append(path.root)

####################################################################################################################################

names = Struct()
names.paramonte = "ParaMonte"
names.paradram = "ParaDRAM"
names.paradise = "ParaDISE"
names.paranest = "ParaNest"
names.paratemp = "ParaTemp"

####################################################################################################################################

platform = Struct()
platform.arch = "x86" if "32" in _platform.architecture()[0] else "x64"
platform.name = sys.platform.lower()
platform.isWin32 = True if platform.name=="win32" else False
platform.isLinux = True if platform.name=="linux" else False
platform.isMacOS = True if platform.name=="darwin" else False
platform.osname = "windows" if platform.isWin32 else platform.name

from datetime import datetime as _dt
dt = _dt.now()
platform.systemInfoFilePrefix = os.path.join(path.auxil, ".systemInfo_");
platform.systemInfoFilePath = platform.systemInfoFilePrefix + "{:04d}".format(dt.year) + "{:02d}".format(dt.month) + "{:02d}".format(dt.day);

if not os.path.isfile(platform.systemInfoFilePath):

    import glob
    fileList = glob.glob(platform.systemInfoFilePrefix + "*")
    for filePath in fileList:
        try:
            os.remove(filePath)
        except:
            pass

    cmd = None
    if platform.isWin32: cmd = "systeminfo > " + platform.systemInfoFilePath
    if platform.isLinux: cmd = "uname -a >>  " + platform.systemInfoFilePath + "; lscpu >> " + platform.systemInfoFilePath
    if platform.isMacOS: cmd = "uname -a >>  " + platform.systemInfoFilePath + "; sysctl -a | grep machdep.cpu >> " + platform.systemInfoFilePath

    if cmd is not None:
        err = os.system(cmd)
        if err != 0:
            platform.systemInfoFilePath = None
            warn( msg   = "Failed to get the system information. skipping..."
                , methodName = names.paramonte
                , marginTop = 1
                , marginBot = 1
                )

####################################################################################################################################

website = Struct()

website.home = Struct()
website.home.url = "https://www.cdslab.org/paramonte"
website.home.install = Struct()
website.home.install._url = website.home.url + "/notes/installation"

# installation Linux

website.home.overview = Struct()
website.home.overview._url = website.home.url + "/notes/overview"
website.home.overview.preface = Struct()
website.home.overview.changes = Struct()
website.home.overview.preface.url = website.home.overview._url + "/preface"
website.home.overview.changes.kernel = Struct()
website.home.overview.changes.matlab = Struct()
website.home.overview.changes.python = Struct()
website.home.overview.changes.kernel.url = website.home.overview._url + "/paramonte-kernel-release-notes"
website.home.overview.changes.matlab.url = website.home.overview._url + "/paramonte-matlab-release-notes"
website.home.overview.changes.python.url = website.home.overview._url + "/paramonte-python-release-notes"

# installation Linux

website.home.install.linux = Struct()
website.home.install.linux.url = website.home.install._url + "/linux"

# installation Windows

website.home.install.windows = Struct()
website.home.install.windows.url = website.home.install._url + "/windows"

# installation Python

website.home.install.python = Struct()
website.home.install.python.url = website.home.install._url + "/python"

# installation MATLAB

website.home.install.matlab = Struct()
website.home.install.matlab.url = website.home.install._url + "/matlab"

# installation macOS

website.home.install.macos = Struct()
website.home.install.macos.url = website.home.install._url + "/macos"
website.home.install.macos.prereqs = Struct()
website.home.install.macos.prereqs.url = website.home.install.macos.url + "/#the-compile-time-and-runtime-prerequisites"
website.home.install.macos.prereqs.cmd = Struct()
website.home.install.macos.prereqs.cmd.url = website.home.install.macos.url + "/#prereqs-install"

# Python examples

website.home.examples = Struct()
website.home.examples._url = website.home.url + "/notes/examples"
website.home.examples.python = Struct()
website.home.examples.python.jupyter = Struct()
website.home.examples.python.postprocess = Struct()
website.home.examples.python.jupyter.url = website.home.examples._url + "/python/jupyter"
website.home.examples.python.postprocess.url = website.home.examples._url + "/python/postprocess"

# Python API

website.home.api = Struct()
website.home.api._url = website.home.url + "/notes/api"
website.home.api.python = Struct()
website.home.api.python.url = website.home.api._url + "/python/autoapi/paramonte"

# ParaDRAM

website.home.usage = Struct()
website.home.usage._url = website.home.url + "/notes/usage"
website.home.usage.paradram = Struct()
website.home.usage.paradram._url = website.home.usage._url + "/paradram"
website.home.usage.paradram.quickstart = Struct()
website.home.usage.paradram.quickstart.url = website.home.usage.paradram._url + "/interface"
website.home.usage.paradram.input = Struct()
website.home.usage.paradram.input.url = website.home.usage.paradram._url + "/input"
website.home.usage.paradram.specifications = Struct()
website.home.usage.paradram.specifications.url = website.home.usage.paradram._url + "/specifications"
website.home.usage.paradram.restart = Struct()
website.home.usage.paradram.restart.url = website.home.usage.paradram._url + "/restart"
website.home.usage.paradram.output = Struct()
website.home.usage.paradram.output.url = website.home.usage.paradram._url + "/output"

# GitHub issues

website.github = Struct()
website.github.url = "https://github.com/cdslaborg/paramonte"
website.github.issues = Struct()
website.github.issues.url = "https://github.com/cdslaborg/paramonte/issues"
website.github.release = Struct()
website.github.release.url = website.github.url + "/releases"
website.github.release.latest = Struct()
website.github.release.latest.url = website.github.release.url + "/latest"
website.github.archive = Struct()
website.github.archive._url = website.github.url + "/archive"
website.github.archive.master = Struct()
website.github.archive.master.zip = Struct()
website.github.archive.master.tar = Struct()
website.github.archive.master.zip.url = website.github.archive._url + "/master.zip"
website.github.archive.master.tar.url = website.github.archive._url + "/master.tar.gz"

# GitHub examples

website.github.examples = Struct()
website.github.examples.url = "https://github.com/cdslaborg/paramontex"

# Intel MPI

website.intel = Struct()
website.intel.mpi = Struct()
website.intel.mpi.home = Struct()
website.intel.mpi.home.url = "https://software.intel.com/en-us/mpi-library"

# Intel Windows

website.intel.mpi.windows = Struct()
website.intel.mpi.windows.url = "https://software.intel.com/en-us/get-started-with-mpi-for-windows"

# OpenMPI

website.openmpi = Struct()
website.openmpi.home = Struct()
website.openmpi.home.url = "https://www.open-mpi.org/"

####################################################################################################################################

def cite(): print(website.home.overview.preface.url + "/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work")
#citation.kernel = Struct()
#citation.kernel.paradram = Struct()
#citation.kernel.paradram.bib = """
#@article{2020arXiv200809589S,
#           author = {{Shahmoradi}, Amir and {Bagheri}, Fatemeh},
#            title = "{ParaDRAM: A Cross-Language Toolbox for Parallel High-Performance Delayed-Rejection Adaptive Metropolis Markov Chain Monte Carlo Simulations}",
#          journal = {arXiv e-prints},
#         keywords = {Computer Science - Computational Engineering, Finance, and Science, Astrophysics - Instrumentation and Methods for Astrophysics, Physics - Data Analysis, Statistics and Probability, Statistics - Computation, Statistics - Machine Learning},
#             year = 2020,
#            month = aug,
#              eid = {arXiv:2008.09589},
#            pages = {arXiv:2008.09589},
#    archivePrefix = {arXiv},
#           eprint = {2008.09589},
#     primaryClass = {cs.CE},
#           adsurl = {https://ui.adsabs.harvard.edu/abs/2020arXiv200809589S},
#          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
#}
#"""

####################################################################################################################################

from _Version import Version
version = Struct()
for versionType in ["interface","kernel"]: setattr(version,versionType,Version(path.auxil,versionType))

####################################################################################################################################

