#!/usr/bin/env python
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
####       https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md
####
####################################################################################################################################
####################################################################################################################################

# ATTN: This code must be executed from the base directory where it exists. Navigate the the directory containing this file, then: python getReadme.py
# This code parses the code coverage percentage reports from the web and builds the latest version of shield.html based on imported percentages.

import sys
import urllib.request
urlBase = "https://cdslaborg.github.io/paramonte-codecov/kernel/"

def genShield():

    shields = """
<a href="https://github.com/cdslaborg/paramonte#license" target="_blank"><img src="https://img.shields.io/github/license/cdslaborg/paramonte?color=orange&style=flat-square" alt="GitHub" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/v/release/cdslaborg/paramonte?color=orange&label=kernel%20release&style=flat-square" alt="GitHub release (latest by date)" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/release-date/cdslaborg/paramonte?color=orange&style=flat-square" alt="GitHub Release Date" /></a>
<a href="https://pypi.org/project/paramonte/" target="_blank"><img src="https://img.shields.io/pypi/v/paramonte?color=orange&label=pypi%20release&style=flat-square" alt="PyPI - release version" /></a>
<a href="https://travis-ci.com/cdslaborg/paramonte" target="_blank"><img src="https://travis-ci.com/cdslaborg/paramonte.svg?branch=main&style=flat-square" alt="Build Status" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/pypi/status/paramonte?style=flat-square" alt="PyPI - Status" /></a>
<a href="https://lgtm.com/projects/g/cdslaborg/paramonte/?mode=list" target="_blank"><img src="https://img.shields.io/lgtm/grade/python/github/cdslaborg/paramonte?label=code%20quality&style=flat-square&color=brightgreen" alt="LGTM Grade" /></a>
<a href="https://cdslaborg.github.io/paramonte-codecov/kernel/serial/" target="_blank"><img src="https://img.shields.io/badge/code%20coverage-serial%20kernel-brightgreen?style=flat-square" alt="kernel code coverage - serial" /></a>
<a href="https://cdslaborg.github.io/paramonte-codecov/kernel/mpi/"    target="_blank"><img src="https://img.shields.io/badge/code%20coverage-MPI%20kernel-brightgreen?style=flat-square" alt="kernel code coverage - MPI" /></a>
<a href="https://cdslaborg.github.io/paramonte-codecov/kernel/caf/"    target="_blank"><img src="https://img.shields.io/badge/code%20coverage-Coarray%20kernel-brightgreen?style=flat-square" alt="kernel code coverage - Coarray" /></a>
<a href="https://www.openhub.net/p/paramonte" target="_blank"><img src="https://img.shields.io/badge/Open%20Hub-stats?color=brightgreen&label=stats&message=Open%20Hub&style=flat-square" alt="stats - Open Hub" /></a>
<a href="https://github.com/cdslaborg/paramonte/issues" target="_blank"><img src="https://img.shields.io/github/issues/cdslaborg/paramonte?style=flat-square" alt="GitHub issues" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/main/src/interface" target="_blank"><img src="https://img.shields.io/badge/available%20in-C%20%2F%20C%2B%2B%20%2F%20Fortran%20%2F%20MATLAB%20%2F%20Python-brightgreen?style=flat-square" alt="supported languages" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/traffic" target="_blank"><img src="https://img.shields.io/github/downloads/cdslaborg/paramonte/total?color=brightgreen&label=kernel%20downloads&style=flat-square" alt="GitHub All Releases" /></a>
<a href="https://libraries.io/pypi/paramonte" target="_blank"><img src="https://img.shields.io/pypi/dm/paramonte?color=brightgreen&label=pypi%20downloads&style=flat-square" alt="PyPI - Downloads" /></a>
<a href="https://pepy.tech/project/paramonte" target="_blank"><img src="https://static.pepy.tech/personalized-badge/paramonte?period=month&units=international_system&left_color=grey&right_color=brightgreen&left_text=pypi%20downloads" alt="PyPI - Downloads" /></a>
<a href="https://pepy.tech/project/paramonte" target="_blank"><img src="https://img.shields.io/badge/dynamic/json?style=flat-square&labelColor=brightgreen&color=grey&maxAge=86400&label=pypi%20downloads&query=%24.total_downloads&url=https%3A%2F%2Fapi.pepy.tech%2Fapi%2Fprojects%2Fparamonte" alt="PyPI - Downloads" /></a>
<a href="https://www.mathworks.com/matlabcentral/fileexchange/78946-paramonte" target="_blank"><img src="https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg" alt="View ParaMonte on File Exchange" /></a>
<a href="https://github.com/cdslaborg/paramonte/" target="_blank"><img src="https://img.shields.io/github/repo-size/cdslaborg/paramonte?style=flat-square" alt="GitHub repo size" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/main/src/interface" target="_blank"><img src="https://img.shields.io/github/languages/count/cdslaborg/paramonte?style=flat-square" alt="GitHub language count" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/contributors" target="_blank"><img src="https://img.shields.io/github/commit-activity/y/cdslaborg/paramonte?style=flat-square" alt="GitHub commit activity" /></a>
<a href="https://github.com/cdslaborg/paramonte/commits/main" target="_blank"><img src="https://img.shields.io/github/last-commit/cdslaborg/paramonte?color=blue&style=flat-square" alt="GitHub last commit" /></a>
<a href="https://zenodo.org/record/4076479#.X4Stte17ng4" target="_blank"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.4076479.svg" alt="citations and references" /></a>
<a href="https://www.cdslab.org/paramonte/notes/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work" target="_blank"><img src="https://img.shields.io/badge/reference-%20%09arXiv%3A1209.4647-blueviolet?style=flat-square" alt="citations and references" /></a>
<a href="https://ascl.net/2008.016" target="_blank"><img src="https://img.shields.io/badge/ascl-2008.016-blue.svg?colorB=262255" alt="ascl:2008.016" /></a>
<a style="border-width:0" href="https://doi.org/10.21105/joss.02741"><img src="https://joss.theoj.org/papers/10.21105/joss.02741/status.svg?style=flat-square" alt="DOI badge" ></a>
<br><br>
<a href="https://twitter.com/intent/tweet?text=ParaMonte%20-%20Plain%20Powerfull%20Parallel%20Monte%20Carlo%20Library:&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" target="_blank"><img src="https://img.shields.io/twitter/url?style=social&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" alt="Twitter" /></a>
"""
#<a href="https://joss.theoj.org/papers/f964b6e22c71515c310fbe3843ad4513"><img src="https://joss.theoj.org/papers/f964b6e22c71515c310fbe3843ad4513/status.svg"></a>

    parDict =   { "serial"  : "serial"
                , "mpi"     : "MPI"
                , "caf"     : "Coarray"
                }

    for par in parDict.keys():
        try:
            request = urllib.request.Request(urlBase + par)
            response = urllib.request.urlopen(request)
            html = response.read().decode("utf8")
            coverage = float(html.split('headerCovTableEntryHi">')[1].split("%")[0])
            search = "code%20coverage-" + parDict[par] + "%20kernel"
            #substitute = "code%20coverage%20%2d%20" + parDict[par] + "%20kernel-" + str(coverage) + "%25"
            substitute = search + "%20:%20" + str(coverage) + "%25"
            print(search)
            print(substitute)
            shields = shields.replace(search, substitute)
            print(parDict[par] + " code coverage {}%".format(coverage))
        except Exception as errmsg:
            print   ( "WARNING: Failed to fetch the ParaMonte Kernel code coverage percentage from the web.\n"
                    + "WARNING: Here is the full exception message:\n\n"
                    + str(errmsg)
                    + "\n\n"
                    + "WARNING: skipping the code coverage import into the shields html.")
    shields = shields.strip("\n")
    with open("shields.html", "w") as file: file.write(shields)
    return shields

if __name__ == "__main__": genShield()