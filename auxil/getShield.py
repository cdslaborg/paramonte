#!/usr/bin/env python
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

# ATTN: This code must be executed from the base directory where it exists. Navigate the directory containing this file, then: python getReadme.py
# This code parses the code coverage percentage reports from the web and builds the latest version of shield.html based on imported percentages.

import sys
import urllib.request

url_codecov_base = "https://www.cdslab.org/paramonte/codecov"
url_codecov = {"fortran" : url_codecov_base + "/fortran/2"}

# <a href="https://github.com/cdslaborg/paramonte#license" target="_blank"><img src="https://img.shields.io/github/license/cdslaborg/paramonte?color=orange&style=flat-square" alt="GitHub" /></a>
# <a href="https://www.cdslab.org/paramonte/codecov/fortran/1/mpi/" target="_blank"><img src="https://img.shields.io/badge/Fortran%20code%20coverage-MPI-brightgreen?style=flat-square" alt="Fortran code coverage - MPI" /></a>
# <a href="https://www.cdslab.org/paramonte/codecov/fortran/1/caf/" target="_blank"><img src="https://img.shields.io/badge/Fortran%20code%20coverage-Coarray-brightgreen?style=flat-square" alt="Fortran code coverage - Coarray" /></a>

def getShield():

    shields = """
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/release-date/cdslaborg/paramonte?color=orange&style=flat-square" alt="GitHub Release Date" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/v/release/cdslaborg/paramonte?color=purple&label=release%20version&style=flat-square" alt="GitHub release (latest by date)" /></a>
<a href="https://pypi.org/project/paramonte/" target="_blank"><img src="https://img.shields.io/pypi/v/paramonte?color=orange&label=pypi%20release&style=flat-square" alt="PyPI - release version" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/pypi/status/paramonte?style=flat-square" alt="PyPI - Status" /></a>
<a href="https://www.cdslab.org/paramonte/codecov/fortran/2/serial/" target="_blank"><img src="https://img.shields.io/badge/Fortran%20code%20coverage-serial-brightgreen?style=flat-square" alt="Fortran code coverage - serial" /></a>
<a href="https://github.com/cdslaborg/paramonte/issues" target="_blank"><img src="https://img.shields.io/github/issues/cdslaborg/paramonte?style=flat-square" alt="GitHub issues" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/main/src" target="_blank"><img src="https://img.shields.io/badge/available%20in-C%20%2F%20C%2B%2B%20%2F%20Fortran%20%2F%20MATLAB%20%2F%20Python-brightgreen?style=flat-square" alt="supported languages" /></a>
<a href="https://www.openhub.net/p/paramonte" target="_blank"><img src="https://img.shields.io/badge/Open%20Hub-stats?color=brightgreen&label=stats&message=Open%20Hub&style=flat-square" alt="stats - Open Hub" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/traffic" target="_blank"><img src="https://img.shields.io/github/downloads/cdslaborg/paramonte/total?color=brightgreen&label=GitHub%20downloads&style=flat-square" alt="GitHub All Releases" /></a>
<a href="https://libraries.io/pypi/paramonte" target="_blank"><img src="https://img.shields.io/pypi/dm/paramonte?color=brightgreen&label=PyPI%20downloads&style=flat-square" alt="PyPI - Downloads" /></a>
<a href="https://pypistats.org/packages/paramonte" target="_blank"><img src="https://img.shields.io/badge/stats-green?style=flat-square&label=PyPI&labelColor=grey&color=brightgreen" alt="PyPI stats" /></a>
<a href="https://www.mathworks.com/matlabcentral/fileexchange/78946-paramonte" target="_blank"><img src="https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg" alt="View ParaMonte on File Exchange" /></a>
<a href="https://github.com/cdslaborg/paramonte/" target="_blank"><img src="https://img.shields.io/github/repo-size/cdslaborg/paramonte?style=flat-square" alt="GitHub repo size" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/main/src" target="_blank"><img src="https://img.shields.io/github/languages/count/cdslaborg/paramonte?style=flat-square" alt="GitHub language count" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/contributors" target="_blank"><img src="https://img.shields.io/github/commit-activity/y/cdslaborg/paramonte?style=flat-square" alt="GitHub commit activity" /></a>
<a href="https://github.com/cdslaborg/paramonte/commits/main" target="_blank"><img src="https://img.shields.io/github/last-commit/cdslaborg/paramonte?color=blue&style=flat-square" alt="GitHub last commit" /></a>
<a href="https://zenodo.org/record/4076479#.X4Stte17ng4" target="_blank"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.4076479.svg" alt="citations and references" /></a>
<a href="https://www.cdslab.org/paramonte/generic/latest/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work" target="_blank"><img src="https://img.shields.io/badge/reference-%20%09arXiv%3A1209.4647-blueviolet?style=flat-square" alt="citations and references" /></a>
<a href="https://ascl.net/2008.016" target="_blank"><img src="https://img.shields.io/badge/ascl-2008.016-blue.svg?colorB=262255" alt="ascl:2008.016" /></a>
<a style="border-width:0" href="https://doi.org/10.21105/joss.02741"><img src="https://joss.theoj.org/papers/10.21105/joss.02741/status.svg?style=flat-square" alt="DOI badge" ></a>
<br><br>
<a href="https://twitter.com/intent/tweet?text=ParaMonte%20-%20Plain%20Powerfull%20Parallel%20Monte%20Carlo%20Library:&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" target="_blank"><img src="https://img.shields.io/twitter/url?style=social&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" alt="Twitter" /></a>
"""
#<a href="https://joss.theoj.org/papers/f964b6e22c71515c310fbe3843ad4513"><img src="https://joss.theoj.org/papers/f964b6e22c71515c310fbe3843ad4513/status.svg"></a>

    parDict =   { "serial"  : "serial"
               #, "mpi"     : "MPI"
               #, "caf"     : "Coarray"
                }
    for par in parDict.keys():
        try:
            request = urllib.request.Request(url_codecov["fortran"] + "/" + par)
            response = urllib.request.urlopen(request)
            html = response.read().decode("utf8")
            for id in ["headerCovTableEntryLow", "headerCovTableEntryMed", "headerCovTableEntryHi"]:
                ids = id + '">'
                try: 
                    coverage = float(html.split(ids)[1].split("%")[0])
                    break
                except Exception as errmsg:
                    errmsg = str(errmsg)
                    coverage = None
            if coverage is None:
                print   ( "WARNING: Failed to fetch the ParaMonte Fortran code coverage percentage from the web.\n"
                        + "WARNING: Here is the full exception message:\n\n"
                        + errmsg
                        + "\n\n"
                        + "WARNING: skipping the code coverage import into the shields html.")
            else:
                search = "Fortran%20code%20coverage-" + parDict[par]
                #substitute = "code%20coverage%20%2d%20" + parDict[par] + "%20fortran-" + str(coverage) + "%25"
                substitute = search + "%20:%20" + str(coverage) + "%25"
                print(search)
                print(substitute)
                shields = shields.replace(search, substitute)
                print(parDict[par] + " code coverage {}%".format(coverage))
        except Exception as errmsg:
            print   ( "WARNING: Failed to fetch the ParaMonte Fortran code coverage percentage from the web.\n"
                    + "WARNING: Here is the full exception message:\n\n"
                    + str(errmsg)
                    + "\n\n"
                    + "WARNING: skipping the code coverage import into the shields html.")
    shields = shields.strip("\n")
    with open("shields.html", "w") as file: file.write(shields)
    return shields

if __name__ == "__main__": getShield()