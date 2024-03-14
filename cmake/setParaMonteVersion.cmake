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

#   This script defines the variable `ParaMonteVersion` with the value from the file `VERSION.md`.

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Parse version from VERSION.md file so that more info can be added and easier to get from scripts
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

file(STRINGS "${paramonte_src_${lang}_dir}/VERSION.md" first_line LIMIT_COUNT 1)
string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+(-rc[0-9]+)?" ParaMonteVersion "${first_line}")
if(NOT (ParaMonteVersion MATCHES "[0-9]+\\.[0-9]+\\.[0-9]+(-rc[0-9]+)?"))
    message(FATAL_ERROR
            "\n"
            "${pmfatal} Could not extract version from git, falling back on VERSION.md, line 3."
            "\n"
            )
endif()
string(REGEX REPLACE "-rc[0-9]+$" ".0" ParaMonteVersion "${ParaMonteVersion}")
message(NOTICE "${pmattn} ParaMonte::${lang} version: ${ParaMonteVersion}")
string(REPLACE "." ";" ParaMonteVersionList "${ParaMonteVersion}")
list(GET ParaMonteVersionList 0 ParaMonteVersionMajor)
list(GET ParaMonteVersionList 1 ParaMonteVersionMinor)
list(GET ParaMonteVersionList 2 ParaMonteVersionPatch)
