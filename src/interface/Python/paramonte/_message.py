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
####   we ask you to acknowledge the ParaMonte library's usage
####   in your work (education/research/industry/development/...)
####   by citing the ParaMonte library as described on this page:
####
####       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
####
####################################################################################################################################
####################################################################################################################################

def informUser  ( msg
                , prefix = ""
                , newLine = "\n"
                , marginTop = 0
                , marginBot = 0
                , end = "\n"
                ):
    for i in range(marginTop): print()
    print(prefix + msg.replace(newLine,newLine+prefix), end=end)
    for i in range(marginBot): print()
    return None

def note( msg
        , methodName = ""
        , newLine = "\n"
        , marginTop = 0
        , marginBot = 0
        , end = "\n"
        ):
    informUser  ( msg
                , prefix = methodName + " - NOTE: "
                , newLine = newLine
                , marginTop = marginTop
                , marginBot = marginBot
                , end = end
                )
    return None

def warn( msg
        , methodName = ""
        , newLine = "\n"
        , marginTop = 0
        , marginBot = 0
        , end = "\n"
        ):
    informUser  ( msg
                , prefix = methodName + " - WARNING: "
                , newLine = newLine
                , marginTop = marginTop
                , marginBot = marginBot
                , end = end
                )
    return None

def abort   ( msg
            , methodName = ""
            , newLine = "\n"
            , marginTop = 0
            , marginBot = 0
            , end = "\n"
            ):
    informUser  ( msg
                , prefix = methodName + " - FATAL: "
                , newLine = newLine
                , marginTop = marginTop
                , marginBot = marginBot
                , end = end
                )
    import sys as _sys
    _sys.exit()
    return None
