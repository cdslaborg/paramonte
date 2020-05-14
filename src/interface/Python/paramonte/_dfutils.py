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

import numpy as _np
import typing as _tp
import pandas as _pd
import _message as msg
from collections.abc import Iterable as _Iterable

####################################################################################################################################

def nam2num ( nameList  : _tp.Union[ _pd.core.indexes.base.Index , _tp.List[str] ]
            , names     : _tp.List[str]
            ):
    nml = list(nameList)
    return [ nml.index(name) for name in names ]
    #if isinstance(nameList,list):
    #    return [ nameList.index(name) for name in names ]
    #else:
    #    
    #    return [ nameList.get_loc(name) for name in names ]

####################################################################################################################################

def getColNamesIndex( dfColumns     : _pd.core.indexes.base.Index
                    , columns       : _tp.Optional[ _tp.Union[ str, range , _tp.List[int] , _tp.List[str] ] ] = None
                    ):

    allColumnsNeeded = False
    if columns is None:
        allColumnsNeeded = True
    elif isinstance(columns,int):
        columns = [columns,]
    elif isinstance(columns,_Iterable):
        if len(columns)==0: allColumnsNeeded = True
    else:
        raise Exception ( "The input argument 'columns' must be a list whose elements are all\n"
                        + "    1.   string-valued, each representing the name of the column from the dataframe, or,\n"
                        + "    2.   integer-valued, each representing the index of the column from the dataframe.\n"
                        + "You have entered:\n"
                        + str(columns)
                        )

    if allColumnsNeeded:
        colnames = dfColumns
        colindex = range(len(colnames))
    elif isinstance(columns, str):
        colnames = [columns]
        colindex = nam2num( dfColumns , [columns,] )
    elif all(isinstance(element, str) for element in columns):
        colnames = columns
        colindex = nam2num( dfColumns , columns )
    elif all(isinstance(element,int) for element in columns):
        colindex = columns
        colnames = dfColumns[colindex]
    #elif all(isinstance(element,float) for element in columns):
    #    colindex = columns
    #    colnames = None
    else:
        raise Exception ( "The input argument 'columns' must be a list whose elements are all\n"
                        + "    1.   string-valued, each representing the name of the column from the dataframe, or,\n"
                        + "    2.   integer-valued, each representing the index of the column from the dataframe.\n"
                        + "You have entered:\n"
                        + str(columns)
                        )
    return colnames, colindex

####################################################################################################################################
