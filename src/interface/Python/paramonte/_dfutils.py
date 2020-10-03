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

#import typing as tp
from collections.abc import Iterable

class Struct: pass

####################################################################################################################################
#### nam2num
####################################################################################################################################

def nam2num ( nameList  # : tp.Union[ pd.core.indexes.base.Index , tp.List[str] ]
            , names     # : tp.List[str]
            ):
    nml = list(nameList)
    return [ nml.index(name) for name in names ]
    #if isinstance(nameList,list):
    #    return [ nameList.index(name) for name in names ]
    #else:
    #    
    #    return [ nameList.get_loc(name) for name in names ]

####################################################################################################################################
#### getColNamesIndex
####################################################################################################################################

def getColNamesIndex( dfColumns         # : pd.core.indexes.base.Index
                    , columns    = None # : tp.Optional[ tp.Union[ str, range , tp.List[int] , tp.List[str] ] ] = None
                    ):

    allColumnsNeeded = False
    if columns is None:
        allColumnsNeeded = True
    elif isinstance(columns,int):
        columns = [columns,]
    elif isinstance(columns,Iterable):
        if len(columns)==0: allColumnsNeeded = True
    else:
        raise Exception ( "\nThe input \"columns\" argument must be a list whose elements are all\n"
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
        raise Exception ( "\nThe input argument 'columns' must be a list whose elements are all\n"
                        + "    1.   string-valued, each representing the name of the column from the dataframe, or,\n"
                        + "    2.   integer-valued, each representing the index of the column from the dataframe.\n"
                        + "You have entered:\n"
                        + str(columns)
                        )
    return colnames, colindex

####################################################################################################################################
#### getMaxLogFunc
####################################################################################################################################

def getMaxLogFunc( dataFrame # : pd.DataFrame
                 , column = "SampleLogFunc" # : tp.Optional[ str ] = "SampleLogFunc"
                 ):
    """

    Returns a structure with components containing the properties of the 
    input ``dataFrame``, corresponding to the mode (maximum) of the data 
    in the column of the dataFrame that is identified by ``column``. 

        **Parameters**

            dataFrame

                A Pandas dataframe containing the output sample data 
                from a ParaMonte simulation.

            column

                A string containing the name of the column of the input 
                ``dataFrame`` that contains values of the objective 
                function (or its logarithm).

        **Returns**

            maxLogFunc

                A structure with the following components:

                    idrow

                        The ID of the row of ``dataFrame`` corresponding 
                        to the mode of ``column``, 

                    value

                        The value of ``column`` at maximum, 

                    dfrow

                        The entire row of ``dataFrame`` corresponding to 
                        the mode of ``column``, 

                    state

                        The location (state) within the domain of the objective 
                        function where the maximum of ``column`` occurs.


    """

    _offset = dataFrame.columns.get_loc(column) + 1
    maxLogFunc = Struct()
    maxLogFunc.idrow = dataFrame[[column]].idxmax().values[0]
    maxLogFunc.value = dataFrame[[column]].iat[maxLogFunc.idrow,0]
    maxLogFunc.dfrow = dataFrame.iloc[maxLogFunc.idrow,:]
    maxLogFunc.state = dataFrame.iloc[maxLogFunc.idrow,_offset:]
    return maxLogFunc

####################################################################################################################################
