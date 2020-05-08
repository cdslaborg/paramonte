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

import _message as _msg
import numpy as _np
import sys as _sys
import os as _os

class _struct:
    pass

####################################################################################################################################
#### Frozen Struct
####################################################################################################################################

class _FrozenClass(object):
    __isfrozen = False
    def __setattr__(self, key, value):
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError ( "\n{} is a frozen class.\n".format(self)
                            + "The requested attribute '{}' does not exist in the object.\n".format(key)
                            + "You cannot add new attributes to an object of frozen class."
                            )
        object.__setattr__(self, key, value)

    def _freeze(self):
        self.__isfrozen = True


####################################################################################################################################
#### object size estimation
####################################################################################################################################

def getSize(obj, seen=None):
    """Recursively finds size of objects"""
    size = _sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([getSize(v, seen) for v in obj.values()])
        size += sum([getSize(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += getSize(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([getSize(i, seen) for i in obj])
    return size

####################################################################################################################################
#### _Timer class
####################################################################################################################################

import numpy as _np
import time as _time
class Timer:

    def __init__(self,_methodName):
        self._methodName = _methodName
        self.start = 0
        self.end = 0

    def tic(self,msg=None,**noteArgs):
        if msg is not None: 
            if "methodName" not in noteArgs.keys(): noteArgs["methodName"] = self._methodName
            if "end" not in noteArgs.keys(): noteArgs["end"] = ""
            _msg.note( msg = msg, **noteArgs )
        self.start = _time.time()

    def toc(self,msg=""):
        self.end = _time.time()
        if msg=="": msg = "done in " + str(_np.round(self.end-self.start,decimals=6)) + " seconds." 
        if msg is not None: print( msg )

####################################################################################################################################
#### check application installation status
####################################################################################################################################

def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

####################################################################################################################################
#### get file list
####################################################################################################################################

def getFileList(file,fileType,methodName,_mpiDisabled):

    suffix = "_" + fileType + ".txt"
    if _os.path.isfile(file): # check if the input path is a full path to a file
        FileList = [file]
        pattern = file
        if suffix != file[-len(suffix):]:
            _msg.warn   ( msg   = "The name of the input file: \n\n"
                                + "    " + file + "\n\n"
                                + "does not end with the expected suffix '" + suffix + "' for a " + fileType + " file type.\n"
                        , methodName = methodName
                        , marginTop = 1
                        , marginBot = 1
                        )
    elif _os.path.isdir(file): # ensure the input path is not a directory
        _msg.abort  ( msg   = "file='" + file + "' cannot point to a directory.\n"
                            + "Provide a string as the value of file that points to a unique " + fileType + " file or\n"
                            + "to the unique name (including path) of the simulation name shared among its output files.\n"
                    , methodName = methodName
                    , marginTop = 1
                    , marginBot = 1
                    )
    else:

        # search for files matching the input pattern

        import glob
        if file[-1:]=="*":
            pattern = file
        else:
            pattern = file + "*" # + suffix

        _ = glob.glob(pattern)
        FileList = []
        for filename in _:
            if suffix in filename: FileList.append(filename)
        if len(FileList)==0:
            _msg.abort  ( msg   = "Failed to detect any " + fileType + " files with the requested pattern: \n\n"
                                + "    " + pattern + "\n\n"
                                + "Provide a string, as the value of the input argument 'file', that either \n\n"
                                + "    - points to one or more " + fileType + " files, or, \n"
                                + "    - represents the unique name of a ParaMonte simulation. \n"
                                + "      This unique-name is the common prefix in the names of \n"
                                + "      the output files of a ParaMonte simulation."
                        , methodName = methodName
                        , marginTop = 1
                        , marginBot = 1
                        )
        else:
            pattern += suffix

    if _mpiDisabled:
        _msg.note   ( msg = str(len(FileList)) + ' files detected matching the pattern: "' + pattern + '"'
                    , methodName = methodName
                    , marginTop = 0
                    , marginBot = 0
                    )
    return FileList