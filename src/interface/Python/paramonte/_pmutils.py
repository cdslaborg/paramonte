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

import _message as err
import numpy as np
import time
import sys
import os

class Struct: pass

newline = chr(10)
creturn = chr(13)

####################################################################################################################################
#### Frozen Struct
####################################################################################################################################

class FrozenClass(object):
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
#### isNumericString
####################################################################################################################################

def isNumericString(string):
    try:
        np.double(string)
        return True
    except ValueError:
        return False

####################################################################################################################################
#### object size estimation
####################################################################################################################################

def getSize(obj, seen=None):
    """Recursively finds size of objects"""
    size = sys.getsizeof(obj)
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
#### check application installation status
####################################################################################################################################

def which(program):
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
#### convert to list
####################################################################################################################################

def getList(x):
    """
    Convert any input object to list.
    """
    if isinstance(x, list):
        return x
    elif isinstance(x, str):
        return [x]
    try:
        return list(x)
    except TypeError:
        return [x]

####################################################################################################################################
#### get file list
####################################################################################################################################

def getFileList(file, fileSuffix, methodName, reportEnabled):

    fullSuffix = "_" + fileSuffix + ".txt"

    if os.path.isfile(file):

        # check if the input path is a full path to a file

        FileList = [file]
        pattern = file
        if fullSuffix != file[-len(fullSuffix):]:
            err.warn( msg   = "The name of the input file: \n\n"
                            + "    " + file + "\n\n"
                            + "does not end with the expected suffix '" + fullSuffix + "' for a " + fileSuffix + " file type.\n"
                    , methodName = methodName
                    , marginTop = 1
                    , marginBot = 1
                    )

    else:

        if os.path.isdir(file):

            # ensure the input path is not a directory

            err.warn( msg   = "file='" + file + "' points to a directory.\n"
                            + "Now searching inside the folder for a " + fileSuffix + " file..."
                    , methodName = methodName
                    , marginTop = 1
                    , marginBot = 1
                    )
            pattern = os.path.join(file, "*" + fullSuffix)

        else:

            # search for files matching the input pattern

            if "*" in file: # file[-1:]=="*":
                pattern = file
            else:
                pattern = file + "*" # + fullSuffix

        import glob
        _ = glob.glob(pattern)

        FileList = []
        for filePath in _:
            if filePath.endswith(fullSuffix):
                FileList.append(filePath)

        if not pattern.endswith(fullSuffix): pattern += fullSuffix

        if len(FileList)==0:
            err.abort   ( msg   = "Failed to detect any " + fileSuffix + " files with the requested pattern: \n\n"
                                + "    " + pattern + "\n\n"
                                + "Provide a string, as the value of the input argument ``file``, that either \n\n"
                                + "    - points to one or more " + fileSuffix + " files, or, \n"
                                + "    - represents the unique name of a ParaMonte simulation. \n"
                                + "      This unique-name is the common prefix in the names of \n"
                                + "      the output files of a ParaMonte simulation.\n\n"
                                + "Most importantly, ensure the requested file is in ASCII format.\n"
                                + "The binary-format chain or restart output files cannot be parsed.\n"
                                + "You can request ASCII-format output files by setting the\n"
                                + "appropriate simulation specifications of the " + methodName + " sampler,\n\n"
                                + "    spec.restartFileFormat = \"ascii\"\n"
                                + "    spec.chainFileFormat = \"ascii\""
                        , methodName = methodName
                        , marginTop = 1
                        , marginBot = 1
                        )
        elif reportEnabled:
            err.note( msg = str(len(FileList)) + ' files detected matching the pattern: "' + pattern + '"'
                    , methodName = methodName
                    , marginTop = 1
                    , marginBot = 1
                    )

    return FileList

####################################################################################################################################
#### Timer class
####################################################################################################################################

class Timer:

    def __init__(self):
        self.start = 0
        self.total = 0
        self.delta = 0
        self.last = 0
        self.tic()

    def tic(self):
        self.start = time.time()
        self.total = self.start
        self.last = self.total
        self.delta = self.total - self.last

    def toc(self):
        self.last = self.total
        self.total = time.time()
        self.delta = self.total - self.last

####################################################################################################################################
#### Progress report
####################################################################################################################################

class Progress:

    def __init__( self
                , msg = None
                , methodName = ""
                , reportEnabled = True
                , end = "\n"
                ):
        self._methodName = methodName
        self._reportEnabled = reportEnabled
        self._oldFraction = 0.0
        self._clockCounter = None
        if self._reportEnabled and msg is not None:
            err.note( msg = msg
                    , methodName = self._methodName
                    , marginTop = 0
                    , marginBot = 0
                    , end = end
                    )
        self.timer = Timer()

    ################################################################################################################################

    def note( self
            , msg = None
            , end = "\n"
            , pre = False
            ):
        self.timer.toc()
        if self._reportEnabled:
            if msg is None: msg = "done in " + str(np.round(self.timer.delta,decimals=6)) + " seconds."
            if pre:
                err.note( msg = msg
                        , methodName = self._methodName
                        , marginTop = 0
                        , marginBot = 0
                        , end = end
                        )
            else:
                print( msg, end = end )
        self.timer.toc()

    ################################################################################################################################

    # dynamic progress bar
    def updateBar   ( self
                    , fraction
                    , progressFraction = 0.05
                    ):
        if self._reportEnabled:
            if fraction<1:
                end = ""
            else:
                end = "\n"
            if fraction > self._oldFraction + progressFraction:
                self._oldFraction = fraction
                print( ".", end=end);

    ################################################################################################################################

    # dynamic clock tic
    def updateClock ( self
                    , fraction
                    ):
        if self._reportEnabled:
            chars = "|/-\\"
            if self._clockCounter is None:
                self._clockCounter = 0
                backspace = "\r"
            else:
                self._clockCounter = self._clockCounter%3 + 1
                backspace = 5*"\r"
            if fraction<1:
                clockTick = chars[self._clockCounter]; 
                end = ""
            else:
                self._clockCounter = None;
                clockTick = chars[0];
                end = "\n"
            print( backspace + clockTick + "{:3.0f}%".format(100*fraction), end=end);

####################################################################################################################################
#### getRandomFilePrefix
####################################################################################################################################

def getRandomFilePrefix(prefix = ""):
    from datetime import datetime as dt
    now = dt.now()
    return  prefix \
            + "{:04d}".format(now.year) + "{:02d}".format(now.month) + "{:02d}".format(now.day) + "_"  \
            + "{:02d}".format(now.hour) + "{:02d}".format(now.minute) + "{:02d}".format(now.second) + "_" \
            + "{:03d}".format(round(now.microsecond/1000))

####################################################################################################################################
#### getLogIntSpace
####################################################################################################################################

def getLogIntSpace(base, logskip, lowerLim, upperLim):
    """
    return a set of unique integer spacings linearly-spaced in the 
    logarithmic scale in the input given base, between the input limits.
    """
    if base<=1: raise ValueError("The input argument \"base\" must be a positive real number. You have entered: " + str(base))
    if logskip<=0: raise ValueError("The input argument \"logskip\" must be a positive real number. You have entered: " + str(logskip))
    if lowerLim<=0: raise ValueError("The input argument \"lowerLim\" must be a positive real number. You have entered: " + str(lowerLim))
    if upperLim<=lowerLim: raise ValueError ( "The input argument \"upperLim\" must be a positive real number larger than "
                                            + "the input argument \"lowerLim\". You have entered: " + str(upperLim))
    return np.unique( np.int32( base**( np.arange( np.log(lowerLim)/np.log(base), np.log(upperLim)/np.log(base), logskip) ) ) );

####################################################################################################################################
