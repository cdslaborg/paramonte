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

#############################################################################
# Given a list of flags, this function will try each, one at a time,
# and choose the first flag that works.  If no flags work, then nothing
# will be set, unless the REQUIRED key is given, in which case an error
# will be given.
# 
# Call is:
# SET_COMPILE_FLAG(FLAGVAR FLAGVAL (Fortran|C|CXX) <REQUIRED> flag1 flag2...)
# 
# For example, if you have the flag CMAKE_C_FLAGS and you want to add
# warnings and want to fail if this is not possible, you might call this
# function in this manner:
# SET_COMPILE_FLAGS(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}" C REQUIRED
#                   "-Wall"     # GNU
#                   "-warn all" # Intel
#                  )
# The optin "-Wall" will be checked first, and if it works, will be
# appended to the CMAKE_C_FLAGS variable.  If it doesn't work, then
# "-warn all" will be tried.  If this doesn't work then checking will
# terminate because REQUIRED was given.  
#
# The reasong that the variable must be given twice (first as the name then
# as the value in quotes) is because of the way CMAKE handles the passing
# of variables in functions; it is difficult to extract a variable's
# contents and assign new values to it from within a function.
#############################################################################

INCLUDE(${CMAKE_ROOT}/Modules/CheckCCompilerFlag.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CheckCXXCompilerFlag.cmake)

FUNCTION(SET_COMPILE_FLAG FLAGVAR FLAGVAL LANG)

    # Do some up front setup if Fortran
    IF(LANG STREQUAL "Fortran")
        # Create a list of error messages from compilers
        SET(FAIL_REGEX
            "unrecognized source type"              # Intel, example: "ifort: command line warning #10161: unrecognized source type 'full'; object file assumed"
            "unrecognized option"                   # Intel, example: "ifort: LINK : warning LNK4044: unrecognized option '/F999999999'; ignored"
           #"command line warning"                  # Intel, example: "ifort: command line warning #10161: unrecognized source type 'full'; object file assumed"
            "command line error"                    # Intel, example: "ifort: command line error: option '/g' is ambiguous"
            "ignoring option"                       # Intel, example: "ifort: command line warning #10157: ignoring option '/W'; argument is of wrong type"
            "ignoring unknown option"               # Intel
            "no action performed"                   # Intel, example: "-check all" on Windows
            "invalid argument"                      # Intel
            "unrecognized .*option"                 # GNU
            "[Uu]nknown switch"                     # Portland Group
            "ignoring unknown option"               # MSVC
            "warning D9002"                         # MSVC, any lang
            "[Uu]nknown option"                     # HP
            "[Ww]arning: [Oo]ption"                 # SunPro
            "command option .* is not recognized"   # XL
           )
    ENDIF(LANG STREQUAL "Fortran")

    # Make a variable holding the flags.  Filter out REQUIRED if it is there
    SET(FLAG_REQUIRED FALSE)
    SET(FLAG_FOUND FALSE)
    UNSET(FLAGLIST)
    FOREACH (var ${ARGN})
        STRING(TOUPPER "${var}" UP)
        IF(UP STREQUAL "REQUIRED")
            SET(FLAG_REQUIRED TRUE)
        ELSE()
            SET(FLAGLIST ${FLAGLIST} "${var}")
        ENDIF(UP STREQUAL "REQUIRED")
    ENDFOREACH (var ${ARGN})

    # Now, loop over each flag
    FOREACH(flag ${FLAGLIST})

        UNSET(FLAG_WORKS)
        # Check the flag for the given language
        IF(LANG STREQUAL "C")
            CHECK_C_COMPILER_FLAG("${flag}" FLAG_WORKS)
        ELSEIF(LANG STREQUAL "CXX")
            CHECK_CXX_COMPILER_FLAG("${flag}" FLAG_WORKS)
        ELSEIF(LANG STREQUAL "Fortran")
            # There is no nice function to do this for FORTRAN, so we must manually
            # create a test program and check if it compiles with a given flag.
            SET(TESTFILE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}")
            SET(TESTFILE "${TESTFILE}/CMakeTmp/testFortranFlags.f90")
            FILE(WRITE "${TESTFILE}"
"
program dummyprog
    i = 5
end program dummyprog
")
            message ( STATUS "testing flag ${flag}" )
            TRY_COMPILE(FLAG_WORKS ${CMAKE_BINARY_DIR} ${TESTFILE}
                COMPILE_DEFINITIONS "${flag}" OUTPUT_VARIABLE OUTPUT)
            
            # Check that the output message doesn't match any errors
            FOREACH(rx ${FAIL_REGEX})
                IF("${OUTPUT}" MATCHES "${rx}")
                    SET(FLAG_WORKS FALSE)
                ENDIF("${OUTPUT}" MATCHES "${rx}")
            ENDFOREACH(rx ${FAIL_REGEX})

        ELSE()
            MESSAGE(FATAL_ERROR "Unknown language in SET_COMPILE_FLAGS: ${LANG}")
        ENDIF(LANG STREQUAL "C")

        # If this worked, use these flags, otherwise use other flags
        IF(FLAG_WORKS)
            # Append this flag to the end of the list that already exists
            SET(${FLAGVAR} "${FLAGVAL} ${flag}" CACHE STRING
                 "Set the ${FLAGVAR} flags" FORCE)
            SET(FLAG_FOUND TRUE)
            MESSAGE(STATUS "${BoldGreen}Adding flag ${flag}" )
            BREAK() # We found something that works, so exit
        ELSE()
            MESSAGE(STATUS "${BoldRed}Discarding flag ${flag}" )
        ENDIF(FLAG_WORKS)

    ENDFOREACH(flag ${FLAGLIST})

    # Raise an error if no flag was found
    IF(FLAG_REQUIRED AND NOT FLAG_FOUND)
        MESSAGE(FATAL_ERROR "No compile flags were found")
    ENDIF(FLAG_REQUIRED AND NOT FLAG_FOUND)

ENDFUNCTION()
