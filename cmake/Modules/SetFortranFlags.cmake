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

######################################################
# Determine and set the Fortran compiler flags 
######################################################

####################################################################
# Make sure that the default build type is RELEASE if not specified.
####################################################################
INCLUDE(${CMAKE_MODULE_PATH}/SetCompileFlag.cmake)

# Make sure the build type is uppercase
STRING(TOUPPER "${BTYPE}" BT)

IF(BT STREQUAL "RELEASE")
    SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
      FORCE)
ELSEIF(BT STREQUAL "DEBUG")
    SET (CMAKE_BUILD_TYPE DEBUG CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
      FORCE)
ELSEIF(BT STREQUAL "TESTING")
    SET (CMAKE_BUILD_TYPE TESTING CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
      FORCE)
ELSEIF(NOT BT)
    SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
      FORCE)
    MESSAGE(STATUS "CMAKE_BUILD_TYPE not given, defaulting to RELEASE")
ELSE()
    MESSAGE(FATAL_ERROR "CMAKE_BUILD_TYPE not valid, choices are DEBUG, RELEASE, or TESTING")
ENDIF(BT STREQUAL "RELEASE")

#########################################################
# If the compiler flags have already been set, return now
#########################################################

IF(CMAKE_Fortran_FLAGS_RELEASE AND CMAKE_Fortran_FLAGS_TESTING AND CMAKE_Fortran_FLAGS_DEBUG)
    RETURN ()
ENDIF(CMAKE_Fortran_FLAGS_RELEASE AND CMAKE_Fortran_FLAGS_TESTING AND CMAKE_Fortran_FLAGS_DEBUG)

########################################################################
# Determine the appropriate flags for this compiler for each build type.
# For each option type, a list of possible flags is given that work
# for various compilers.  The first flag that works is chosen.
# If none of the flags work, nothing is added (unless the REQUIRED 
# flag is given in the call).  This way unknown compiles are supported.
#######################################################################

#####################
### GENERAL FLAGS ###
#####################

#message(STATUS "CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
#unset(CMAKE_Fortran_FLAGS)

# Don't add underscores in symbols for C-compatability
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-fno-underscoring")

# Ensure standard semantics
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran REQUIRED   "/standard-semantics"   # Intel Windows
                                    "-standard-semantics"   # Intel
                                    "-std=gnu"              # GNU
                 )

# ensure IFPORT module of Intel ifort is called if compiler is ifort
if( ${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel|INTEL|intel" )
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
            Fortran REQUIRED
            "/DINTEL"                       # Intel
            )
endif()

# add coarray flags if needed
if (NOT CAF STREQUAL "")
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                    Fortran REQUIRED
                    "/Qcoarray=${CAF}"      # Intel Windows
                    "-coarray=${CAF}"       # Intel
                    "-fcoarray=${CAF}"      # GNU
                    )
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                    Fortran REQUIRED
                    "/DCAF"      # Intel Windows
                    "-DCAF"      # Intel/GNU
                    )
    
endif()

# add Fortran preprocessor flags
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran REQUIRED
                "/fpp"  # Intel Windows
                "-fpp"  # Intel
                "-cpp"  # GNU
                )

# add Fortran unlimited line length flags: Intel's default is virtually unlimited
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran
                "-ffree-line-length-none"   # GNU
                )

# add Fortran shared library generation flags
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran
                "/libs:dll /dll"                # Intel Windows
                "-shared"                       # Intel/GNU
                "-dynamiclib -arch_only i386 -noall_load -weak_references_mismatches non-weak"  # Intel Mac
                )

# Intel-provided libraries to be linked in statically
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran
                "-static-intel"                 # Intel Windows
                                                # Intel: none
                                                # GNU: none
                )

# Enforce ifort linker to search for unresolved references
# in a multithreaded run-time library.
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran
                "/threads"                      # Intel Windows
                "-threads"                      # Intel
                                                # GNU: none
                )

# Request compiler to generate position-independent code.
# Option -fpic must be used when building shared objects.
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran
                                                # Intel Windows: none
                "-fpic"                         # Intel Linux/Mac/GNU
                "-fPIC"                         # GNU
                )

# There is some bug where -march=native doesn't work on Mac
IF(APPLE)
    SET(GNUNATIVE "-mtune=native")
ELSE()
    SET(GNUNATIVE "-march=native")
ENDIF()
# Optimize for the host's architecture
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran
                "/QxHost"                       # Intel Windows
                "-xHost"                        # Intel
                ${GNUNATIVE}                    # GNU
                "-ta=host"                      # PGI
                )

# set link time stack size of the executable to ~4GB
#SET_COMPILE_FLAG(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}"
SET_COMPILE_FLAG(CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS}"
                Fortran
                "/F0x100000000"                     # Intel Windows
                                                    # Intel 
                                                    # GNU
                                                    # PGI
                )



###################
### DEBUG FLAGS ###
###################

# NOTE: debugging symbols (-g or /debug:full) are already on by default

#unset(CMAKE_Fortran_FLAGS_DEBUG)

if (BTYPE STREQUAL "debug")

    # generate full debug information
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                    Fortran REQUIRED
                    "/debug:full"               # Intel Windows
                    "-debug full"               # Intel
                    )

    # generate full debug information
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                    Fortran REQUIRED
                    "/Zi"                       # Intel Windows
                    "-g"                        # Intel/GNU
                    )

    # Disable optimizations
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                    Fortran REQUIRED
                    "/Od"                       # Intel Windows
                    "-O0"                       # All compilers not on Windows
                    )

    # Perform run-time checks on whether array subscript and substring references are within declared bounds (same as the -check bounds option)
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                    Fortran
                    "/CB"                       # Intel Windows
                    "-CB"                       # Intel
                    "-fcheck=bounds"            # GNU (New style)
                    "-fbounds-check"            # GNU (Old style)
                    "-Mbounds"                  # PGI
                    )
    
    # initialize variables to zero or to various numeric exceptional values
    # snan: Determines whether the compiler initializes to signaling NaN all uninitialized variables of intrinsic type
    # REAL or COMPLEX that are saved, local, automatic, or allocated variables.
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                    Fortran
                    "/Qinit:snan"               # Intel Windows
                    "-init=snan"                # Intel
                    "-finit-real=snan"          # GNU: no equivalent found
                    )
    
    # Turn on all warnings 
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                    Fortran
                    "/warn:all"                 # Intel Windows
                    "-warn all"                 # Intel
                    "-Wall"                     # GNU
                                                # PGI (on by default)
                    )
    
    # Tell the compiler to generate an interface block for each routine in a source file
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                    Fortran
                    "/gen-interfaces"           # Intel Windows
                    "-gen-interfaces"           # Intel
                                                # GNU: unknown
                                                # PGI: unknown
                    )
    
    # Traceback
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                    Fortran
                    "/traceback"                # Intel Windows
                    "-traceback"                # Intel/PGI
                    "-fbacktrace"               # GNU (gfortran)
                    "-ftrace=full"              # GNU (g95)
                    )
    
    # Check all
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                    Fortran
                    "/check:all"                # Intel Windows
                    "-check all"                # Intel
                    "-fcheck=all"               # GNU
                    )
    
    # Check array bounds
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                    Fortran
                    "/check:bounds"             # Intel Windows
                    "-check bounds"             # Intel
                    "-fcheck=bounds"            # GNU
                    "-fbounds-check"            # G95
                    "-Mbounds"                  # PGI
                    )
    
    # Floating-point invalid, divide-by-zero, and overflow exceptions are enabled
    # throughout the application when the main program is compiled with the value 0
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                    Fortran
                    "/fpe-all:0"                            # Intel Windows
                    "-fpe-all=0"                            # Intel
                    "-ffpe-trap=zero,overflow,underflow"    # GNU
                    )
    
    # Floating-point invalid, divide-by-zero, and overflow exceptions are enabled
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                    Fortran
                    "/Qdiag-error-limit:10"     # Intel Windows
                    "-diag-error-limit=10"      # Intel Linux
                    "-fmax-errors=10"           # GNU
                    )
    
    # Initializes stack local variables to an unusual  value  to  aid error  detection.
    # Normally,  these  local variables should be initialized in the application.
    # The option sets any  uninitialized  local  variables  that  are
    # allocated on the stack to a value that is typically interpreted
    # as a very large integer or an invalid address.   References  to
    # these  variables  are then likely to cause run-time errors that
    # can help you detect coding errors.
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                    Fortran
                    "/Qtrapuv"                  # Intel Windows
                    "-ftrapuv"                  # Intel Linux
                                                # GNU: unknown
                    )
    
    # Enable warning options for usages of language features which may be problematic.
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                    Fortran
                                                # Intel Windows: unknown
                                                # Intel Linux: unknown
                    "-Wextra"                   # GNU
                    )
    
    # Produce a warning when "suspicious" code constructs are encountered. 
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                    Fortran
                                                # Intel Windows: unknown
                                                # Intel Linux: unknown
                    "-Wsurprising"              # GNU
                    )
    
    # Warn about implicit conversions that are likely to change the value of the expression after conversion. Implied by -Wall. 
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                    Fortran
                                                # Intel Windows: unknown
                                                # Intel Linux: unknown
                    "-Wconversion"              # GNU
                    )
    
    # Warn about possible aliasing of dummy arguments. Specifically,
    # it warns if the same actual argument is associated with a dummy argument
    # with INTENT(IN) and a dummy argument with INTENT(OUT) in a call with an explicit interface. 
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                    Fortran
                                                # Intel Windows: unknown
                                                # Intel Linux: unknown
                    "-Waliasing"                # GNU
                    )
    
    # Issue warnings for uses of extensions to Fortran
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                    Fortran
                                                # Intel Windows: unknown
                                                # Intel Linux: unknown
                    "-pedantic"                 # GNU
                    )

endif()

######################
#### TESTING FLAGS ###
######################

if(BTYPE STREQUAL "Testing")

    # Optimizations
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_TESTING "${CMAKE_Fortran_FLAGS_TESTING}"
                    Fortran REQUIRED
                    "-O0"                       # All compilers not on Windows
                    "/O0"                       # Intel Windows
                    )

endif()

#####################
### RELEASE FLAGS ###
#####################

# NOTE: agressive optimizations (-O3) are already turned on by default

#unset(CMAKE_Fortran_FLAGS_RELEASE)
#message( STATUS "CMAKE_Fortran_FLAGS_RELEASE: " ${CMAKE_Fortran_FLAGS_RELEASE} )
#set(CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS} FORCE)
#message( STATUS "CMAKE_Fortran_FLAGS_RELEASE: " ${CMAKE_Fortran_FLAGS_RELEASE} )

if (BTYPE STREQUAL "release")

    # enable optimizations
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                    Fortran REQUIRED
                    "/O3"                           # Intel Windows
                    "-O3"                           # Intel/GNU
                    )

    # Unroll loops: The compiler uses default heuristics when unrolling loops (Intel)
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                    Fortran
                    "/Qunroll"                      # Intel Windows
                    "-unroll"                       # Intel
                    "-funroll-loops"                # GNU
                    "-Munroll"                      # PGI
                    )

    # Unroll loops aggressively
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                    Fortran
                    "/Qunroll-aggressive"           # Intel Windows
                    "-unroll-aggressive"            # Intel
                                                    # GNU: unknown
                    )

    # Inline functions
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                    Fortran
                    "/inline:all /Ob2"                      # Intel Windows
                    "-finline-functions -inline-level=2"    # Intel
                    "-finline-functions"                    # Intel/GNU
                    "-Minline"                              # PGI
                    )
    
    # Interprocedural (link-time) optimizations
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                    Fortran
                    "/Qipo"                         # Intel Windows
                    "-ipo"                          # Intel
                    "-flto"                         # GNU
                    "-Mipa"                         # PGI
                    )
    
    # Single-file optimizations
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                    Fortran
                    "/Qip"                          # Intel Windows
                    "-ip"                           # Intel
                    )

    # This flag is added by default at O3 level and not needed.
    ## Vectorize code: Vectorization is enabled if option O2 or higher is in effect (Intel).
    #SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
    #                Fortran
    #                "/Qvec"                         # Intel Windows
    #                "-vec"                          # Intel
    #                "-ftree-vectorize"              # GNU
    #                "-Mvect"                        # PGI
    #                )
    
    # optimization report: (Intel: 2 is medium level of report details)
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                    Fortran
                    "/Qopt-report:2"                    # Intel Windows
                    "-qopt-report=2"                    # Intel
                    "-fopt-info-all=GFortranOptReport.txt"  # GNU
                    )

endif()
