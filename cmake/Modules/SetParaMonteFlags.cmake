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

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# set path to the Python interpreter
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set( Python_PATH "" CACHE STRING "Path to the Python interpreter" )

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# set up MATLAB
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# unset(MATLAB_FPP_FLAGS)
# unset(MATLAB_LINKER_FLAGS)
# #unset(MATLAB_EXPORT_VERSION)
# unset(MATLAB_PLATFORM_SUBDIR)
#
# if (WIN32)
# else()
#     if (DEFINED MATLAB_ROOT_DIR)
#         if(APPLE)
#             set(MATLAB_PLATFORM_SUBDIR maci64)
#         else()
#             set(MATLAB_PLATFORM_SUBDIR glnxa64)
#         endif()
#         #if (intel_compiler)
#         #    set(MATLAB_FPP_FLAGS -DMATLAB_MEX_FILE)
#         #elseif(gnu_compiler)
#             set(MATLAB_FPP_FLAGS -DMATLAB_MEX_FILE -D_GNU_SOURCE -DMEXPRINT_ENABLED)
#         #endif()
#         set(MATLAB_LINKER_FLAGS
#         -Wl,--no-undefined
#         -Wl,--as-needed -Wl,-rpath-link,"${MATLAB_ROOT_DIR}/bin/${MATLAB_PLATFORM_SUBDIR}"
#         #-L"$(MATLAB_ROOT_DIR)/bin/${MATLAB_PLATFORM_SUBDIR}"
#         -Wl,-rpath-link,"$(MATLAB_ROOT_DIR)/extern/bin/${MATLAB_PLATFORM_SUBDIR}"
#         -Wl,--version-script,"$(MATLAB_ROOT_DIR)/extern/lib/${MATLAB_PLATFORM_SUBDIR}/fortran_exportsmexfileversion.map"
#         )
#         #set(MATLAB_EXPORT_VERSION -Wl,--version-script,"$(MATLAB_ROOT_DIR)/extern/lib/${MATLAB_PLATFORM_SUBDIR}/fortran_exportsmexfileversion.map")
#     endif()
# endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# unset flags
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

message(STATUS "${pmattn} Default CMAKE_C_FLAGS:       ${CMAKE_C_FLAGS}")
message(STATUS "${pmattn} Default CMAKE_CXX_FLAGS:     ${CMAKE_CXX_FLAGS}")
message(STATUS "${pmattn} Default CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")

unset(FL_FLAGS)
unset(FPP_FLAGS)
unset(CCL_FLAGS)
unset(FCL_FLAGS)
unset(CCL_FLAGS_GNU)
unset(FCL_FLAGS_INTEL)
unset(CMAKE_Fortran_FLAGS)

set(FPP_FLAGS "${FPP_FLAGS}" ${USER_PREPROCESSOR_MACROS})

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# set preprocesssing compiler flags
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# to save the intermediate files via ifort: FPP /Qsave_temps <original file> <intermediate file>
# Will be used to pass any user-defined macro definition to the intel compiler for preprocessing the files.
# A macro is denoted by /define:MACRO_NAME=VALUE, the value of which can be dropped, and if so, it will be given a default value of 1.
# Example: /define:INTEL, will create a macro INTEL with a default value of 1.
# set( INTEL_Fortran_PREPROCESSOR_MACROS "" CACHE STRING "Intel Fortran compiler preprocessor definitions" )

if (intel_compiler)
    if (WIN32)
        set(FPP_FLAGS "${FPP_FLAGS}"
        /fpp /DINTEL_COMPILER_ENABLED
        )
    else()
        set(FPP_FLAGS "${FPP_FLAGS}"
        -fpp -DINTEL_COMPILER_ENABLED
        )
    endif()
elseif (gnu_compiler)
    set(FPP_FLAGS "${FPP_FLAGS}"
    -cpp -DGNU_COMPILER_ENABLED
    )
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# set C-Fortran preprocessor interoperability flag
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (CFI_ENABLED)
    set(FPP_FLAGS "${FPP_FLAGS}"
    -DCFI_ENABLED
    )
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ParaMonte Version Preprocessor Flag: This flag is not used anymore as it is too aggressive. Changing its value causes an entire rebuild.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (DEFINED fppParaMonteVersion)
    set(FPP_FLAGS "${FPP_FLAGS}"
    -DPARAMONTE_VERSION=\"'${fppParaMonteVersion}'\"
    )
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ParaMonte OS Preprocessor Flag
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if(WIN32)
    set(FPP_FLAGS "${FPP_FLAGS}"
    -DOS_IS_WINDOWS
    )
elseif(APPLE)
    set(FPP_FLAGS "${FPP_FLAGS}"
    -DOS_IS_DARWIN
    )
elseif(UNIX AND NOT APPLE)
    set(FPP_FLAGS "${FPP_FLAGS}"
    -DOS_IS_LINUX
    )
endif()

if(DEFINED OS_IS_WSL)
    set(FPP_FLAGS "${FPP_FLAGS}"
    -DOS_IS_WSL
    )
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# set language preprocessor flag
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (${INTERFACE_LANGUAGE} STREQUAL "c" OR ${INTERFACE_LANGUAGE} STREQUAL "C")
    set(FPP_FLAGS "${FPP_FLAGS}" -DC_ENABLED)
elseif (${INTERFACE_LANGUAGE} STREQUAL "c++" OR ${INTERFACE_LANGUAGE} STREQUAL "C++")
    set(FPP_FLAGS "${FPP_FLAGS}" -DCPP_ENABLED)
elseif (${INTERFACE_LANGUAGE} MATCHES "[fF][oO][rR][tT][rR][aA][nN]")
    set(FPP_FLAGS "${FPP_FLAGS}" -DFORTRAN_ENABLED)
elseif (${INTERFACE_LANGUAGE} MATCHES "[jJ][aA][vV][aA]")
    set(FPP_FLAGS "${FPP_FLAGS}" -DJAVA_ENABLED)
elseif (${INTERFACE_LANGUAGE} MATCHES "[jJ][uU][lL][iI][aA]")
    set(FPP_FLAGS "${FPP_FLAGS}" -DJULIA_ENABLED)
elseif (${INTERFACE_LANGUAGE} MATCHES "[mM][aA][tT][hH][eE][mM][aA][tT][iI][cC][aA]")
    set(FPP_FLAGS "${FPP_FLAGS}" -DMATHEMATICA_ENABLED)
elseif (${INTERFACE_LANGUAGE} MATCHES "[mM][aA][tT][lL][aA][bB]")
    set(FPP_FLAGS "${FPP_FLAGS}" -DMATLAB_ENABLED)
elseif (${INTERFACE_LANGUAGE} MATCHES "[pP][yY][tT][hH][oO][nN]")
    set(FPP_FLAGS "${FPP_FLAGS}" -DPYTHON_ENABLED)
elseif (${INTERFACE_LANGUAGE} MATCHES "[rR]")
    set(FPP_FLAGS "${FPP_FLAGS}" -DR_ENABLED)
else()
    message (FATAL_ERROR
            "\n"
            "${pmfatal}\n"
            "${pmfatal} Unrecognized interface language detected.\n"
            "${pmfatal} INTERFACE_LANGUAGE: ${INTERFACE_LANGUAGE}\n"
            "${pmfatal} possible values are: C/Fortran/MATLAB/Python.\n"
            "\n"
            )
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: set preprocessor build flags
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (CMAKE_BUILD_TYPE MATCHES "Debug|DEBUG|debug")
    set(FPP_FLAGS "${FPP_FLAGS}" -DDEBUG_ENABLED)
elseif (CMAKE_BUILD_TYPE MATCHES "Testing|TESTING|testing")
    set(FPP_FLAGS "${FPP_FLAGS}" -DTESTING_ENABLD)
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: set the testing type
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (BASIC_TEST_ENABLED)
    set(FPP_FLAGS "${FPP_FLAGS}" -DBASIC_TEST_ENABLED)
endif()

if (SAMPLER_TEST_ENABLED)
    set(FPP_FLAGS "${FPP_FLAGS}" -DSAMPLER_TEST_ENABLED)
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: set code coverage preprocessor flag
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (CODECOV_ENABLED)
    set(FPP_FLAGS "${FPP_FLAGS}" -DCODECOV_ENABLED)
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: set shared library preprocessor flag
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (pmlib_shared)
    set(FPP_FLAGS "${FPP_FLAGS}" -DDLL_ENABLED)
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: set set default Fortran compiler flags in different build modes.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

option(MT_ENABLED "Enabled multithreaded library linking" ON)

if (intel_compiler)

    if (WIN32)
        set(FCL_FLAGS "${FCL_FLAGS}"
        /standard-semantics   # determines whether the current Fortran Standard behavior of the compiler is fully implemented.
        /F0x1000000000        # specify the stack reserve amount for the program.
        /nologo               # no logo
        )
        if (MT_ENABLED)
            set(FCL_FLAGS "${FCL_FLAGS}"
            /threads
            )
        endif()
    else()
        set(FCL_FLAGS "${FCL_FLAGS}"
        -standard-semantics   # determines whether the current Fortran Standard behavior of the compiler is fully implemented.
        -nologo               # no logo
        )
        if (MT_ENABLED)
            set(FCL_FLAGS "${FCL_FLAGS}"
            -threads
            )
        endif()
    endif()

    if (CODECOV_ENABLED)
        set(FCL_FLAGS "${FCL_FLAGS}" -prof-gen=srcpos)
        set(CCL_FLAGS "${CCL_FLAGS}" -prof-gen=srcpos)
        set(FL_FLAGS "${FL_FLAGS}" -prof-gen=srcpos)
    endif()

elseif (gnu_compiler)

    # -std=legacy is required to bypass the new gfortran 10 error on argument-mismatch.
    # The MPICH 3.2 library mpi_bcast still requires argument-mismatch, which causes gfortran to break the compilation.
    # The problem still persists in debug mode. Therefore, when gfortran is 10, debug mode is disabled.
    #set(FCL_FLAGS -std=gnu -ffree-line-length-none -fallow-argument-mismatch CACHE STRING "GNU Fortran default compiler flags" )

    set(FCL_FLAGS "${FCL_FLAGS}"
    -ffree-line-length-none
    -std=legacy
    )

    set(CCL_FLAGS "${CCL_FLAGS}"
    -ffree-line-length-none
    )

    if (MT_ENABLED)
        set(FCL_FLAGS "${FCL_FLAGS}" -pthread)
        set(CCL_FLAGS "${CCL_FLAGS}" -pthread)
    endif()

    if (CODECOV_ENABLED)
        set(FCL_FLAGS "${FCL_FLAGS}"
        # -fcf-protection=full
        -ftest-coverage
        -fprofile-arcs
        --coverage
        )
        set(CCL_FLAGS "${CCL_FLAGS}"
        # -fcf-protection=full
        -ftest-coverage
        -fprofile-arcs 
        --coverage 
        )
        set(FL_FLAGS "${FL_FLAGS}"
        -fcf-protection=full
        -ftest-coverage 
        -fprofile-arcs 
        -static-libgcc 
        --coverage 
        -lgcov 
        )
    endif()

endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# set C/Fortran compiler/linker build flags
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (intel_compiler)

    if (${CMAKE_BUILD_TYPE} MATCHES "[Dd][eE][bB][uU][gG]")

        if (WIN32)

            set(CCL_FLAGS_INTEL "${CCL_FLAGS_INTEL}"
            /debug:full
            /Zi
            /Od
            /Wall
            /traceback
            /Qcheck-pointers:rw                 # check bounds for memory access through pointers.
            /Qcheck-pointers-undimensioned      # check bounds for memory access through arrays that are declared without dimensions.
            /Qdiag-error-limit:10
            /Qtrapuv
            )

            set(FCL_FLAGS_INTEL "${FCL_FLAGS_INTEL}"
            /debug:full
            /stand:f08                          # issue compile-time messages for nonstandard language elements.
            /Zi
            /CB
            /Od
            /Qinit:snan,arrays
            /warn:all
            /gen-interfaces
            /traceback
            /check:all
            /check:bounds
            /fpe-all:0
            /Qdiag-error-limit:10
            /Qdiag-disable:5268
            /Qdiag-disable:7025
            /Qtrapuv
            )

        else()

            set(CCL_FLAGS_INTEL "${CCL_FLAGS_INTEL}"
            -debug full
            -g3
            -O0
            -Wall
            -traceback
            -check-pointers=rw                    # check bounds for memory access through pointers.
            -check-pointers-undimensioned         # check bounds for memory access through arrays that are declared without dimensions.
            -diag-error-limit=10
            -ftrapuv
            )

            set(FCL_FLAGS_INTEL "${FCL_FLAGS_INTEL}"
            -stand f08                # issue compile-time messages for nonstandard language elements.
            -debug full               # generate full debug information
            -g3                       # generate full debug information
            -O0                       # disable optimizations
            -CB                       # Perform run-time bound-checks on array subscript and substring references (same as the -check bounds option)
            -init:snan,arrays         # initialize arrays and scalars to NaN
            -warn all                 # enable all warning
            -gen-interfaces           # generate interface block for each routine in the source files
            -traceback                # trace back for debugging
            -check all                # check all
            -check bounds             # check array bounds
            -fpe-all=0                # Floating-point invalid, divide-by-zero, and overflow exceptions are enabled
            -diag-error-limit=10      # max diagnostic error count
            -diag-disable=5268        # Extension to standard: The text exceeds right hand column allowed on the line.
            -diag-disable=7025        # This directive is not standard Fxx.
            -diag-disable=10346       # optimization reporting will be enabled at link time when performing interprocedural optimizations.
            -ftrapuv                  # Initializes stack local variables to an unusual value to aid error detection.
            )

        endif()

    elseif (${CMAKE_BUILD_TYPE} MATCHES "[tT][eE][sS][tT][iI][nN][gG]")

        if (WIN32)

            set(CCL_FLAGS_INTEL "${CCL_FLAGS_INTEL}"
            /Od  # disable optimization
            )

            set(FCL_FLAGS_INTEL "${FCL_FLAGS_INTEL}"
            /Od  # disable optimization
            )

        else()

            set(CCL_FLAGS_INTEL "${CCL_FLAGS_INTEL}"
            -O0  # disable optimization
            )

            set(FCL_FLAGS_INTEL "${FCL_FLAGS_INTEL}"
            -O0  # disable optimization
            )

        endif()

    elseif (${CMAKE_BUILD_TYPE} MATCHES "[rR][eE][lL][eE][aA][sS][eE]")

        if(WIN32)

            set(CCL_FLAGS_INTEL "${CCL_FLAGS_INTEL}"
            /O3
            /Qip
            /Qipo
            /Qunroll
            /Qunroll-aggressive
            /Ob2
            /Qinline-dllimport
            # /Qparallel  # Tells the auto-parallelizer to generate multithreaded code for loops that can be safely executed in parallel.
            )
            # /Qipo-c generates a single object file, containing all object files.
            # INTEL_Fortran_RELEASE_FLAGS=/fast /O3 /Qip /Qipo /Qunroll /Qunroll-aggressive /inline:all /Ob2 /Qparallel /Qinline-dllimport

            set(FCL_FLAGS_INTEL "${FCL_FLAGS_INTEL}"
            /O3                         # set the optimizations level
            /Qip                        # determines whether additional interprocedural optimizations for single-file compilation are enabled.
            /Qipo                       # enable interprocedural optimization between files.
            /Qunroll                    # [:n] set the maximum number of times to unroll loops (no number n means automatic).
            /Qunroll-aggressive         # use more aggressive unrolling for certain loops.
            /Qvec                       # enable vectorization.
            /Qopt-report:2              # generate optimization report. Level 2 is the default. Use 5 for the greatest details and 0 for no report.
           #/Qguide-vec:4               # enable guidance for auto-vectorization, causing the compiler to generate messages suggesting ways to improve optimization (default=4, highest).
           #/Qparallel                  # generate multithreaded code for loops that can be safely executed in parallel.
           #/Qipo-c:                    # Tells the compiler to optimize across multiple files and generate a single object file ipo_out.obj without linking
                                        # info at: https://software.intel.com/en-us/Fortran-compiler-developer-guide-and-reference-ipo-c-qipo-c
            )

        else()

            set(CCL_FLAGS_INTEL "${CCL_FLAGS_INTEL}"
            -O3
            -ip
            -ipo
            -unroll
            -unroll-aggressive
            -inline-level=2
            # -parallel   # Tells the auto-parallelizer to generate multithreaded code for loops that can be safely executed in parallel.
            )

            set(FCL_FLAGS_INTEL "${FCL_FLAGS_INTEL}"
            -O3                         # set the optimizations level
            -ip                         # determines whether additional interprocedural optimizations for single-file compilation are enabled.
           #-ipo                        # enable interprocedural optimization between files.
            -unroll                     # [=n] set the maximum number of times to unroll loops (no number n means automatic).
            -unroll-aggressive          # use more aggressive unrolling for certain loops.
            -finline-functions          # enables function inlining for single file compilation.
            -diag-disable=10346         # optimization reporting will be enabled at link time when performing interprocedural optimizations.
            -diag-disable=10397         # optimization reporting will be enabled at link time when performing interprocedural optimizations.
            -qopt-report=2              # generate optimization report. Level 2 is the default. Use 5 for the greatest details and 0 for no report.
           #-guide-vec=4                # enable guidance for auto-vectorization, causing the compiler to generate messages suggesting ways to improve optimization (default=4, highest).
           #-parallel                   # generate multithreaded code for loops that can be safely executed in parallel. This option requires MKL libraries.
           #-qopt-subscript-in-range    # assumes there are no "large" integers being used or being computed inside loops. A "large" integer is typically > 2^31.
            )

        endif()

    else()

        message (FATAL_ERROR
                "\n"
                "${pmfatal}\n"
                "${pmfatal} Unrecognized built type: ${CMAKE_BUILD_TYPE}.\n"
                "${pmfatal} possible values are: debug/testing/release.\n"
                "\n"
                )

    endif() # build type

elseif(gnu_compiler)

    if (${CMAKE_BUILD_TYPE} MATCHES "[Dd][eE][bB][uU][gG]")

        set(CCL_FLAGS_GNU "${CCL_FLAGS_GNU}"
        -g                                    # generate full debug information
        -O0                                   # disable optimizations
        -fcheck=all                           # enable the generation of run-time checks
        -ffpe-trap=zero,overflow,underflow    # Floating-point invalid, divide-by-zero, and overflow exceptions are enabled
        -finit-real=snan                      # initialize REAL and COMPLEX variables with a signaling NaN
        -fbacktrace                           # trace back for debugging
        --pedantic                            # issue warnings for uses of extensions to the Fortran standard
        -fmax-errors=10                       # max diagnostic error count
        -Wall                                 # enable all warnings:
                                              # -Waliasing, -Wampersand, -Wconversion, -Wsurprising, -Wc-binding-type, -Wintrinsics-std, -Wtabs, -Wintrinsic-shadow,
                                              # -Wline-truncation, -Wtarget-lifetime, -Winteger-division, -Wreal-q-constant, -Wunused, -Wundefined-do-loop
        )

        set(FCL_FLAGS_GNU "${FCL_FLAGS_GNU}"
        -g3                                 # generate full debug information
        -O0                                 # disable optimizations
       #-fsanitize=undefined                # enable UndefinedBehaviorSanitizer for undefined behavior detection.
       #-fsanitize=address                  # enable AddressSanitizer, for memory error detection, like out-of-bounds and use-after-free bugs.
       #-fsanitize=leak                     # enable LeakSanitizer for memory leak detection.
        -fcheck=all                         # enable the generation of run-time checks
        -ffpe-trap=zero,overflow,underflow  # Floating-point invalid, divide-by-zero, and overflow exceptions are enabled
        -finit-real=snan                    # initialize REAL and COMPLEX variables with a signaling NaN
        -fbacktrace                         # trace back for debugging
       #--pedantic                          # issue warnings for uses of extensions to the Fortran standard. Gfortran10 with MPICH 3.2 in debug mode crashes with this flag at mpi_bcast. Excluded until MPICH upgraded.
        -fmax-errors=10                     # max diagnostic error count
        -Wno-maybe-uninitialized            # avoid warning of no array pre-allocation.
        -Wall                               # enable all warnings:
                                            # -Waliasing, -Wampersand, -Wconversion, -Wsurprising, -Wc-binding-type, -Wintrinsics-std, -Wtabs, -Wintrinsic-shadow,
                                            # -Wline-truncation, -Wtarget-lifetime, -Winteger-division, -Wreal-q-constant, -Wunused, -Wundefined-do-loop
                                            # gfortran10 crashes and cannot compile MPI ParaMonte with mpich in debug mode. Therefore -wall is disabled for now, until MPICH upgrades interface.
        #-Waliasing
        #-Wampersand
        #-Wconversion
        #-Wsurprising
        #-Wc-binding-type
        #-Wintrinsics-std
        #-Wtabs
        #-Wintrinsic-shadow
        #-Wline-truncation
        #-Wtarget-lifetime
        #-Winteger-division
        #-Wreal-q-constant
        #-Wunused
        #-Wundefined-do-loop
        )

    elseif (${CMAKE_BUILD_TYPE} MATCHES "[tT][eE][sS][tT][iI][nN][gG]") # TESTING

        set(CCL_FLAGS_GNU "${CCL_FLAGS_GNU}"
        -O0  # disable optimization
        )

        set(FCL_FLAGS_GNU "${FCL_FLAGS_GNU}"
        #-static-libgfortran -static-libgcc
        -O0  # disable optimization
        )

    elseif (${CMAKE_BUILD_TYPE} MATCHES "[rR][eE][lL][eE][aA][sS][eE]") # RELEASE

        set(CCL_FLAGS_GNU "${CCL_FLAGS_GNU}"
        -O3                       # set the optimizations level
        -flto                     # enable interprocedural optimization between files.
        -funroll-loops            # [=n] set the maximum number of times to unroll loops (no number n means automatic).
        -finline-functions        # consider all functions for inlining, even if they are not declared inline.
        -ftree-vectorize          # perform vectorization on trees. enables -ftree-loop-vectorize and -ftree-slp-vectorize.
        )

        set(FCL_FLAGS_GNU "${FCL_FLAGS_GNU}" -fopt-info-all=GFortranOptReport.txt)
        if(APPLE)
            # GNU 9.3 -O results in runtime crashes of the examples on Mac (except -O0)
            set(FCL_FLAGS_GNU "${FCL_FLAGS_GNU}"
            -fauto-inc-dec
            -fbranch-count-reg
            -fcombine-stack-adjustments
            -fcompare-elim
            -fcprop-registers
            -fdce
            -fdefer-pop
            #-fdelayed-branch
            -fdse
            -fforward-propagate
            -fguess-branch-probability
            -fif-conversion
            -fif-conversion2
            -finline-functions-called-once
            -fipa-profile
            -fipa-pure-const
            -fipa-reference
            -fipa-reference-addressable
            -fmerge-constants
            -fmove-loop-invariants
            -fomit-frame-pointer
            -freorder-blocks
            -fshrink-wrap
            -fshrink-wrap-separate
            -fsplit-wide-types
            -fssa-backprop
            -fssa-phiopt
            -ftree-bit-ccp
            -ftree-ccp
            -ftree-ch
            -ftree-coalesce-vars
            -ftree-copy-prop
            -ftree-dce
            -ftree-dominator-opts
            -ftree-dse
            -ftree-forwprop
            -ftree-fre
            -ftree-phiprop
            -ftree-pta
            -ftree-scev-cprop
            -ftree-sink
            -ftree-slsr
            -ftree-sra
            -ftree-ter
            -funit-at-a-time
            -falign-functions  -falign-jumps
            -falign-labels  -falign-loops
            -fcaller-saves
            -fcode-hoisting
            -fcrossjumping
            -fcse-follow-jumps  -fcse-skip-blocks
            -fdelete-null-pointer-checks
            -fdevirtualize  -fdevirtualize-speculatively
            -fexpensive-optimizations
            -fgcse  -fgcse-lm
            -fhoist-adjacent-loads
            -finline-functions
            -finline-small-functions
            -findirect-inlining
            -fipa-bit-cp  -fipa-cp  -fipa-icf
            -fipa-ra  -fipa-sra  -fipa-vrp
            -fisolate-erroneous-paths-dereference
            -flra-remat
            -foptimize-sibling-calls
            -foptimize-strlen
            -fpartial-inlining
            -fpeephole2
            -freorder-blocks-algorithm=stc
            -freorder-blocks-and-partition  -freorder-functions
            -frerun-cse-after-loop
            -fschedule-insns  -fschedule-insns2
            -fsched-interblock  -fsched-spec
            -fstore-merging
            -fstrict-aliasing
            -fthread-jumps
            -ftree-builtin-call-dce
            -ftree-pre
            -ftree-switch-conversion  -ftree-tail-merge
            -ftree-vrp
            -fgcse-after-reload
            -fipa-cp-clone
            -floop-interchange
            -floop-unroll-and-jam
            -fpeel-loops
            -fpredictive-commoning
            -fsplit-paths
            -ftree-loop-distribute-patterns
            -ftree-loop-distribution
            -ftree-loop-vectorize
            -ftree-partial-pre
            -ftree-slp-vectorize
            -funswitch-loops
            -fvect-cost-model
            -fversion-loops-for-strides
            )

        else()

            set(FCL_FLAGS_GNU "${FCL_FLAGS_GNU}"
            #-static-libgfortran -static-libgcc
            -O3                       # set the optimizations level
            #-flto                    # enable interprocedural optimization between files.
            -funroll-loops            # [=n] set the maximum number of times to unroll loops (no number n means automatic).
            -finline-functions        # consider all functions for inlining, even if they are not declared inline.
            -ftree-vectorize          # perform vectorization on trees. enables -ftree-loop-vectorize and -ftree-slp-vectorize.
            )

        endif()

    else()

        message (FATAL_ERROR
                "\n"
                "${pmfatal}\n"
                "${pmfatal} Unrecognized built type: ${CMAKE_BUILD_TYPE}.\n"
                "${pmfatal} possible values are: debug/testing/release.\n"
                "\n"
                )

    endif() # build type

endif() # compiler suite

if (gnu_compiler)
    set(FCL_FLAGS "${FCL_FLAGS}" "${FCL_FLAGS_GNU}")
    set(CCL_FLAGS "${CCL_FLAGS}" "${CCL_FLAGS_GNU}")
elseif (intel_compiler)
    set(FCL_FLAGS "${FCL_FLAGS}" "${FCL_FLAGS_INTEL}")
    set(CCL_FLAGS "${CCL_FLAGS}" "${CCL_FLAGS_INTEL}")
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# set shared library Fortran linker flags
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (pmlib_shared)

    if (gnu_compiler)

        set(FCL_FLAGS "${FCL_FLAGS}"
        -shared
        -fPIC
        )

    else(intel_compiler)

        if(WIN32)
            set(FCL_FLAGS "${FCL_FLAGS}"
            #/threads " # these flags are actually included by default in recent ifort implementations
            /libs:dll
            )
        elseif(UNIX)
            set(FCL_FLAGS "${FCL_FLAGS}"
            -fpic # Request compiler to generate position-independent code.
            )
            if (APPLE)
                set(FCL_FLAGS "${FCL_FLAGS}"
                -noall_load
                # -weak_references_mismatches non-weak -threads -arch_only i386
                )
            endif()
        endif()

    endif()

endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: set up coarray flags
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

message( STATUS "${pmattn} setting up Coarray Fortran (CAF) parallelization model. Options: single, shared, distributed" )
message( STATUS "${pmattn} requested CAF: ${CAFTYPE}" )

if( "${CAFTYPE}" STREQUAL "single" OR
    "${CAFTYPE}" STREQUAL "shared" OR
    "${CAFTYPE}" STREQUAL "distributed" )
    set(FPP_FLAGS "${FPP_FLAGS}" -DCAF_ENABLED)
    message( STATUS "${pmattn} enabling Coarray Fortran syntax via preprocesor flag -DCAF_ENABLED" )
    set(CAF_ENABLED ON CACHE BOOL "Enable Coarray Fortran parallelism" FORCE)
    set(FOR_COARRAY_NUM_IMAGES 2 CACHE STRING "The default number of coarray images/processes" FORCE)
    if (intel_compiler)
        set(FCL_FLAGS "${FCL_FLAGS}" 
        -coarray=${CAFTYPE}
        )
    endif()
else()
    message( STATUS "${pmattn} ignoring Coarray Fortran parallelization." )
    set(CAF_ENABLED OFF CACHE BOOL "Coarray Fortran parallelism" FORCE)
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: set non-coarray parallelization flags and definitions to be passed to the preprocessors
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (MPI_ENABLED)
    if (CAF_ENABLED)
        message (FATAL_ERROR
                "\n"
                "${pmfatal}\n"
                "${pmfatal} Coarray parallelism cannot be currently mixed with with MPI.\n"
                "${pmfatal} MPI_ENABLED: ${MPI_ENABLED}\n"
                "${pmfatal} CAF_ENABLED: ${CAF_ENABLED}\n"
                "${pmfatal} CAFTYPE: ${CAFTYPE}\n"
                "${pmfatal} set MPI_ENABLED and CAFTYPE to appropriate values in the ParaMonte CMAKE CACHE file and rebuild.\n"
                )
    else()
        set(FPP_FLAGS "${FPP_FLAGS}" -DMPI_ENABLED)
    endif()
endif()

if (OMP_ENABLED)
    set(FPP_FLAGS "${FPP_FLAGS}" -DOMP_ENABLED)
    if (intel_compiler)
        if(WIN32)
            set(FCL_FLAGS "${FCL_FLAGS}"
            /Qopenmp
            )
        else()
            set(FCL_FLAGS "${FCL_FLAGS}"
            -qopenmp
            )
        endif()
    elseif (gnu_compiler)
        set(FCL_FLAGS "${FCL_FLAGS}"
        -fopenmp
        )
    endif()
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: set memory allocation type
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (HEAP_ARRAY_ENABLED)
    if (intel_compiler)
        if (WIN32)
            set(FCL_FLAGS "${FCL_FLAGS}"
            /heap-arrays:10
            )
        else()
            set(FCL_FLAGS "${FCL_FLAGS}"
            -heap-arrays=10
            )
        endif()
    elseif(gnu_compiler)
        set(FCL_FLAGS "${FCL_FLAGS}"
        -fmax-stack-var-size=10
        # -frecursive overwrites -fmax-stack-var-size=10 and causes all allocations to happen on stack.
        )
        # @todo:
        # Strangely, -frecursive flag causes deadlock with sampler tests in Coarray mode.
        # Therefore, the -frecursive flag is reversed to -fmax-stack-var-size=10 until the behavior is understood.
        # The use of -frecursive was based on the GFortran-10 warning message about unsafe storage move from stack to static storage.
        # See also this thread: https://comp.lang.fortran.narkive.com/WApl1KMt/gfortran-stack-size-warning#post4
        # Perhaps adding `non_recursive` to functions would fix this warning message.
    endif()
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: report build spec
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

message( STATUS "${pmattn} operating system / platform: ${CMAKE_SYSTEM_NAME} / ${CMAKE_SYSTEM_PROCESSOR}" )
message( STATUS "${pmattn} C-Fortran compiler suite / version: ${CMAKE_Fortran_COMPILER_ID} / ${CMAKE_Fortran_COMPILER_VERSION}" )
message( STATUS "${pmattn} C-Fortran interoperation enabled - CFI_ENABLED: ${CFI_ENABLED}" )
message( STATUS "${pmattn} Coarray parallelization enabled - CAF_ENABLED: ${CAF_ENABLED}" )
message( STATUS "${pmattn} MPI parallelization enabled - MPI_ENABLED: ${MPI_ENABLED}" )
message( STATUS "${pmattn} library test enabled - TEST_RUN_ENABLED: ${ParaMonteTest_RUN_ENABLED}" )
message( STATUS "${pmattn} library build type: ${CMAKE_BUILD_TYPE}" )
message( STATUS "${pmattn} library type: ${LTYPE}" )
message( STATUS "${pmattn} pmlib_shared: ${pmlib_shared}" )
message( STATUS "${pmattn} C/C++ compiler/linker flags: ${CCL_FLAGS}" )
message( STATUS "${pmattn} Fortran compiler/linker flags - FCL_FLAGS: ${FCL_FLAGS}" )
message( STATUS "${pmattn} preprocessor flags - FPP_FLAGS: ${FPP_FLAGS}" )

#set(CMAKE_Fortran_FLAGS "${FPP_FLAGS}" CACHE STRING "Fortran compiler/linker flags" FORCE)
#string(REPLACE ";" " " CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
