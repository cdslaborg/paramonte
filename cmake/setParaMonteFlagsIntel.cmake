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

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set preprocesssing compiler flags.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# to save the intermediate files via ifort: FPP /Qsave_temps <original file> <intermediate file>
# Will be used to pass any user-defined macro definition to the intel compiler for preprocessing the files.
# A macro is denoted by /define:MACRO_NAME=VALUE, the value of which can be dropped, and if so, it will be given a default value of 1.
# Example: /define:INTEL, will create a macro INTEL with a default value of 1.
# set( INTEL_Fortran_PREPROCESSOR_MACROS "" CACHE STRING "Intel Fortran compiler preprocessor definitions" )

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Set the preprocessed source file generation flag.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (FPP_ENABLED)
    if (WIN32)
        set(FC_FLAGS "${FC_FLAGS}" /Qsave-temps) # enable preprocessed source file generation.
    else()
        set(FC_FLAGS "${FC_FLAGS}" -save-temps) # enable preprocessed source file generation.
    endif()
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Set the default Fortran compiler flags for all build modes.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (WIN32)
    set(FCL_FLAGS "${FCL_FLAGS}"
    /nologo                 # no logo.
    /F0x1000000000          # specify the stack reserve amount for the program.
    /integer-size:32        # Use 32-bits representation as default integer/logical kind.
    /standard-semantics     # determines whether the current Fortran Standard behavior of the compiler is fully implemented.
    /Qdiag-disable=5268     # Extension to standard: The text exceeds right hand column allowed on the line.
    /Qdiag-disable=7025     # This directive is not standard Fxx.
    /Qdiag-disable=10346    # optimization reporting will be enabled at link time when performing interprocedural optimizations.
    /Qdiag-disable:10448    # Disable ifort deprecation message.
    )
    if (MT_ENABLED)
        set(FCL_FLAGS "${FCL_FLAGS}"
        /threads
        )
    endif()
    #if (${csid_is_intel} AND ${MPI_ENABLED})
    #    #set(FCL_FLAGS "${MPI_LINK_FLAGS}" "-static_mpi")
    #    set(MPI_LINK_FLAGS "${MPI_LINK_FLAGS}" "-static_mpi")
    #    set(MPI_COMPILER_FLAGS "${MPI_COMPILER_FLAGS}" "-static_mpi")
    #    set(MPI_Fortran_LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}" "-static_mpi")
    #    set(MPI_Fortran_COMPILE_OPTIONS "${MPI_Fortran_COMPILE_OPTIONS}" "-static_mpi")
    #endif()
else()
    set(FCL_FLAGS "${FCL_FLAGS}"
    -nologo                 # no logo
    -standard-semantics     # determines whether the current Fortran Standard behavior of the compiler is fully implemented.
    -integer-size 32        # Use 32-bits representation as default integer/logical kind.
    -diag-disable=7025      # This directive is not standard Fxx.
    -diag-disable=5268      # Extension to standard: The text exceeds right hand column allowed on the line.
    -diag-disable=10346     # optimization reporting will be enabled at link time when performing interprocedural optimizations.
    -diag-disable=10448     # Disable ifort deprecation message.
    -static-intel
    )
    if (MT_ENABLED)
        set(FCL_FLAGS "${FCL_FLAGS}"
        -threads
        )
    endif()
endif()

if (${codecov_enabled})
    set(FCL_FLAGS "${FCL_FLAGS}" -prof-gen=srcpos)
    set(CCL_FLAGS "${CCL_FLAGS}" -prof-gen=srcpos)
    set(FL_FLAGS "${FL_FLAGS}" -prof-gen=srcpos)
endif()

if (${perfprof_enabled})
    message(FATAL_ERROR
            "\n"
            "${pmfatal} Performance profiling with Intel compiler is currently unsupported.\n"
            "\n"
            )
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set C/Fortran compiler/linker debug build flags.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (WIN32)
    #set(CCL_FLAGS_DEBUG "${CCL_FLAGS_DEBUG}"
    #/debug:full
    #/Zi
    #/Od
    #/Wall
    #/traceback
    #/Qcheck-pointers:rw                 # check bounds for memory access through pointers.
    #/Qcheck-pointers-undimensioned      # check bounds for memory access through arrays that are declared without dimensions.
    #/Qdiag-error-limit:10
    #/Qtrapuv
    #)
    set(FCL_FLAGS_DEBUG "${FCL_FLAGS_DEBUG}"
    /debug:full
    /stand:f18
    /Zi
    /CB
    /Od
    /Qfp-stack-check
    /Qinit:snan,arrays
    /warn:all
    /gen-interfaces
    /traceback
    /check:all
    /check:bounds
    /fpe:0
   #/fpe-all:0
    /Qdiag-disable:5268
    /Qdiag-disable:7025
    /Qdiag-disable=10346
    /Qdiag-error-limit:10
    /Qtrapuv
    )
else()
    set(FCL_FLAGS_DEBUG "${FCL_FLAGS_DEBUG}"
    -stand f18                          #   Tells the compiler to issue compile-time messages for nonstandard language elements.
    -debug full                         #   generate full debug information
    -g3                                 #   generate full debug information
    -O0                                 #   disable optimizations
    -CB                                 #   Perform run-time bound-checks on array subscript and substring references (same as the -check bounds option)
   #-fp-stack-check                     #   only ifort: Generate extra code after every function call to ensure that the floating-point stack is in the expected state. This feature is only available for ifort.
    -init:snan,arrays                   #   initialize arrays and scalars to NaN
    -warn all                           #   enable all warning: -warn declarations enables warnings about implicit typing.
    -gen-interfaces                     #   generate interface block for each routine in the source files
    -traceback                          #   trace back for debugging
    -check all                          #   check all
    -check bounds                       #   check array bounds
    -fpe0                               #   Ignore underflow (yield 0.0); Abort on other IEEE exceptions.
   #-fpe-all=0                          #   Floating-point invalid, divide-by-zero, and overflow exceptions are enabled
    -diag-disable=5268                  #   Extension to standard: The text exceeds right hand column allowed on the line.
    -diag-disable=7025                  #   This directive is not standard Fxx.
    -diag-disable=10346                 #   optimization reporting will be enabled at link time when performing interprocedural optimizations.
    -diag-error-limit=10                #   max diagnostic error count
    -ftrapuv                            #   Initializes stack local variables to an unusual value to aid error detection.
    )
    #set(CCL_FLAGS_INTEL "${CCL_FLAGS_INTEL}"
    #-debug full
    #-g3
    #-O0
    #-Wall
    #-traceback
    #-check-pointers=rw                  #   check bounds for memory access through pointers.
    #-check-pointers-undimensioned       #   check bounds for memory access through arrays that are declared without dimensions.
    #-diag-error-limit=10
    #-ftrapuv
    #)
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set C/Fortran compiler/linker release build flags.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if(WIN32)
    #   Do NOT add interprocedural optimization flags to release modes in this file.
    #   These are generic flags that also apply to tests.
    #   Specifying IPO/LTO for tests lengthens tests compilation to hours.
    # INTEL_Fortran_RELEASE_FLAGS=/fast /O3 /Qip /Qipo /Qunroll /Qunroll-aggressive /inline:all /Ob2 /Qparallel /Qinline-dllimport
    # set(CCL_FLAGS_INTEL "${CCL_FLAGS_INTEL}"
    # /O3
    # /Ob2
    # /Qip
    # /Qipo
    # /Qunroll
    # /Qinline-dllimport
    # /Qunroll-aggressive
    # # /Qparallel  # Tells the auto-parallelizer to generate multithreaded code for loops that can be safely executed in parallel.
    # # /Qipo-c generates a single object file, containing all object files.
    # )
    set(FCL_FLAGS_RELEASE "${FCL_FLAGS_RELEASE}"
    /O3                             # Enable O3 optimization.
    /Qftz                           # Flushes subnormal results to zero. It may improve performance if the subnormal values are not critical to your application's behavior.
    /Qunroll                        # [:n] set the maximum number of times to unroll loops (no number n means automatic).
    /Qunroll-aggressive             # use more aggressive unrolling for certain loops.
    /Qinline-forceinline            # Instructs the compiler to force inlining of functions suggested for inlining whenever the compiler is capable doing so.
    /Qdiag-disable:10448            # Disables "ifort: remark #10448: Intel(R) Fortran Compiler Classic (ifort) is now deprecated and will be discontinued late 2024."
   #/Qvec                           # enable vectorization. This is now the default as of Intel 2019. Use /Qax to specify multiple target archs.
   #/Qguide-vec:4                   # enable guidance for auto-vectorization, causing the compiler to generate messages suggesting ways to improve optimization (default=4, highest).
   #/Qparallel                      # generate multithreaded code for loops that can be safely executed in parallel.
   #/Qipo-c:                        # Tells the compiler to optimize across multiple files and generate a single object file ipo_out.obj without linking
                                    # info at: https://software.intel.com/en-us/Fortran-compiler-developer-guide-and-reference-ipo-c-qipo-c
    )
else()
    set(FCL_FLAGS_RELEASE "${FCL_FLAGS_RELEASE}"
    -O3                             # set the optimizations level
    -unroll                         # [=n] set the maximum number of times to unroll loops (no number n means automatic).
   #-unroll-aggressive              # use more aggressive unrolling for certain loops. This is available only in fort and so excluded as of Dec 2023 because ifort is being phased out.
    -diag-disable=10346             # optimization reporting will be enabled at link time when performing interprocedural optimizations.
    -diag-disable=10397             # optimization reporting will be enabled at link time when performing interprocedural optimizations.
   #-guide-vec=4                    # enable guidance for auto-vectorization, causing the compiler to generate messages suggesting ways to improve optimization (default=4, highest).
   #-parallel                       # generate multithreaded code for loops that can be safely executed in parallel. This option requires MKL libraries.
   #-qopt-subscript-in-range        # assumes there are no "large" integers being used or being computed inside loops. A "large" integer is typically > 2^31.
    -ftz                            # Flushes subnormal results to zero.
    )
    #set(CCL_FLAGS_INTEL "${CCL_FLAGS_INTEL}"
    #-O3
    #-ip
    #-ipo
    #-unroll
    #-unroll-aggressive
    #-inline-level=2
    ## -parallel
    #)
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set C/Fortran compiler/linker testing build flags. Ideally, this should bet set to release + debug flags.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (WIN32)
    #set(CCL_FLAGS_TESTING "${CCL_FLAGS_TESTING}"
    #/Od  # disable optimization
    #)
    set(FCL_FLAGS_TESTING "${FCL_FLAGS_TESTING}"
    /Od # disable optimization
    )
else()
    set(FCL_FLAGS_TESTING "${FCL_FLAGS_TESTING}"
    -O0 # disable optimization
    )
    #set(CCL_FLAGS_TESTING "${CCL_FLAGS_TESTING}"
    #-O0  # disable optimization
    #)
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Set up coarray flags.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if ("${CAF_ENABLED}")
    message(NOTICE "${pmattn} Setting up Coarray Fortran (CAF) parallelization model. Options: single, shared, distributed" )
    message(NOTICE "${pmattn} Coarray Fortran (CAF) parallelization model enabled: ${CAF_ENABLED}" )
    message(NOTICE "${pmattn} Requested CAF parallelism: ${CAF_TYPE}")
    if ("${CAF_TYPE}" STREQUAL "single" OR
        "${CAF_TYPE}" STREQUAL "shared" OR
        "${CAF_TYPE}" STREQUAL "distributed")
        message(NOTICE "${pmattn} Enabling Intel Coarray Fortran parallelism: ${CAF_TYPE}")
        if (WIN32)
            set(FCL_FLAGS "${FCL_FLAGS}" /Qcoarray=${CAF_TYPE})
        else()
            set(FCL_FLAGS "${FCL_FLAGS}" -coarray=${CAF_TYPE})
        endif()
    else()
        message(FATAL_ERROR "${pmfatal} Internal error occurred. Unrecognized Coarray Fortran parallelization model: CAF_TYPE=${CAF_TYPE}\n")
    endif()
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Set the OpenMP parallelization flags.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (OMP_ENABLED)
    if(WIN32)
        set(FCL_FLAGS "${FCL_FLAGS}"
        /Qparallel          # requires openmp
        /Qopenmp            # Qiopenmp
        )
        if (NOT "${build}" STREQUAL "debug")
            set(FCL_FLAGS "${FCL_FLAGS}"
            /Qopt-matmul    # requires openmp
            )
        endif()
    else()
        set(FCL_FLAGS "${FCL_FLAGS}"
        -qopenmp-link=static
        -parallel           # requires openmp
        -qopenmp            # fiopenmp
        )
        if (NOT "${build}" STREQUAL "debug")
            set(FCL_FLAGS "${FCL_FLAGS}"
            -qopt-matmul    # requires openmp
            )
        endif()
    endif()
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Set memory allocation type. This flag is not recognized by the linker on macOS. So it must be used only at compile step.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (HEAP_ARRAY_ENABLED)
    if (WIN32)
        set(FC_FLAGS "${FC_FLAGS}"
        /heap-arrays#:10 # in KB # WARNING: Do NOT set the minimum heap array size. Otherwise, the various routines will seg-fault spectacularly at runtime.
        )
    else()
        set(FC_FLAGS "${FC_FLAGS}"
        -heap-arrays#=10 # in KB # WARNING: Do NOT set the minimum heap array size. Otherwise, the various routines will seg-fault spectacularly at runtime.
        )
    endif()
endif()
