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
#: Set the default Fortran compiler flags for all build modes.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# -std=legacy is required to bypass the new gfortran 10 error on argument-mismatch.
# The MPICH 3.2 library mpi_bcast still requires argument-mismatch, which causes gfortran to break the compilation.
# The problem still persists in debug mode. Therefore, when gfortran is 10, debug mode is disabled.
#set(FCL_FLAGS -std=gnu -ffree-line-length-none -fallow-argument-mismatch CACHE STRING "GNU Fortran default compiler flags" )

set(FCL_FLAGS "${FCL_FLAGS}"
    -ffree-line-length-none
    -fimplicit-none
    -std=legacy
    )

set(CCL_FLAGS "${CCL_FLAGS}"
    -ffree-line-length-none
    )

if (MT_ENABLED)
    set(FCL_FLAGS "${FCL_FLAGS}" -pthread)
    set(CCL_FLAGS "${CCL_FLAGS}" -pthread)
endif()

if (${codecov_enabled})
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

if (${perfprof_enabled})
    #set(CCL_FLAGS "${CCL_FLAGS}" -pg -g) # enable profiling + enable debugging which is required for line-by-line profiling with gprof
    set(FCL_FLAGS "${FCL_FLAGS}" -pg -g) # enable profiling + enable debugging which is required for line-by-line profiling with gprof
    set(FL_FLAGS "${FL_FLAGS}" -pg -g)   # enable profiling + enable debugging which is required for line-by-line profiling with gprof
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Set the preprocessed source file generation flag.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (FPP_ENABLED)
    # cwd option is important for consistency with Intel ifort behavior both
    # through CMake to store all temporary files directly in paramonte_bld_obj_dir.
    set(FC_FLAGS "${FC_FLAGS}" -cpp -save-temps=cwd) # enable preprocessed source file generation.
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set C/Fortran compiler/linker debug build flags.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#set(CCL_FLAGS_DEBUG "${CCL_FLAGS_DEBUG}"
#   -g                                  # generate full debug information
#   -O0                                 # disable optimizations
#   -fcheck=all                         # enable the generation of run-time checks
#   -ffpe-trap=invalid,zero,overflow    # ,underflow : Floating-point invalid, divide-by-zero, and overflow exceptions are enabled
#   -finit-real=snan                    # initialize REAL and COMPLEX variables with a signaling NaN
#   -fbacktrace                         # trace back for debugging
#   -pedantic                           # issue warnings for uses of extensions to the Fortran standard
#   -fmax-errors=10                     # max diagnostic error count
#   -Wall                               # enable all warnings:
#                                       # -Waliasing, -Wampersand, -Wconversion, -Wsurprising, -Wc-binding-type, -Wintrinsics-std, -Wtabs, -Wintrinsic-shadow,
#                                       # -Wline-truncation, -Wtarget-lifetime, -Winteger-division, -Wreal-q-constant, -Wunused, -Wundefined-do-loop
#   )

set(FCL_FLAGS_DEBUG "${FCL_FLAGS_DEBUG}"
    -g3                                 # generate full debug information
    -O0                                 # disable optimizations
   #-fsanitize=undefined                # enable UndefinedBehaviorSanitizer for undefined behavior detection.
   #-fsanitize=address                  # enable AddressSanitizer, for memory error detection, like out-of-bounds and use-after-free bugs.
   #-fsanitize=leak                     # enable LeakSanitizer for memory leak detection.
    -fcheck=all                         # enable the generation of run-time checks
    -ffpe-trap=invalid,zero,overflow    # ,underflow : Floating-point invalid, divide-by-zero, and overflow exceptions are enabled
    -ffpe-summary=all                   # Specify a list of floating-point exceptions, whose flag status is printed to ERROR_UNIT when invoking STOP and ERROR STOP.
                                        # Can be either ‘none’, ‘all’ or a comma-separated list of the following exceptions:
                                        # ‘invalid’, ‘zero’, ‘overflow’, ‘underflow’, ‘inexact’ and ‘denormal’
    -finit-integer=-2147483647          # initilize all integers to negative infinity
    -finit-real=snan                    # initialize REAL and COMPLEX variables with a signaling NaN
    -fbacktrace                         # trace back for debugging
   #-pedantic                           # issue warnings for uses of extensions to the Fortran standard. Gfortran10 with MPICH 3.2 in debug mode crashes with this flag at mpi_bcast. Excluded until MPICH upgraded.
    -fmax-errors=10                     # max diagnostic error count
    -Wno-maybe-uninitialized            # avoid warning of no array pre-allocation.
    -Wall                               # enable all warnings:
                                        # -Waliasing, -Wampersand, -Wconversion, -Wsurprising, -Wc-binding-type, -Wintrinsics-std, -Wtabs, -Wintrinsic-shadow,
                                        # -Wline-truncation, -Wtarget-lifetime, -Winteger-division, -Wreal-q-constant, -Wunused, -Wundefined-do-loop
                                        # gfortran10 crashes and cannot compile MPI ParaMonte with mpich in debug mode. Therefore -wall is disabled for now, until MPICH upgrades interface.
   #-Wconversion-extra                  # Warn about implicit conversions between different types and kinds. This option does not imply -Wconversion.
   #-Werror=conversion                  # Turn all implicit conversions into an error. This is important to avoid inadvertent implicit change of precision in generic procedures of various kinds, due to the use of `RK` to represent different kinds.
   #-Werror=conversion-extra            # Turn all implicit conversions into an error. This is too aggressive and as such currently deactivated. For example, it yields an error on the multiplication of integer with real.
    -Wno-surprising                     # -Wsurpring yields many false positives like "Array x at (1) is larger than limit set by '-fmax-stack-var-size='".
    -fno-unsafe-math-optimizations
    -fsignaling-nans
    -frounding-math
    #-Waliasing
    #-Wampersand
    #-Wconversion
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

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set C/Fortran compiler/linker testing build flags.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#set(CCL_FLAGS_TESTING "${CCL_FLAGS_TESTING}"
#    -O0  # disable optimization
#    )

set(FCL_FLAGS_TESTING "${FCL_FLAGS_TESTING}"
    #-static-libgfortran -static-libgcc
    -O0  # disable optimization
    )

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set C/Fortran compiler/linker release build flags.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#set(CCL_FLAGS_TESTING "${CCL_FLAGS_TESTING}"
#    -O3                         # set the optimizations level
#    -flto                       # enable interprocedural optimization between files.
#    -funroll-loops              # [=n] set the maximum number of times to unroll loops (no number n means automatic).
#    -finline-functions          # consider all functions for inlining, even if they are not declared inline.
#    -ftree-vectorize            # perform vectorization on trees. enables -ftree-loop-vectorize and -ftree-slp-vectorize.
#    )

if(APPLE)
    # GNU 9.3 -O results in runtime crashes of the examples on Mac (except -O0)
    set(FCL_FLAGS_RELEASE "${FCL_FLAGS_RELEASE}"
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
    set(FCL_FLAGS_RELEASE "${FCL_FLAGS_RELEASE}"
        #-static-libgfortran -static-libgcc
        -ftree-vectorize        # perform vectorization on trees. enables -ftree-loop-vectorize and -ftree-slp-vectorize.
        -funroll-loops          # [=n] set the maximum number of times to unroll loops (no number n means automatic).
        -O3                     # set the optimizations level
        )
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Set up coarray flags.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

message(NOTICE "${pmattn} Setting up Coarray Fortran (CAF) parallelization model. Options: single, shared, distributed" )
message(NOTICE "${pmattn} Coarray Fortran (CAF) parallelization model enabled: ${CAF_ENABLED}" )
message(NOTICE "${pmattn} Requested CAF parallelism: ${CAF_TYPE}")

if ("${CAF_ENABLED}")
    if ("${CAF_TYPE}" STREQUAL "single")
        message(NOTICE "${pmattn} Enabling Intel Coarray Fortran parallelism: ${CAF_TYPE}")
        set(FCL_FLAGS "${FCL_FLAGS}" -coarray="${CAF_TYPE}")
    elseif (NOT ("${CAF_TYPE}" STREQUAL "shared" OR "${CAF_TYPE}" STREQUAL "distributed"))
        message(FATAL_ERROR "${pmfatal} Internal error occurred. Unrecognized Coarray Fortran parallelization model: CAF_TYPE=${CAF_TYPE}\n")
    endif()
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Set the OpenMP parallelization flags.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (OMP_ENABLED)
    set(FCL_FLAGS "${FCL_FLAGS}"
        #-ftree-parallelize-loops=
        -floop-parallelize-all
        -frecursive
        -fopenmp
        )
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Set memory allocation type. This flag is not recognized by the linker on macOS. So it must be used only at compile step.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#   \todo
#   Strangely, -frecursive flag causes deadlock with sampler tests in Coarray mode.
#   Therefore, the -frecursive flag is reversed to -fmax-stack-var-size=10 until the behavior is understood.
#   The use of -frecursive was based on the GFortran-10 warning message about unsafe storage move from stack to static storage.
#   See also this thread: https://comp.lang.fortran.narkive.com/WApl1KMt/gfortran-stack-size-warning#post4
#   Perhaps adding `non_recursive` to functions would fix this warning message.

if (HEAP_ARRAY_ENABLED)
    set(FC_FLAGS "${FC_FLAGS}"
        -fmax-stack-var-size=10 # in bytes
        # -frecursive overwrites -fmax-stack-var-size=10 and causes all allocations to happen on stack.
        )
else()
    set(FC_FLAGS "${FC_FLAGS}"
        -fstack-arrays # applies to only the local variables
        # -frecursive overwrites -fmax-stack-var-size=10 and causes all allocations to happen on stack.
        # -fautomatic
        )
endif()
