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

# This Bash file configures the flags required for building ParaMonte library, tests, and examples on non-Windows Operating Systems.

FILE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

#echo " "
#echo "Configuring ParaMonte Build..."
#echo "Configuration File: " $FILE_DIR

if [ -z ${ParaMonteTest_RUN_ENABLED+x} ]
then
    ParaMonteTest_RUN_ENABLED=false
fi
export ParaMonteTest_RUN_ENABLED

if [ -z ${ParaMonteExample_RUN_ENABLED+x} ]
then
    ParaMonteExample_RUN_ENABLED=false
fi
export ParaMonteExample_RUN_ENABLED

echo " "

####################################################################################################################################

#    set interface language:   c, fortran, matlab, python (all must be lower case)

if [ -z ${INTERFACE_LANGUAGE+x} ]; then
    ############################################################
    INTERFACE_LANGUAGE=c
    ############################################################
    export INTERFACE_LANGUAGE
fi
echo "    interface language: $INTERFACE_LANGUAGE"

####################################################################################################################################

# unset compiler suite: intel, gnu

if [ -z ${PMCS+x} ]
then
    ############################################################
    #unset PMCS
    ############################################################
    export PMCS
fi
echo "                 Compiler suite: $PMCS"

####################################################################################################################################

# define build type:
#            relwithdebinfo:   Release with Debug info.
#                   release:   Full-blown highly optimized production-level library build.
#                     debug:   Used only for development step. No optimizations performed.

if [ -z ${BTYPE+x} ]
then
    ############################################################
    BTYPE=release
    ############################################################
    export BTYPE
fi
echo "                     Build type: $BTYPE"

####################################################################################################################################

# define linking type:
#                   dynamic:   Use this flag when you have R/Python/MATLAB/Julia code to which you need to link the ParaMonte library dynamically, using DLL files. 
#                    static:   Use this flag when you have Fortran/C/C++ code to which you want to link the ParaMonte library statically.
#                              You can also link dynamically your Fortran/C/C++ codes using DLL files by specifying LTYPE=dynamic flag instead.

if [ -z ${LTYPE+x} ]; then
    ############################################################
    LTYPE=dynamic
    ############################################################
    export LTYPE
fi
echo "         ParaMonte library type: $LTYPE"

####################################################################################################################################

# set Coarray Fortran (CAF) parallelization model:
#                      none:   No Coarray Parallelization will be performed. Serial library will be built.   
#                    single:   Coarray Parallelism will be invoked. However, only one image will perform the tasks, and there is no parallel communications.
#                    shared:   Coarray Parallelism will be invoked. It causes the underlying Intel® Message Passing Interface (MPI) parallelization
#                              to run on one node with multiple cores or processors with shared memory.
#               distributed:   This option requires a special license to be installed for Intel compiler suite, and is only available on Linux systems, although you can specify it here.
#                              On the Linux systems, it causes the underlying Intel® MPI Library parallelization to run in a multi-node environment (multiple CPUs with distributed memory).

if [ -z ${CAFTYPE+x} ]; then
    ############################################################
    CAFTYPE=none
    ############################################################
    export CAFTYPE
fi
echo "       Coarray parallelism type: $CAFTYPE"

####################################################################################################################################

# set MPI parallelization model flags:

if [ -z ${MPI_ENABLED+x} ]; then
    ############################################################
    MPI_ENABLED=true
    ############################################################
    export MPI_ENABLED
fi
echo "        MPI parallelism enabled: $MPI_ENABLED"

####################################################################################################################################

# set OpenMP parallelization model flags:

if [ -z ${OMP_ENABLED+x} ]; then
    ############################################################
    OMP_ENABLED=false
    ############################################################
    export OMP_ENABLED
fi
echo "     OpenMP parallelism enabled: $OMP_ENABLED"

####################################################################################################################################

# set interoperability mode:
#                      true:   When you are calling ParaMonte library from any language other than Fortran.
#                     false:   When your objective function to be sampled by ParaMonte library is written in Fortran (as opposed to C).
#                              It is also fine to set this flag to true for Fortran application. However, the syntax of the Fortran objective function that
#                              is passed to ParaMonte library will have to conform to the rules of C-Fortran interoperability standard,
#                              as given in the abstract interface in the ParaMonte source file: ParaMonteLogFunc_mod.f90

if [ -z ${CFI_ENABLED+x} ]; then
    ############################################################
    CFI_ENABLED=true
    ############################################################
    export CFI_ENABLED
fi
echo "    C-Fortran interface enabled: $CFI_ENABLED"

####################################################################################################################################

# set Fortran/C array allocation resource.
#                      true:   All automatic and stack arrays will be allocated on the heap. Use this when calling ParaMonte from other languages.
#                     false:   All automatic and stack arrays will be allocated on the stack. This could lead to stack overflow when calling DLL from other languages.

if [ -z ${HEAP_ARRAY_ENABLED+x} ]; then
    ############################################################
    HEAP_ARRAY_ENABLED=true
    ############################################################
    export HEAP_ARRAY_ENABLED
fi
echo "  Heap array allocation enabled: $HEAP_ARRAY_ENABLED"

####################################################################################################################################

# set number of Fortran Coarray images (that will run in parallel, if Coarray parallel programming is requested by user)

if [ -z ${FOR_COARRAY_NUM_IMAGES+x} ]; then
    ############################################################
    FOR_COARRAY_NUM_IMAGES=3
    ############################################################
    export FOR_COARRAY_NUM_IMAGES
fi
echo "        Default number of cores: $FOR_COARRAY_NUM_IMAGES"

echo " "

####################################################################################################################################
