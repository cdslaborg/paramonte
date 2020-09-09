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

# Find LAPACK (finds BLAS also) if not already found

IF(NOT LAPACK_FOUND)
    ENABLE_LANGUAGE(C) # Some libraries need a C compiler to find
    FIND_PACKAGE(LAPACK REQUIRED)
    # Remember that LAPACK (and BLAS) was found.  For some reason the
    # FindLAPACK routine doesn't place these into the CACHE.
    SET(BLAS_FOUND TRUE CACHE INTERNAL "BLAS was found" FORCE)
    SET(LAPACK_FOUND TRUE CACHE INTERNAL "LAPACK was found" FORCE)
    SET(BLAS_LIBRARIES ${BLAS_LIBRARIES} CACHE INTERNAL "BLAS LIBS" FORCE)
    SET(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} CACHE INTERNAL "LAPACK LIBS" FORCE)
ENDIF(NOT LAPACK_FOUND)
