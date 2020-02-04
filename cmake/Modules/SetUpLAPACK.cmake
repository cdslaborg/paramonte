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
