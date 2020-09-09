!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!!   MIT License
!!!!
!!!!   ParaMonte: plain powerful parallel Monte Carlo library.
!!!!
!!!!   Copyright (C) 2012-present, The Computational Data Science Lab
!!!!
!!!!   This file is part of the ParaMonte library.
!!!!
!!!!   Permission is hereby granted, free of charge, to any person obtaining a 
!!!!   copy of this software and associated documentation files (the "Software"), 
!!!!   to deal in the Software without restriction, including without limitation 
!!!!   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
!!!!   and/or sell copies of the Software, and to permit persons to whom the 
!!!!   Software is furnished to do so, subject to the following conditions:
!!!!
!!!!   The above copyright notice and this permission notice shall be 
!!!!   included in all copies or substantial portions of the Software.
!!!!
!!!!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
!!!!   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
!!!!   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
!!!!   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
!!!!   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
!!!!   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
!!!!   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!!!!
!!!!   ACKNOWLEDGMENT
!!!!
!!!!   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
!!!!   As per the ParaMonte library license agreement terms, if you use any parts of 
!!!!   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
!!!!   work (education/research/industry/development/...) by citing the ParaMonte 
!!!!   library as described on this page:
!!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Keep in mind that Fortran is case-insensitive, except for character values and string values. 
! So, feel free to use upper-case or lower-case for the Fortran syntax and entity names.
! The ParaMonte library uses camelCase convention for variable naming.

#if defined IS_COMPATIBLE_COMPILER

    program main

    ! This is the Object-Oriented-Programming (OOP) style interface to the ParaMonte routines.
    ! This is more flexible but less portable, as its compilation requires the same compiler 
    ! brand and version with which the ParaMonte library has been built.

    use paramonte, only: ParaDRAM
    use LogFunc_mod, only: getLogFunc, NDIM

    implicit none
    type(ParaDRAM) :: pd

    call pd%runSampler  ( ndim = NDIM &
                        , getLogFunc = getLogFunc &
                        , inputFile = "./paramonte.in" &    ! this is optional argument
                        ! You can also specify simulation specifications as input arguments, like 
                        ! the following. This is possible only from the OOP interface to ParaDRAM.
                        , greedyAdaptationCount = 0 &       ! this is optional argument
                        , description = "an example run" &  ! this is optional argument
                        ! More optional arguments can appear here. 
                        ! See the ParaDRAM routine's list of input arguments.
                        )

    end program main

#else

    ! This is the default simple procedural interface to the ParaMonte routines.
    ! The first two arguments to the sampler routines are mandatory.
    ! The inputFile argument is optional.

    program main

    use paramonte, only: runParaDRAM
    use LogFunc_mod, only: getLogFunc, NDIM

    implicit none

    call runParaDRAM( NDIM &
                    , getLogFunc &
                    , "./paramonte.in" & ! this is optional argument
                    )

    end program main

#endif