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
!  
!  Description:
!       - The Fortran module source file containing the prototypes of the ParaMonte library routines to be called from Fortran.
!  Prototypes:
!       - runParaDRAM   : The procedural interface to the ParaDRAM sampler routine. This is the basic default interface.
!       - ParaDRAM_type : The ParaDRAM sampler class (derived type).
!       - ParaDRAM      : An alias for ParaDRAM_type().
!  Author:
!       - Computational Data Science Lab, Monday 9:03 AM, May 16 2016, ICES, UT Austin
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module paramonte

    ! ParaMonte default integer, real, and complex kinds are defined by IK, RK, and CK below.
    ! Use IK and RK from this module for integer and real kind specifications in your Fortran codes
    ! to ensure the highest level of consistency with the integer and real kinds of the ParaMonte library.

#if defined IS_COMPATIBLE_COMPILER

    ! Use the identifiers here to access the Object-Oriented interface to the ParaMonte routines.
    ! To enable the object-oriented interface, you must pass the IS_COMPATIBLE_COMPILER 
    ! preprocessor flag to the compiler at the time of compiling this module.

    use ParaDRAM_mod, only: IK, RK, CK
    use ParaDRAM_mod, only: getLogFunc_proc
    use ParaDRAM_mod, only: getLogFunc_interface => getLogFunc_proc
    use ParaDRAM_mod, only: ParaDRAM => ParaDRAM_type
    use ParaDRAM_mod, only: ParaDRAM_type

#else

    ! This is the procedural interface to the ParaDRAM sampler routine of the ParaMonte library.
    ! By default, if the preprocessor flag IS_COMPATIBLE_COMPILER is not passed to the compiler,
    ! the following procedural interface will be used.

    use, intrinsic :: iso_fortran_env, only: IK => int32, RK => real64
    implicit none

    ! The Fortran objective function interface (getLogFunc). Here, `proc` stands for the procedure interface.

    abstract interface
        function getLogFunc_proc(ndim,Point) result(logFunc)
            import :: IK, RK
            integer(IK) , intent(in)    :: ndim
            real(RK)    , intent(in)    :: Point(ndim)
            real(RK)                    :: logFunc
        end function getLogFunc_proc
    end interface

    interface
        subroutine runParaDRAM  ( ndim          &
                                , getLogFunc    &
                                , inputFile     &
                                )
            import :: IK, getLogFunc_proc
            implicit none
            integer(IK) , intent(in)            :: ndim
            procedure(getLogFunc_proc)          :: getLogFunc
            character(*), intent(in), optional  :: inputFile
        end subroutine runParaDRAM
    end interface

#endif

end module paramonte