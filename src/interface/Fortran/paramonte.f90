!***********************************************************************************************************************************
!***********************************************************************************************************************************
!
!   ParaMonte: plain powerful parallel Monte Carlo library.
!
!   Copyright (C) 2012-present, The Computational Data Science Lab
!
!   This file is part of the ParaMonte library.
!
!   ParaMonte is free software: you can redistribute it and/or modify it
!   under the terms of the GNU Lesser General Public License as published
!   by the Free Software Foundation, version 3 of the License.
!
!   ParaMonte is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with the ParaMonte library. If not, see,
!
!       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
!
!   ACKNOWLEDGMENT
!
!   As per the ParaMonte library license agreement terms,
!   if you use any parts of this library for any purposes,
!   we ask you to acknowledge the use of the ParaMonte library
!   in your work (education/research/industry/development/...)
!   by citing the ParaMonte library as described on this page:
!
!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!
!***********************************************************************************************************************************
!***********************************************************************************************************************************
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
!**********************************************************************************************************************************
!**********************************************************************************************************************************

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

    interface
        subroutine runParaDRAM  ( ndim          &
                                , getLogFunc    &
                                , inputFile     &
                                )
            implicit none
            integer(IK) , intent(in)            :: ndim
            procedure(getLogFunc_proc)          :: getLogFunc
            character(*), intent(in), optional  :: inputFile
        end subroutine runParaDRAM
    end interface

    ! The Fortran objective function interface (getLogFunc). Here, `proc` stands for the procedure interface.

    abstract interface
        function getLogFunc_proc(ndim,Point) result(logFunc)
            import :: IK, RK
            integer(IK) , intent(in)    :: ndim
            real(RK)    , intent(in)    :: Point(ndim)
            real(RK)                    :: logFunc
        end function getLogFunc_proc
    end interface

#endif

end module paramonte