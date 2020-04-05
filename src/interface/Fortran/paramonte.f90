!**********************************************************************************************************************************
!**********************************************************************************************************************************
!
!  ParaMonte: plain powerful parallel Monte Carlo library.
!
!  Copyright (C) 2012-present, The Computational Data Science Lab
!
!  This file is part of the ParaMonte library. 
!
!  ParaMonte is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, version 3 of the License.
!
!  ParaMonte is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!**********************************************************************************************************************************
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
    ! to ensure the highest level of consistency with the ParaMonte library's integer and real kinds. 

#if defined IS_COMPATIBLE_COMPILER

    use ParaDRAM_mod, only: IK, RK, CK
    use ParaDRAM_mod, only: getLogFunc_proc
    use ParaDRAM_mod, only: getLogFunc_interface => getLogFunc_proc
    use ParaDRAM_mod, only: ParaDRAM => ParaDRAM_type
    use ParaDRAM_mod, only: ParaDRAM_type

#else

    use, intrinsic :: iso_fortran_env, only: IK => int32, RK => real64
    implicit none

    ! The procedural interface to the ParaDRAM sampler routine of the ParaMonte library.
    ! For the object oriented interface, you must pass the IS_COMPATIBLE_COMPILER 
    ! preprocessor flag to the compiler at the time of comipling this module.

    interface
        subroutine runParaDRAM  ( ndim          &
                                , getLogFunc    &
                                , inputFile     &
                                ) bind(C, name="runParaDRAM")
            implicit none
            integer(IK) , intent(in)    :: ndim
            character(*), intent(in)    :: inputFile
            procedure(getLogFunc_proc)  :: getLogFunc
        end subroutine runParaDRAM
    end interface

    ! The Fortran objective function interface (getLogFunc). Here, `proc` stands for the procedure interface.

    abstract interface
        function getLogFunc_proc(ndim,Point) result(logFunc) bind(C)
            import :: IK, RK
            integer(IK) , intent(in)    :: ndim
            real(RK)    , intent(in)    :: Point(ndim)
            real(RK)                    :: logFunc
        end function getLogFunc_proc
    end interface

#endif

end module paramonte