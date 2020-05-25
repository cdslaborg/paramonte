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

! This module defines the abstract interface of the objective function to be called by ParaMonte routines to sample it.

module ParaMonteLogFunc_mod

    use Constants_mod, only : IK, RK

    ! Fortran interface to the objective function implementation (getLogFunc). Here, `proc` stands for the procedure interface.

    abstract interface
#if defined CFI_ENABLED
        function getLogFunc_proc(ndim,Point) result(logFunc) bind(C)
#else
        function getLogFunc_proc(ndim,Point) result(logFunc)
#endif
            import :: IK, RK
#if defined CFI_ENABLED
            integer(IK), intent(in), value  :: ndim
#else
            integer(IK), intent(in)         :: ndim
#endif
            real(RK), intent(in)            :: Point(ndim)
            real(RK)                        :: logFunc
        end function getLogFunc_proc
    end interface

    ! Fortran interface for the gradient of the objective function

    abstract interface
        function getGradLogFunc_proc(ndim,Point) result(GradLogFunc)
            import :: IK, RK
            integer(IK), intent(in) :: ndim
            real(RK)   , intent(in) :: Point(ndim)
            real(RK)                :: GradLogFunc(ndim)
        end function getGradLogFunc_proc
    end interface


    ! C-interoperable interface for the gradient of the objective function

#if defined CFI_ENABLED
    abstract interface
        subroutine getGradLogFunc4C_proc(ndim,Point,GradLogFunc) bind(C)
            import :: IK, RK
            integer(IK), intent(in), value  :: ndim
            real(RK)   , intent(in)         :: Point(ndim)
            real(RK)   , intent(out)        :: GradLogFunc(ndim)
        end subroutine getGradLogFunc4C_proc
    end interface
#endif

end module ParaMonteLogFunc_mod
