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
