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

module Test_EconomicsToolbox_mod

    use EconomicsToolbox_mod
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_EconomicsToolbox

    type(Test_type) :: Test

    type :: TestFuncRosenBrock1D_type
        ! MATLAB reference code:
        ! function f = objectivefcn1(x)
        !     f = exp(-x^2) * cos(x) *sin(2*x);
        ! end
        ! CMD: x0 = [0.25];
        ! CMD: options = optimset('tolx',1.e-10)
        ! CMD: [x, f] = fminsearch(@objectivefcn1,x0,options)
        real(RK) :: XMIN_REF = -0.471354350447655_RK
        real(RK) :: FMIN_REF = -0.577293243101421_RK
    contains
        procedure, nopass :: get => getTestFuncRosenBrock1D
    end type TestFuncRosenBrock1D_type

    type :: TestFuncRosenBrock2D_type
        ! MATLAB reference code:
        ! function f = objectivefcn1(x)
        !     f = exp(-(x(1)-x(2))^2 - 2*x(1)^2)*cos(x(2))*sin(2*x(2));
        ! end
        ! CMD: x0 = [0.25,-0.25];
        ! CMD: options = optimset('tolx',1.e-10)
        ! CMD: [x, f] = fminsearch(@objectivefcn1,x0,options)
        real(RK) :: XMIN_REF(2) = [-0.169552635123540_RK, -0.508657910611664_RK]
        real(RK) :: FMIN_REF = -0.625285429865811_RK
    contains
        procedure, nopass :: get => getTestFuncRosenBrock2D
    end type TestFuncRosenBrock2D_type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_EconomicsToolbox()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call test_minimize()
        call Test%finalize()
#ifdef CAF_ENABLED
        sync all
#endif
    end subroutine test_EconomicsToolbox

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_minimize()

        use Constants_mod, only: RK, IK, NEGINF_RK, POSINF_RK
        implicit none

        real(RK) :: xmin(2), fmin
        type(TestFuncRosenBrock2D_type) :: TestFuncRosenBrock2D

        if (Test%Image%isFirst) call Test%testing("powell()")

        xmin = [0.25_RK,-0.25_RK]
        call powell ( p = xmin &
                    , fret = fmin &
                    , minimum = [NEGINF_RK, NEGINF_RK] &
                    , maximum = [POSINF_RK, POSINF_RK] &
                    , func = getTestFuncRosenBrock2D &
                    )

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "xmin    ", xmin
            write(Test%outputUnit,"(*(g0,:,', '))") "XMIN_REF", TestFuncRosenBrock2D%XMIN_REF
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "fmin    ", fmin
            write(Test%outputUnit,"(*(g0,:,', '))") "FMIN_REF", TestFuncRosenBrock2D%FMIN_REF
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if

        Test%assertion = all( abs(xmin-TestFuncRosenBrock2D%XMIN_REF) / abs(TestFuncRosenBrock2D%XMIN_REF) < 1.e-6_RK )
        call Test%verify()
        Test%assertion = abs(fmin-TestFuncRosenBrock2D%FMIN_REF) / abs(TestFuncRosenBrock2D%FMIN_REF) < 1.e-6_RK
        call Test%verify()

    end subroutine test_minimize

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getTestFuncRosenBrock1D(x)
        use Constants_mod, only: RK
        implicit none
        real(RK), intent(in)    :: x
        real(RK)                :: getTestFuncRosenBrock1D
        getTestFuncRosenBrock1D = exp(-x**2) * cos(x) * sin(2*x)
    end function getTestFuncRosenBrock1D

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getTestFuncRosenBrock2D(x)
        use Constants_mod, only: RK
        implicit none
        real(RK), intent(in)    :: x(:)
        real(RK)                :: getTestFuncRosenBrock2D
        getTestFuncRosenBrock2D = exp(-(x(1)-x(2))**2 - 2*x(1)**2) * cos(x(2)) * sin(2*x(2))
    end function getTestFuncRosenBrock2D

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_EconomicsToolbox_mod