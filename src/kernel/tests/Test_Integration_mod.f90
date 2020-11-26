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

!>  \brief This module contains tests of the module [Integration_mod](@ref integration_mod).
!>  @author Amir Shahmoradi

module Test_Integration_mod

    use Integration_mod
    !use QuadPackDouble_mod, only: dqagi
    use Constants_mod, only: IK, RK
    use Err_mod, only: Err_type
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_Integration

    type(Test_type) :: Test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_Integration()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_doQuadRombOpen_1, "test_doQuadRombOpen_1")
        call Test%run(test_doQuadRombOpen_2, "test_doQuadRombOpen_2")
        call Test%run(test_doQuadRombOpen_3, "test_doQuadRombOpen_3")
        call Test%run(test_doQuadRombOpen_4, "test_doQuadRombOpen_4")
        call Test%run(test_doQuadRombClosed_1, "test_doQuadRombClosed_1")
        call Test%finalize()
    end subroutine test_Integration

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! integrate the Gaussian PDF via the midexp routine on the open interval \f$[0, +\infty)\f$.
    function test_doQuadRombOpen_1() result(assertion)

        use Constants_mod, only: RK, IK
        use Integration_mod

        implicit none

        logical                 :: assertion
        real(RK)                :: integral, relativeError, difference 
        real(RK), parameter     :: integral_ref = 0.5_RK
        real(RK), parameter     :: tolerance = 1.e-10_RK
        integer(IK)             :: numFuncEval, ierr

        call doQuadRombOpen ( getFunc = getTestFuncOpenInterval_1 &
                            , integrate = midexp &
                            , lowerLim = 0._RK &
                            , upperLim = huge(0._RK) &
                            , maxRelativeError = 0.1*tolerance &
                            , nRefinement = 10_IK &
                            , integral = integral &
                            , relativeError = relativeError &
                            , numFuncEval = numFuncEval &
                            , ierr = ierr &
                            )

        if (ierr/=0_IK) then
            assertion = .false.
            if (Test%isDebugMode .and. .not. assertion) then
                write(Test%outputUnit,"(*(g0))") "ierr = ", ierr, " /= 0"
            end if
            return
        end if

        difference = abs(integral - integral_ref) / abs(integral_ref)
        assertion = difference < tolerance
        assertion = assertion .and. relativeError <= tolerance

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "integral_ref  = ", integral_ref
            write(Test%outputUnit,"(*(g0))") "integral      = ", integral
            write(Test%outputUnit,"(*(g0))") "difference    = ", difference
            write(Test%outputUnit,"(*(g0))") "relativeError = ", relativeError
            write(Test%outputUnit,"(*(g0))")
        end if

    end function test_doQuadRombOpen_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_doQuadRombOpen_2() result(assertion)

        use Constants_mod, only: RK, IK
        use Integration_mod

        implicit none

        logical                 :: assertion
        real(RK)                :: relativeError, difference
        real(RK)                :: integral
        real(RK), parameter     :: integral_ref = 0.078649603525143_RK
        real(RK), parameter     :: tolerance = 1.e-10_RK
        integer(IK)             :: numFuncEval, ierr

        call doQuadRombOpen ( getFunc = getTestFuncOpenInterval_1 &
                            , integrate = midinf &
                            , lowerLim = 1._RK &
                            , upperLim = huge(0._RK) &
                            , maxRelativeError = 0.1*tolerance &
                            , nRefinement = 10_IK &
                            , integral = integral &
                            , relativeError = relativeError &
                            , numFuncEval = numFuncEval &
                            , ierr = ierr &
                            )

        if (ierr/=0_IK) then
            assertion = .false.
            if (Test%isDebugMode .and. .not. assertion) then
                write(Test%outputUnit,"(*(g0))") "ierr = ", ierr, " /= 0"
            end if
            return
        end if

        difference = abs(integral - integral_ref) / abs(integral_ref)
        assertion = difference < tolerance
        assertion = assertion .and. relativeError <= tolerance

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "integral_ref  = ", integral_ref
            write(Test%outputUnit,"(*(g0))") "integral      = ", integral
            write(Test%outputUnit,"(*(g0))") "difference    = ", difference
            write(Test%outputUnit,"(*(g0))") "relativeError = ", relativeError
            write(Test%outputUnit,"(*(g0))")
        end if

    end function test_doQuadRombOpen_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_doQuadRombOpen_3() result(assertion)

        use Constants_mod, only: RK, IK
        use Integration_mod

        implicit none

        logical                 :: assertion
        real(RK)                :: relativeError, difference
        real(RK)                :: integral
        real(RK), parameter     :: integral_ref = 0.078649603525143_RK
        real(RK), parameter     :: tolerance = 1.e-4_RK
        integer(IK)             :: numFuncEval, ierr

        call doQuadRombOpen ( getFunc = getTestFuncOpenInterval_1 &
                            , integrate = midinf &
                            , lowerLim = -huge(0._RK) &
                            , upperLim = -1._RK &
                            , maxRelativeError = 1.e-10_RK &
                            , nRefinement = 10_IK &
                            , integral = integral &
                            , relativeError = relativeError &
                            , numFuncEval = numFuncEval &
                            , ierr = ierr &
                            )

        if (ierr/=0_IK) then
            assertion = .false.
            if (Test%isDebugMode .and. .not. assertion) then
                write(Test%outputUnit,"(*(g0))") "ierr = ", ierr, " /= 0"
            end if
            return
        end if

        difference = abs(integral - integral_ref) / abs(integral_ref)
        assertion = difference < tolerance
        assertion = assertion .and. relativeError <= tolerance

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "integral_ref  = ", integral_ref
            write(Test%outputUnit,"(*(g0))") "integral      = ", integral
            write(Test%outputUnit,"(*(g0))") "difference    = ", difference
            write(Test%outputUnit,"(*(g0))") "relativeError = ", relativeError
            write(Test%outputUnit,"(*(g0))")
        end if

    end function test_doQuadRombOpen_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_doQuadRombOpen_4() result(assertion)

        use Constants_mod, only: RK, IK
        use Integration_mod

        implicit none

        logical                 :: assertion
        real(RK)                :: relativeError, difference
        real(RK)                :: integral
        real(RK), parameter     :: integral_ref = 1._RK
        real(RK), parameter     :: tolerance = 1.e-10_RK
        integer(IK)             :: numFuncEval, ierr

        call doQuadRombOpen ( getFunc = getTestFuncOpenInterval_1 &
                            , integrate = midpnt &
                            , lowerLim = -5._RK &
                            , upperLim = +5._RK &
                            , maxRelativeError = 1.e-10_RK &
                            , nRefinement = 10_IK &
                            , integral = integral &
                            , relativeError = relativeError &
                            , numFuncEval = numFuncEval &
                            , ierr = ierr &
                            )

        if (ierr/=0_IK) then
            assertion = .false.
            if (Test%isDebugMode .and. .not. assertion) then
                write(Test%outputUnit,"(*(g0))") "ierr = ", ierr, " /= 0"
            end if
            return
        end if

        difference = abs(integral - integral_ref) / abs(integral_ref)
        assertion = difference < tolerance
        assertion = assertion .and. relativeError <= tolerance

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "integral_ref  = ", integral_ref
            write(Test%outputUnit,"(*(g0))") "integral      = ", integral
            write(Test%outputUnit,"(*(g0))") "difference    = ", difference
            write(Test%outputUnit,"(*(g0))") "relativeError = ", relativeError
            write(Test%outputUnit,"(*(g0))")
        end if

    end function test_doQuadRombOpen_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_doQuadRombClosed_1() result(assertion)

        use Constants_mod, only: RK, IK
        use Integration_mod

        implicit none

        logical                 :: assertion
        real(RK)                :: relativeError, difference
        real(RK)                :: integral
        real(RK), parameter     :: integral_ref = 1._RK
        real(RK), parameter     :: tolerance = 1.e-10_RK
        integer(IK)             :: numFuncEval, ierr

        call doQuadRombClosed   ( getFunc = getTestFuncOpenInterval_1 &
                                , lowerLim = -5._RK &
                                , upperLim = +5._RK &
                                , maxRelativeError = 0.1*tolerance &
                                , nRefinement = 10_IK &
                                , integral = integral &
                                , relativeError = relativeError &
                                , numFuncEval = numFuncEval &
                                , ierr = ierr &
                                )

        if (ierr/=0_IK) then
            assertion = .false.
            if (Test%isDebugMode .and. .not. assertion) then
                write(Test%outputUnit,"(*(g0))") "ierr = ", ierr, " /= 0"
            end if
            return
        end if

        difference = abs(integral - integral_ref) / abs(integral_ref)
        assertion = difference < tolerance
        assertion = assertion .and. relativeError <= tolerance

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "integral_ref  = ", integral_ref
            write(Test%outputUnit,"(*(g0))") "integral      = ", integral
            write(Test%outputUnit,"(*(g0))") "difference    = ", difference
            write(Test%outputUnit,"(*(g0))") "relativeError = ", relativeError
            write(Test%outputUnit,"(*(g0))")
        end if

    end function test_doQuadRombClosed_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! The Gaussian PDF with sigma = 1/sqrt(2)
    function getTestFuncOpenInterval_1(x) result(funcVal)
        use Constants_mod, only: RK, INVSQRTPI
        implicit none
        real(RK), intent(in)    :: x
        real(RK)                :: funcVal
        funcVal = INVSQRTPI * exp(-x**2)
    end function getTestFuncOpenInterval_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! requires `midsql`.
    ! The MATLAB test function: https://www.mathworks.com/help/matlab/ref/integral.html
    function getTestFuncOpenIntervalMATLAB(x) result(funcVal)
        use Constants_mod, only: RK
        implicit none
        real(RK), intent(in)    :: x
        real(RK)                :: funcVal
        funcVal = exp(-x**2) * log(x)**2
    end function getTestFuncOpenIntervalMATLAB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_Integration_mod