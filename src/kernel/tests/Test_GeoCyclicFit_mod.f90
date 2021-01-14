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
!!!!       https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief This module contains tests of the module [GeoCyclicFit_mod](@ref geocyclicfit_mod).
!>  \author Amir Shahmoradi

module Test_GeoCyclicFit_mod

    use Test_mod, only: Test_type
    use GeoCyclicFit_mod
    implicit none

    private
    public :: test_GeoCyclicFit

    type(Test_type) :: Test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_GeoCyclicFit()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_fitGeoCyclicLogPDF_1, "test_fitGeoCyclicLogPDF_1")
        call Test%finalize()
    end subroutine test_GeoCyclicFit

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_fitGeoCyclicLogPDF_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                     :: assertion
        integer(IK)                 :: i
        integer(IK) , parameter     :: nparam = 2_IK
        integer(IK) , parameter     :: numTrial = 64_IK
        integer(IK) , parameter     :: maxNumTrial = 64_IK
        integer(IK) , parameter     :: SuccessStep(numTrial) = [ (i, i = 1, numTrial) ]
        real(RK)    , parameter     :: LogCount(numTrial) = log( real(  [ 64_IK, 55_IK, 53_IK, 43_IK, 54_IK, 41_IK &
                                                                        , 45_IK, 55_IK, 50_IK, 42_IK, 48_IK, 52_IK &
                                                                        , 38_IK, 52_IK, 56_IK, 54_IK, 45_IK, 54_IK &
                                                                        , 69_IK, 50_IK, 50_IK, 49_IK, 45_IK, 38_IK &
                                                                        , 45_IK, 34_IK, 55_IK, 51_IK, 49_IK, 49_IK &
                                                                        , 47_IK, 58_IK, 37_IK, 54_IK, 50_IK, 59_IK &
                                                                        , 37_IK, 39_IK, 36_IK, 52_IK, 51_IK, 37_IK &
                                                                        , 44_IK, 46_IK, 37_IK, 29_IK, 41_IK, 39_IK &
                                                                        , 50_IK, 39_IK, 46_IK, 42_IK, 54_IK, 54_IK &
                                                                        , 54_IK, 24_IK, 44_IK, 43_IK, 37_IK, 43_IK &
                                                                        , 53_IK, 47_IK, 50_IK, 42_IK ], kind = RK ))
        real(RK)    , parameter     :: successProb = 0.7_RK
        real(RK)    , parameter     :: tolerance = 1.e-6_RK
        real(RK)                    :: xmin_ref(nparam) = [ 0.31952589641887075E-002_RK, 7.992349027030083_RK ]
        real(RK)                    :: Difference(nparam)
        type(GeoCyclicFit_type)     :: GeoCyclicFit

        GeoCyclicFit%PowellMinimum = GeoCyclicFit%fit(maxNumTrial, numTrial, SuccessStep, LogCount)
        assertion = .not. GeoCyclicFit%PowellMinimum%Err%occurred
        if (.not. assertion) return

        Difference = abs(GeoCyclicFit%PowellMinimum%xmin - xmin_ref) / abs(xmin_ref)
        assertion = all( Difference < tolerance )

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "xmin_ref    =", xmin_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "xmin        =", GeoCyclicFit%PowellMinimum%xmin
            write(Test%outputUnit,"(*(g0,:,' '))") "Difference  =", Difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_fitGeoCyclicLogPDF_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_GeoCyclicFit_mod