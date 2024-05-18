!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief 
!>  This module contains tests of the module [pm_distGeomCyclic](@ref pm_distGeomCyclic).
!>
!>  \final
!>
!>  \author 
!>  Amir Shahmoradi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_distGeomCyclic

    use pm_test, only: test_type, LK
    use pm_distGeomCyclic
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)
        call test%run(test_fitGeoCyclicLogPDF_1, SK_"test_fitGeoCyclicLogPDF_1")
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_fitGeoCyclicLogPDF_1() result(assertion)
        use pm_kind, only: IK, RK
        implicit none
        logical(LK)                 :: assertion
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
        type(FitGeoCyclic_type)     :: FitGeoCyclic

        FitGeoCyclic%powellMin = FitGeoCyclic%fit(maxNumTrial, numTrial, SuccessStep, LogCount)
        assertion = .not. FitGeoCyclic%powellMin%err%occurred
        call test%assert(assertion)

        Difference = abs(FitGeoCyclic%powellMin%xmin - xmin_ref) / abs(xmin_ref)
        assertion = all( Difference < tolerance )
        call test%assert(assertion)

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))") "xmin_ref    =", xmin_ref
            write(test%disp%unit,"(*(g0,:,' '))") "xmin        =", FitGeoCyclic%powellMin%xmin
            write(test%disp%unit,"(*(g0,:,' '))") "Difference  =", Difference
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_fitGeoCyclicLogPDF_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_distGeomCyclic