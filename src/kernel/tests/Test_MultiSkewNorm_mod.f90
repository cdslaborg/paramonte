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

!>  \brief This module contains tests of the module [MultiSkewNorm_mod](@ref multiskewnorm_mod).
!>  \author Amir Shahmoradi

module Test_MultiSkewNorm_mod

    use MultiSkewNorm_mod
    use Err_mod, only: Err_type
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_MultiSkewNorm

    type(Test_type) :: Test

    integer(IK) , parameter :: lenRnd = 50_IK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_MultiSkewNorm()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_getRandMSN_1, "test_getRandMSN_1")
        call Test%run(test_getAugChoFac_1, "test_getAugChoFac_1")
        call Test%run(test_getLogProbMSN_1, "test_getLogProbMSN_1")
        call Test%run(test_getLogProbMSN_2, "test_getLogProbMSN_2")
        call Test%finalize()
    end subroutine test_MultiSkewNorm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! \todo
    ! What is the best method of testing for randomness?
    function test_getRandMSN_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK, np = 10000_IK
        real(RK)    , parameter :: Alpha(nd) = [-10._RK, 10._RK]
        real(RK)    , parameter :: MeanVec(nd) = [ 1._RK, 2._RK ]
        real(RK)    , parameter :: CovMat(nd,nd) = reshape( [ 1.5_RK, -0.9_RK, -0.9_RK, 1.5_RK ], shape = shape(CovMat) )
        real(RK)                :: AugChoLow(nd+1,nd+1)
        real(RK)                :: AugChoDia(nd+1)
        integer(IK)             :: ip

        assertion = .true. !< There is really no easy testing of randomness. For now, we trust the giants.

        call getAugChoFac(nd = nd, Alpha = Alpha, CovMat = CovMat, AugChoLow = AugChoLow, AugChoDia = AugChoDia)

        Test%File = Test%openFile()
        do ip = 1, np
            write(Test%File%unit,"(*(g0,:,','))") getRandMSN(nd = nd, MeanVec = MeanVec, AugChoLow = AugChoLow, AugChoDia = AugChoDia)
        end do
        close(Test%File%unit)

    end function test_getRandMSN_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbMSN_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        integer(IK)             :: i
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        real(RK)    , parameter :: MeanVec(nd) = [(real(i**2+1._RK,RK),i=1,nd)]
        real(RK)    , parameter :: Alpha(nd) = [-2._RK,5._RK,2._RK]
        real(RK)    , parameter :: Point(nd) = MeanVec + [(real(i,RK),i=1,nd)]
        real(RK)    , parameter :: InvCovMat(nd,nd) = reshape(  [ 1.500000000000000_RK, 0.000000000000000_RK, -0.50000000000000_RK &
                                                                , 0.000000000000000_RK, 0.500000000000000_RK, 0.000000000000000_RK &
                                                                , -0.50000000000000_RK, 0.000000000000000_RK, 0.500000000000000_RK ] &
                                                                , shape = shape(InvCovMat) )
        real(RK)    , parameter :: logSqrtDetInvCovMat = -0.693147180559945_RK
        real(RK)    , parameter :: logProb_ref = -5.2568155996140193_RK
        real(RK)                :: difference
        real(RK)                :: logProb

        logProb = getLogProbMSN ( nd = nd & ! LCOV_EXCL_LINE
                                , MeanVec = MeanVec & ! LCOV_EXCL_LINE
                                , InvCovMat = InvCovMat & ! LCOV_EXCL_LINE
                                , logSqrtDetInvCovMat = logSqrtDetInvCovMat & ! LCOV_EXCL_LINE
                                , Alpha = Alpha & ! LCOV_EXCL_LINE
                                , point = point & ! LCOV_EXCL_LINE
                                )

        difference = abs( (logProb - logProb_ref) / logProb_ref )
        assertion = difference <= tolerance

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "logProb_ref    ", logProb_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "logProb        ", logProb
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference     ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbMSN_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> When `alpha = 0.`, `getLogProbMSN()` must be equivalent to `getLogProbMVN()`. 
    function test_getLogProbMSN_2() result(assertion)

        use Statistics_mod, only: getLogProbMVN
        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        integer(IK)             :: i
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        real(RK)    , parameter :: MeanVec(nd) = [(real(i**2+1._RK,RK),i=1,nd)]
        real(RK)    , parameter :: Alpha(nd) = [0._RK, 0._RK, 0._RK]
        real(RK)    , parameter :: Point(nd) = MeanVec + [(real(i,RK),i=1,nd)]
        real(RK)    , parameter :: InvCovMat(nd,nd) = reshape(  [ 1.500000000000000_RK, 0.000000000000000_RK, -0.50000000000000_RK &
                                                                , 0.000000000000000_RK, 0.500000000000000_RK, 0.000000000000000_RK &
                                                                , -0.50000000000000_RK, 0.000000000000000_RK, 0.500000000000000_RK ] &
                                                                , shape = shape(InvCovMat) )
        real(RK)    , parameter :: logSqrtDetInvCovMat = -0.693147180559945_RK
        real(RK)                :: logProb_ref
        real(RK)                :: difference
        real(RK)                :: logProb

        logProb_ref = getLogProbMVN ( nd = nd &
                                    , MeanVec = MeanVec &
                                    , InvCovMat = InvCovMat &
                                    , logSqrtDetInvCovMat = logSqrtDetInvCovMat &
                                    , point = point &
                                    )
        logProb = getLogProbMSN ( nd = nd & ! LCOV_EXCL_LINE
                                , MeanVec = MeanVec & ! LCOV_EXCL_LINE
                                , InvCovMat = InvCovMat & ! LCOV_EXCL_LINE
                                , logSqrtDetInvCovMat = logSqrtDetInvCovMat & ! LCOV_EXCL_LINE
                                , Alpha = Alpha & ! LCOV_EXCL_LINE
                                , point = point & ! LCOV_EXCL_LINE
                                )

        difference = abs( (logProb - logProb_ref) / logProb_ref )
        assertion = difference <= tolerance

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "logProb_ref    ", logProb_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "logProb        ", logProb
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference     ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbMSN_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getAugChoFac_1() result(assertion)
        implicit none
        logical :: assertion
        integer(IK)             :: i, j
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-6_RK
        real(RK)    , parameter :: Alpha(nd) = [-2._RK, 1._RK]
        real(RK)    , parameter :: CovMat(nd,nd) = reshape( [ 1.5_RK, -0.9_RK, -0.9_RK, 1.5_RK ], shape = shape(CovMat) )
        real(RK)    , parameter :: AugCovMat_ref(nd,nd) = reshape( [ +1.000000000000000_RK, -1.121171170000000_RK,  0.948683300000000_RK &
                                                                   , -1.121171170000000_RK, +1.500000000000000_RK, -0.900000000000000_RK &
                                                                   , +0.948683300000000_RK, -0.900000000000000_RK,  1.500000000000000_RK &
                                                                   ], shape = shape(AugCovMat_ref) )
        real(RK)    , parameter :: AugChoLow_ref(nd+1,nd+1) = reshape( [ +1.000000000000000_RK, -1.121171170000000_RK, +0.948683300000000_RK &
                                                                       ,                 0._RK, +0.492925154116557_RK, +0.331970004074427_RK &
                                                                       ,                 0._RK,                 0._RK, +0.699854208171913_RK &
                                                                       ], shape = shape(AugChoLow_ref) )
        real(RK)                :: AugChoLow(nd+1,nd+1), AugChoDia(nd+1), difference

        assertion = .true.

        call getAugChoFac(nd = nd, Alpha = Alpha, CovMat = CovMat, AugChoLow = AugChoLow, AugChoDia = AugChoDia)

        do j = 1, nd + 1
            AugChoLow(j,j) = AugChoDia(j)
            do i = j, nd
                difference = abs( (AugChoLow(i,j) - AugChoLow_ref(i,j)) / AugChoLow_ref(i,j))
                assertion = assertion .and. difference <= tolerance
                if (Test%isVerboseMode .and. .not. assertion) then
                    ! LCOV_EXCL_START
                    write(Test%outputUnit,"(*(g0,:,', '))")
                    write(Test%outputUnit,"(*(g0,:,', '))") "i, j               ", i, j
                    write(Test%outputUnit,"(*(g0,:,', '))") "AugChoLow_ref(i,j) ", AugChoLow_ref(i,j)
                    write(Test%outputUnit,"(*(g0,:,', '))") "AugChoLow(i,j)     ", AugChoLow(i,j)
                    write(Test%outputUnit,"(*(g0,:,', '))") "Difference         ", difference
                    write(Test%outputUnit,"(*(g0,:,', '))")
                    ! LCOV_EXCL_STOP
                end if
            end do
        end do


    end function test_getAugChoFac_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_MultiSkewNorm_mod ! LCOV_EXCL_LINE