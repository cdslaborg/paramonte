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

!>  \brief This module contains tests of the module [pm_matrixInv](@ref pm_matrixInv).
!>  \author Amir Shahmoradi

module test_pm_matrixInv

    use pm_matrixTrans ! LCOV_EXCL_LINE
    use pm_err, only: err_type
    use pm_test, only: test_type, LK
    use pm_kind, only: LK

    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
    module function test_placeHolder() result(assertion); logical(LK) :: assertion; end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        implicit none

        test = test_type(MODULE_NAME)

        call test%run(test_genInvMat_1, SK_"test_genInvMat_1")
        call test%run(test_getInvChoLow_1, SK_"test_getInvChoLow_1")
        call test%run(test_getMatInvLow_1, SK_"test_getMatInvLow_1")
        call test%run(test_getMatInvDet_1, SK_"test_getMatInvDet_1")
       !call test%run(test_getMatInvChoUpp_1, SK_"test_getMatInvChoUpp_1")

        call test%run(test_getInvPosDefMat_1, SK_"test_getInvPosDefMat_1")
        call test%run(test_getInvPosDefMat_2, SK_"test_getInvPosDefMat_2")

        call test%run(test_getMatInvFromChoLow_1, SK_"test_getMatInvFromChoLow_1")
        call test%run(test_getInvLowFromChoLow_1, SK_"test_getInvLowFromChoLow_1")

        call test%run(test_getInvPosDefMatSqrtDet_1, SK_"test_getInvPosDefMatSqrtDet_1")
        call test%run(test_getInvPosDefMatSqrtDet_2, SK_"test_getInvPosDefMatSqrtDet_2")
        call test%run(test_getInvPosDefMatSqrtDet_3, SK_"test_getInvPosDefMatSqrtDet_3")

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getInvPosDefMatSqrtDet_1() result(assertion)

        use pm_kind, only: IK, RK
        implicit none

        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape(  [ 1._RK, 0._RK, 1._RK &
                                                                , 0._RK, 2._RK, 0._RK &
                                                                , 1._RK, 0._RK, 3._RK ], shape = shape(PosDefMat) )
        real(RK)    , parameter :: MatInvMat_ref(nd,nd) = reshape(  [ 1.500000000000000_RK, 0.000000000000000_RK, -0.50000000000000_RK &
                                                                    , 0.000000000000000_RK, 0.500000000000000_RK, 0.000000000000000_RK &
                                                                    , -0.50000000000000_RK, 0.000000000000000_RK, 0.500000000000000_RK ] &
                                                                    , shape = shape(MatInvMat_ref) )
        real(RK)    , parameter :: ChoDia_ref(nd) = [ 1.000000000000000_RK, 1.414213562373095_RK, 1.414213562373095_RK ]
        real(RK)    , parameter :: sqrtDetInvPosDefMat_ref = 0.5_RK
        real(RK)                :: MatInvMat(nd,nd), sqrtDetInvPosDefMat
        real(RK), allocatable   :: MatInvMat_diff(:,:), sqrtDetInvPosDefMat_diff

        MatInvMat = PosDefMat

        call setInvPosDefMatSqrtDet(nd = nd, MatInvMat = MatInvMat, sqrtDetInvPosDefMat = sqrtDetInvPosDefMat)

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(MatInvMat_diff)) deallocate(MatInvMat_diff); allocate(MatInvMat_diff, mold = MatInvMat)

        MatInvMat_diff = abs(MatInvMat - MatInvMat_ref)
        sqrtDetInvPosDefMat_diff = abs(sqrtDetInvPosDefMat - sqrtDetInvPosDefMat_ref)

        assertion = all(MatInvMat_diff < tolerance) .and. sqrtDetInvPosDefMat_diff < tolerance

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "MatInvMat_ref  = ", MatInvMat_ref
            write(test%disp%unit,"(*(g0,:,', '))") "MatInvMat      = ", MatInvMat
            write(test%disp%unit,"(*(g0,:,', '))") "MatInvMat_diff = ", MatInvMat_diff
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "sqrtDetInvPosDefMat_ref   = ", sqrtDetInvPosDefMat_ref
            write(test%disp%unit,"(*(g0,:,', '))") "sqrtDetInvPosDefMat_diff  = ", sqrtDetInvPosDefMat
            write(test%disp%unit,"(*(g0,:,', '))") "sqrtDetInvPosDefMat       = ", sqrtDetInvPosDefMat_diff
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getInvPosDefMatSqrtDet_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The output `sqrtDetInvPosDefMat` must be set to a negative value, if the input matrix is non-positive-definite.
    function test_getInvPosDefMatSqrtDet_2() result(assertion)
        use pm_kind, only: IK, RK
        implicit none
        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape(  [ 1._RK, 0._RK, -1._RK &
                                                                , 0._RK, 2._RK, -0._RK &
                                                                , 1._RK, 0._RK, -3._RK ], shape = shape(PosDefMat) )
        real(RK)                :: MatInvMat(nd,nd), sqrtDetInvPosDefMat

        MatInvMat = PosDefMat

        call setInvPosDefMatSqrtDet(nd = nd, MatInvMat = MatInvMat, sqrtDetInvPosDefMat = sqrtDetInvPosDefMat)

        assertion = sqrtDetInvPosDefMat < 0._RK

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "MatInvMat              = ", MatInvMat
            write(test%disp%unit,"(*(g0,:,', '))") "sqrtDetInvPosDefMat    = ", sqrtDetInvPosDefMat
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if

    end function test_getInvPosDefMatSqrtDet_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test with an 1-dimensional input matrix.
    function test_getInvPosDefMatSqrtDet_3() result(assertion)
        use pm_kind, only: IK, RK
        implicit none
        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 1_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: MatInvMat_ref(nd,nd) = reshape( [ 0.5_RK ], shape = shape(MatInvMat_ref) )
        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape( [ 2._RK ], shape = shape(PosDefMat) )
        real(RK)    , parameter :: sqrtDetInvPosDefMat_ref = 0.5_RK
        real(RK)                :: MatInvMat(nd,nd), sqrtDetInvPosDefMat
        real(RK), allocatable   :: MatInvMat_diff(:,:), sqrtDetInvPosDefMat_diff

        MatInvMat = PosDefMat

        call setInvPosDefMatSqrtDet(nd = nd, MatInvMat = MatInvMat, sqrtDetInvPosDefMat = sqrtDetInvPosDefMat)
        if (allocated(MatInvMat_diff)) deallocate(MatInvMat_diff); allocate(MatInvMat_diff, mold = MatInvMat)

        MatInvMat_diff = abs(MatInvMat - MatInvMat_ref)
        sqrtDetInvPosDefMat_diff = abs(sqrtDetInvPosDefMat - sqrtDetInvPosDefMat_ref)

        assertion = all(MatInvMat_diff < tolerance) .and. sqrtDetInvPosDefMat_diff < tolerance

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "MatInvMat_ref  = ", MatInvMat_ref
            write(test%disp%unit,"(*(g0,:,', '))") "MatInvMat      = ", MatInvMat
            write(test%disp%unit,"(*(g0,:,', '))") "MatInvMat_diff = ", MatInvMat_diff
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "sqrtDetInvPosDefMat_ref   = ", sqrtDetInvPosDefMat_ref
            write(test%disp%unit,"(*(g0,:,', '))") "sqrtDetInvPosDefMat_diff  = ", sqrtDetInvPosDefMat
            write(test%disp%unit,"(*(g0,:,', '))") "sqrtDetInvPosDefMat       = ", sqrtDetInvPosDefMat_diff
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP


    end function test_getInvPosDefMatSqrtDet_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getMatInvFromChoLow_1() result(assertion)

        use pm_kind, only: IK, RK
        implicit none

        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: choLow(nd,nd) = reshape( [ 1.000000000000000_RK, 0.000000000000000_RK, 1.000000000000000_RK &
                                                            , 0.000000000000000_RK, 2.000000000000000_RK, 0.000000000000000_RK &
                                                            , 1.000000000000000_RK, 0.000000000000000_RK, 3.000000000000000_RK ] &
                                                            , shape = shape(choLow) )
        real(RK)    , parameter :: choDia(nd) = [ 1.000000000000000_RK, 1.414213562373095_RK, 1.414213562373095_RK ]
        real(RK)    , parameter :: ChoDia_ref(nd) = [ 1.000000000000000_RK, 1.414213562373095_RK, 1.414213562373095_RK ]
        real(RK)    , parameter :: InvMatFromChoLow_ref(nd,nd) = reshape(  [ 1.500000000000000_RK, 0.000000000000000_RK, -0.50000000000000_RK &
                                                                            , 0.000000000000000_RK, 0.500000000000000_RK, 0.000000000000000_RK &
                                                                            , -0.50000000000000_RK, 0.000000000000000_RK, 0.500000000000000_RK ] &
                                                                            , shape = shape(InvMatFromChoLow_ref) )
        real(RK)    , parameter :: sqrtDetInvPosDefMat_ref = 0.5_RK
        real(RK)                :: InvMatFromChoLow(nd,nd)
        real(RK), allocatable   :: InvMatFromChoLow_diff(:,:)

        InvMatFromChoLow = getMatInvFromChoLow(nd = nd, choLow = choLow, choDia = choDia)

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(InvMatFromChoLow_diff)) deallocate(InvMatFromChoLow_diff); allocate(InvMatFromChoLow_diff, mold = InvMatFromChoLow)

        InvMatFromChoLow_diff = abs(InvMatFromChoLow - InvMatFromChoLow_ref)

        assertion = all(InvMatFromChoLow_diff < tolerance)

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "InvMatFromChoLow_ref  = ", InvMatFromChoLow_ref
            write(test%disp%unit,"(*(g0,:,', '))") "InvMatFromChoLow      = ", InvMatFromChoLow
            write(test%disp%unit,"(*(g0,:,', '))") "InvMatFromChoLow_diff = ", InvMatFromChoLow_diff
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getMatInvFromChoLow_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getInvLowFromChoLow_1() result(assertion)

        use pm_matrixSymCopy, only: setMatSymFromMatLow
        use pm_kind, only: IK, RK
        implicit none

        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: choLow(nd,nd) = reshape( [ 1.000000000000000_RK, 0.000000000000000_RK, 1.000000000000000_RK &
                                                            , 0.000000000000000_RK, 2.000000000000000_RK, 0.000000000000000_RK &
                                                            , 1.000000000000000_RK, 0.000000000000000_RK, 3.000000000000000_RK ] &
                                                            , shape = shape(choLow) )
        real(RK)    , parameter :: choDia(nd) = [ 1.000000000000000_RK, 1.414213562373095_RK, 1.414213562373095_RK ]
        real(RK)    , parameter :: ChoDia_ref(nd) = [ 1.000000000000000_RK, 1.414213562373095_RK, 1.414213562373095_RK ]
        real(RK)    , parameter :: InvMatFromChoLow_ref(nd,nd) = reshape(   [ 1.500000000000000_RK, 0.000000000000000_RK, -0.50000000000000_RK &
                                                                            , 0.000000000000000_RK, 0.500000000000000_RK, 0.000000000000000_RK &
                                                                            , -0.50000000000000_RK, 0.000000000000000_RK, 0.500000000000000_RK ] &
                                                                            , shape = shape(InvMatFromChoLow_ref) )
        real(RK)    , parameter :: sqrtDetInvPosDefMat_ref = 0.5_RK
        real(RK)                :: InvMatFromChoLow(nd,nd)
        real(RK), allocatable   :: InvMatFromChoLow_diff(:,:)

        InvMatFromChoLow = getInvLowFromChoLow(nd = nd, choLow = choLow, choDia = choDia)
        call setMatSymFromMatLow(InvMatFromChoLow)

        allocate(InvMatFromChoLow_diff, mold = InvMatFromChoLow)

        InvMatFromChoLow_diff = abs(InvMatFromChoLow - InvMatFromChoLow_ref)

        assertion = all(InvMatFromChoLow_diff < tolerance)

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "InvMatFromChoLow_ref  = ", InvMatFromChoLow_ref
            write(test%disp%unit,"(*(g0,:,', '))") "InvMatFromChoLow      = ", InvMatFromChoLow
            write(test%disp%unit,"(*(g0,:,', '))") "InvMatFromChoLow_diff = ", InvMatFromChoLow_diff
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        deallocate(InvMatFromChoLow_diff)

    end function test_getInvLowFromChoLow_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getInvPosDefMat_1() result(assertion)

        use pm_kind, only: IK, RK
        implicit none

        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape(  [ 1._RK, 0._RK, 1._RK &
                                                                , 0._RK, 2._RK, 0._RK &
                                                                , 1._RK, 0._RK, 3._RK ], shape = shape(PosDefMat) )
        real(RK)    , parameter :: MatInvMat_ref(nd,nd) = reshape(  [ 1.500000000000000_RK, 0.000000000000000_RK, -0.50000000000000_RK &
                                                                    , 0.000000000000000_RK, 0.500000000000000_RK, 0.000000000000000_RK &
                                                                    , -0.50000000000000_RK, 0.000000000000000_RK, 0.500000000000000_RK ] &
                                                                    , shape = shape(MatInvMat_ref) )
        real(RK)                :: MatInvMat(nd,nd)
        real(RK), allocatable   :: MatInvMat_diff(:,:)

        MatInvMat = getInvPosDefMat(nd = nd, PosDefMat = PosDefMat)

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(MatInvMat_diff)) deallocate(MatInvMat_diff); allocate(MatInvMat_diff, mold = MatInvMat)
        MatInvMat_diff = abs(MatInvMat - MatInvMat_ref)

        assertion = all(MatInvMat_diff < tolerance)

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "MatInvMat_ref  = ", MatInvMat_ref
            write(test%disp%unit,"(*(g0,:,', '))") "MatInvMat      = ", MatInvMat
            write(test%disp%unit,"(*(g0,:,', '))") "MatInvMat_diff = ", MatInvMat_diff
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if

    end function test_getInvPosDefMat_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The first element of `MatInvMat` must be set to a negative value, if the input matrix is non-positive-definite.
    function test_getInvPosDefMat_2() result(assertion)

        use pm_kind, only: IK, RK
        implicit none

        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape(  [ 1._RK, 0._RK, -1._RK &
                                                                , 0._RK, 2._RK, -0._RK &
                                                                , 1._RK, 0._RK, -3._RK ], shape = shape(PosDefMat) )
        real(RK)                :: MatInvMat(nd,nd)

        MatInvMat = getInvPosDefMat(nd = nd, PosDefMat = PosDefMat)

        assertion = MatInvMat(1,1) < 0._RK

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "MatInvMat      = ", MatInvMat
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getInvPosDefMat_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getMatInvDet_1() result(assertion)

        use pm_kind, only: IK, RK
        implicit none

        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape(  [ 1._RK, 0._RK, 1._RK &
                                                                , 0._RK, 2._RK, 0._RK &
                                                                , 1._RK, 0._RK, 3._RK ], shape = shape(PosDefMat) )
        real(RK)    , parameter :: MatInvMat_ref(nd,nd) = reshape(  [ 1.500000000000000_RK, 0.000000000000000_RK, -0.50000000000000_RK &
                                                                    , 0.000000000000000_RK, 0.500000000000000_RK, 0.000000000000000_RK &
                                                                    , -0.50000000000000_RK, 0.000000000000000_RK, 0.500000000000000_RK ] &
                                                                    , shape = shape(MatInvMat_ref) )
        real(RK)    , parameter :: MatrixLUP_ref(nd,nd)  = reshape( [ 1.000000000000000_RK, 0.000000000000000_RK, 1.000000000000000_RK &
                                                                    , 0.000000000000000_RK, 2.000000000000000_RK, 0.000000000000000_RK &
                                                                    , 1.000000000000000_RK, 0.000000000000000_RK, 2.000000000000000_RK ] &
                                                                    , shape = shape(MatInvMat_ref) )
        real(RK)    , parameter :: detInvMat_ref = 0.25_RK
        real(RK)                :: MatInvMat(nd,nd), detInvMat, detInvMat_diff
        real(RK), allocatable   :: MatrixLUP(:,:), MatInvMat_diff(:,:), MatrixLUP_diff(:,:)

        MatrixLUP = PosDefMat

        call getMatInvDet(nd = nd, MatrixLUP = MatrixLUP, InverseMatrix = MatInvMat, detInvMat = detInvMat, failed = assertion)
        assertion = .not. assertion
        call test%assert(assertion)

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(MatrixLUP_diff)) deallocate(MatrixLUP_diff); allocate(MatrixLUP_diff, mold = MatrixLUP)
        MatrixLUP_diff = abs(MatrixLUP - MatrixLUP_ref)

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(MatInvMat_diff)) deallocate(MatInvMat_diff); allocate(MatInvMat_diff, mold = MatInvMat)
        MatInvMat_diff = abs(MatInvMat - MatInvMat_ref)

        detInvMat_diff = abs(detInvMat - detInvMat_ref)

        assertion = all(MatInvMat_diff < tolerance) .and. detInvMat_diff < tolerance

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "MatrixLUP_ref  = ", MatrixLUP_ref
            write(test%disp%unit,"(*(g0,:,', '))") "MatrixLUP      = ", MatrixLUP
            write(test%disp%unit,"(*(g0,:,', '))") "MatrixLUP_diff = ", MatrixLUP_diff
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "MatInvMat_ref  = ", MatInvMat_ref
            write(test%disp%unit,"(*(g0,:,', '))") "MatInvMat      = ", MatInvMat
            write(test%disp%unit,"(*(g0,:,', '))") "MatInvMat_diff = ", MatInvMat_diff
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "detInvMat_ref  = ", detInvMat_ref
            write(test%disp%unit,"(*(g0,:,', '))") "detInvMat      = ", detInvMat
            write(test%disp%unit,"(*(g0,:,', '))") "detInvMat_diff = ", detInvMat_diff
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getMatInvDet_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_genInvMat_1() result(assertion)

        use pm_kind, only: IK, RK
        implicit none

        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape(  [ 1._RK, 0._RK, 1._RK &
                                                                , 0._RK, 2._RK, 0._RK &
                                                                , 1._RK, 0._RK, 3._RK ], shape = shape(PosDefMat) )
        real(RK)    , parameter :: InverseMatrix_ref(nd,nd) = reshape(  [ 1.500000000000000_RK, 0.000000000000000_RK, -0.50000000000000_RK &
                                                                        , 0.000000000000000_RK, 0.500000000000000_RK, 0.000000000000000_RK &
                                                                        , -0.50000000000000_RK, 0.000000000000000_RK, 0.500000000000000_RK ] &
                                                                        , shape = shape(InverseMatrix_ref) )
        real(RK)                :: InverseMatrix(nd,nd)
        real(RK), allocatable   :: InverseMatrix_diff(:,:)

        InverseMatrix = getMatInv(nd = nd, Matrix = PosDefMat)

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(InverseMatrix_diff)) deallocate(InverseMatrix_diff); allocate(InverseMatrix_diff, mold = InverseMatrix)

        InverseMatrix_diff = abs(InverseMatrix - InverseMatrix_ref)

        assertion = all(InverseMatrix_diff < tolerance)

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "InverseMatrix_ref  = ", InverseMatrix_ref
            write(test%disp%unit,"(*(g0,:,', '))") "InverseMatrix      = ", InverseMatrix
            write(test%disp%unit,"(*(g0,:,', '))") "InverseMatrix_diff = ", InverseMatrix_diff
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_genInvMat_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getMatInvLow_1() result(assertion)

        use pm_kind, only: IK, RK
        implicit none

        logical(LK)             :: assertion
        integer(IK)             :: i,j
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape(  [ 1.00_RK, +0.25_RK, +1.0_RK &
                                                                , 0.25_RK, +2.00_RK, -0.5_RK &
                                                                , 1.00_RK, -0.50_RK, +3.0_RK ], shape = shape(PosDefMat) )
        ! This is the Cholesky lower of the `PosDefMat`.
        real(RK)    , parameter :: MatLow(nd,nd) = transpose(reshape(   [ 1.00_RK,                 0._RK,                0._RK &
                                                                        , 0.25_RK, +1.391941090707505_RK,                0._RK &
                                                                        , 1.00_RK, -0.538815906080325_RK, 1.307546335452338_RK ], shape = shape(MatLow)))
        real(RK)    , parameter :: InvMatLow_ref(nd,nd) = transpose(reshape([ +1.000000000000000_RK,                0._RK,                0._RK &
                                                                            , -0.179605302026775_RK, 0.718421208107100_RK,                0._RK &
                                                                            , -0.838803309535462_RK, 0.296048226894869_RK, 0.764791252811745_RK ], shape = shape(InvMatLow_ref)))
        real(RK)                :: InvMatLow(nd,nd)
        real(RK), allocatable   :: InvMatLow_diff(:)

        InvMatLow = getLogPDF(nd, MatLow)

        InvMatLow_diff = abs([((InvMatLow(i,j) - InvMatLow_ref(i,j), i = j, nd), j = 1, nd)])
        assertion = all(InvMatLow_diff < tolerance)

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "InvMatLow_ref  = ", InvMatLow_ref
            write(test%disp%unit,"(*(g0,:,', '))") "InvMatLow      = ", InvMatLow
            write(test%disp%unit,"(*(g0,:,', '))") "InvMatLow_diff = ", InvMatLow_diff
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        deallocate(InvMatLow_diff)

    end function test_getMatInvLow_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getInvChoLow_1() result(assertion)

        use pm_kind, only: IK, RK
        implicit none

        logical(LK)             :: assertion
        integer(IK)             :: i,j
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape(  [ 1.00_RK, +0.25_RK, +1.0_RK &
                                                                , 0.25_RK, +2.00_RK, -0.5_RK &
                                                                , 1.00_RK, -0.50_RK, +3.0_RK ], shape = shape(PosDefMat) )
        ! This is the Cholesky lower of the `PosDefMat`.
        real(RK)    , parameter :: MatLow(nd,nd) = transpose(reshape(   [ 0.00_RK,                 0._RK, 0._RK &
                                                                        , 0.25_RK,                 0._RK, 0._RK &
                                                                        , 1.00_RK, -0.538815906080325_RK, 0._RK ], shape = shape(MatLow)))
        real(RK)    , parameter :: DiagMat(nd) = [1.00_RK, 1.391941090707505_RK, 1.307546335452338_RK]
        real(RK)    , parameter :: InvMatLow_ref(nd,nd) = transpose(reshape([ +1.000000000000000_RK,                0._RK,                0._RK &
                                                                            , -0.179605302026775_RK, 0.718421208107100_RK,                0._RK &
                                                                            , -0.838803309535462_RK, 0.296048226894869_RK, 0.764791252811745_RK ], shape = shape(InvMatLow_ref)))
        real(RK)                :: InvMatLow(nd,nd)
        real(RK), allocatable   :: InvMatLow_diff(:)

        InvMatLow = getLogPDF(nd, MatLow, DiagMat)

        InvMatLow_diff = abs([((InvMatLow(i,j) - InvMatLow_ref(i,j), i = j, nd), j = 1, nd)])
        assertion = all(InvMatLow_diff < tolerance)

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "InvMatLow_ref  = ", InvMatLow_ref
            write(test%disp%unit,"(*(g0,:,', '))") "InvMatLow      = ", InvMatLow
            write(test%disp%unit,"(*(g0,:,', '))") "InvMatLow_diff = ", InvMatLow_diff
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        deallocate(InvMatLow_diff)

    end function test_getInvChoLow_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    function test_getMatInvChoUpp_1() result(assertion)
!
!        use pm_kind, only: IK, RK
!        implicit none
!
!        logical(LK)             :: assertion
!        integer(IK) , parameter :: nd = 3_IK
!        real(RK)    , parameter :: tolerance = 1.e-12_RK
!        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape(  [ +2.0_RK, 0.5_RK, -0.3_RK &
!                                                                , +0.5_RK, 1.0_RK, +0.1_RK &
!                                                                , -0.3_RK, 0.1_RK, +5.0_RK &
!                                                                ] , shape = shape(PosDefMat) )
!        real(RK)    , parameter :: InvMatChoUpp_ref(nd,nd) = transpose( reshape([ +0.579558652729384_RK, -0.385983440243132_RK, +0.053396918610710_RK &
!                                                                                , -0.293844367015099_RK, +1.150987224157956_RK, -0.020020030050088_RK &
!                                                                                , +0.040650406504065_RK, -0.040650406504065_RK, +0.203252032520325_RK &
!                                                                                ], shape = shape(InvMatChoUpp_ref) ))
!        real(RK)    , parameter :: InvMatChoDia_ref(nd) = [ 0.761287496764123_RK, 1.001001502504383_RK, 0.447213595499958_RK]
!        real(RK)                :: InvMatChoUpp(nd,nd)
!        real(RK)                :: InvMatChoDia(nd)
!        real(RK), allocatable   :: InvMatChoDia_diff(:)
!        real(RK), allocatable   :: InvMatChoUpp_diff(:,:)
!        real(RK)                :: choLow(nd,nd)
!        real(RK)                :: choDia(nd)
!
!        choLow = PosDefMat
!        call setChoLow(choLow, choDia, nd)
!        call getMatInvChoUpp(nd,choLow,choDia,InvMatChoUpp,InvMatChoDia)
!
!        InvMatChoUpp_diff = abs(InvMatChoUpp - InvMatChoUpp_ref)
!        assertion = all(InvMatChoUpp_diff < tolerance)
!        call test%assert(assertion)
!
!        InvMatChoDia_diff = abs(InvMatChoDia - InvMatChoDia_ref)
!        assertion = assertion .and. all(InvMatChoDia_diff < tolerance)
!        call test%assert(assertion)
!
!        if (test%traceable .and. .not. assertion) then
!            ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0,:,', '))")
!            write(test%disp%unit,"(*(g0,:,', '))") "InvMatChoUpp_ref   = ", InvMatChoUpp_ref
!            write(test%disp%unit,"(*(g0,:,', '))") "InvMatChoUpp       = ", InvMatChoUpp
!            write(test%disp%unit,"(*(g0,:,', '))") "InvMatChoUpp_diff  = ", InvMatChoUpp_diff
!            write(test%disp%unit,"(*(g0,:,', '))") "InvMatChoDia_ref   = ", InvMatChoDia_ref
!            write(test%disp%unit,"(*(g0,:,', '))") "InvMatChoDia       = ", InvMatChoDia
!            write(test%disp%unit,"(*(g0,:,', '))") "InvMatChoDia_diff  = ", InvMatChoDia_diff
!            write(test%disp%unit,"(*(g0,:,', '))")
!            ! LCOV_EXCL_STOP
!        end if
!
!        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
!        deallocate(InvMatChoUpp_diff)
!
!    end function test_getMatInvChoUpp_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_matrixInv ! LCOV_EXCL_LINE