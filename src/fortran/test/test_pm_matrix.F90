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
!>  This module contains tests of the module [pm_matrixDet](@ref pm_matrixDet).
!>
!>  \final
!>
!>  \author 
!>  Amir Shahmoradi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_matrix

    use pm_matrixDet ! LCOV_EXCL_LINE
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

        call test%run(test_getDet_1, SK_"test_getDet_1")

        call test%run(test_isPosDef_1, SK_"test_isPosDef_1")
        call test%run(test_isPosDef_2, SK_"test_isPosDef_2")

        call test%run(test_sortPosDefMat_1, SK_"test_sortPosDefMat_1")
        call test%run(test_getRegresCoef_1, SK_"test_getRegresCoef_1")

        call test%run(test_getMatDetSqrtPosDefMat_1, SK_"test_getMatDetSqrtPosDefMat_1")
        call test%run(test_getMatDetSqrtPosDefMat_2, SK_"test_getMatDetSqrtPosDefMat_2")

        call test%run(test_getMatDetSqrtLogPosDefMat_1, SK_"test_getMatDetSqrtLogPosDefMat_1")
        call test%run(test_getMatDetSqrtLogPosDefMat_2, SK_"test_getMatDetSqrtLogPosDefMat_2")

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getDet_1() result(assertion)

        use pm_kind, only: IK, RK
        implicit none

        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: mat(nd,nd) = reshape([ -1._RK, -0._RK, -1._RK &
                                                        , -0._RK, -2._RK, -0._RK &
                                                        , -1._RK, -0._RK, -3._RK ], shape = shape(mat) )
        real(RK)    , parameter :: determinant_ref = -4._RK
        real(RK)                :: determinant, determinant_diff

        determinant = getMatDet(Matrix = mat)

        determinant_diff = abs(determinant - determinant_ref)

        assertion = determinant_diff < tolerance

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "determinant_ref  = ", determinant_ref
            write(test%disp%unit,"(*(g0,:,', '))") "determinant      = ", determinant
            write(test%disp%unit,"(*(g0,:,', '))") "determinant_diff = ", determinant_diff
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getDet_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getMatDetSqrtPosDefMat_1() result(assertion)

        use pm_kind, only: IK, RK
        implicit none

        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: mat(nd,nd) = reshape(  [ 1._RK, 0._RK, 1._RK &
                                                                , 0._RK, 2._RK, 0._RK &
                                                                , 1._RK, 0._RK, 3._RK ], shape = shape(mat) )
        real(RK)    , parameter :: sqrtDetPosDefMat_ref = 2._RK
        real(RK)                :: sqrtDetPosDefMat, sqrtDetPosDefMat_diff

        sqrtDetPosDefMat = getMatDetSqrtPosDefMat(nd = nd, mat = mat)

        sqrtDetPosDefMat_diff = abs(sqrtDetPosDefMat - sqrtDetPosDefMat_ref)

        assertion = sqrtDetPosDefMat_diff < tolerance

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "sqrtDetPosDefMat_ref  = ", sqrtDetPosDefMat_ref
            write(test%disp%unit,"(*(g0,:,', '))") "sqrtDetPosDefMat      = ", sqrtDetPosDefMat
            write(test%disp%unit,"(*(g0,:,', '))") "sqrtDetPosDefMat_diff = ", sqrtDetPosDefMat_diff
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getMatDetSqrtPosDefMat_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The output `sqrtDetPosDefMat` must be set to a negative value, if the input matrix is non-positive-definite.
    function test_getMatDetSqrtPosDefMat_2() result(assertion)

        use pm_kind, only: IK, RK
        implicit none

        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: mat(nd,nd) = reshape(  [ -1._RK, -0._RK, -1._RK &
                                                                , -0._RK, -2._RK, -0._RK &
                                                                , -1._RK, -0._RK, -3._RK ], shape = shape(mat) )
        real(RK)                :: sqrtDetPosDefMat

        sqrtDetPosDefMat = getMatDetSqrtPosDefMat(nd = nd, mat = mat)

        assertion = sqrtDetPosDefMat < 0._RK

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "sqrtDetPosDefMat      = ", sqrtDetPosDefMat
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if

    end function test_getMatDetSqrtPosDefMat_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getMatDetSqrtLogPosDefMat_1() result(assertion)

        use pm_kind, only: IK, RK
        implicit none

        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: mat(nd,nd) = reshape(  [ 1._RK, 0._RK, 1._RK &
                                                                , 0._RK, 2._RK, 0._RK &
                                                                , 1._RK, 0._RK, 3._RK ], shape = shape(mat) )
        real(RK)    , parameter :: logSqrtDetPosDefMat_ref = log(2._RK)
        real(RK)                :: logSqrtDetPosDefMat, logSqrtDetPosDefMat_diff
        real(RK), allocatable   :: Matrix(:,:)
        logical(LK)             :: failed

        Matrix = mat

        call getMatDetSqrtLogPosDefMat(nd = nd, mat = Matrix, logSqrtDetPosDefMat = logSqrtDetPosDefMat, failed = failed)

        assertion = .not. failed
        call test%assert(assertion)

        logSqrtDetPosDefMat_diff = abs(logSqrtDetPosDefMat - logSqrtDetPosDefMat_ref)

        assertion = logSqrtDetPosDefMat_diff < tolerance

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "logSqrtDetPosDefMat_ref  = ", logSqrtDetPosDefMat_ref
            write(test%disp%unit,"(*(g0,:,', '))") "logSqrtDetPosDefMat      = ", logSqrtDetPosDefMat
            write(test%disp%unit,"(*(g0,:,', '))") "logSqrtDetPosDefMat_diff = ", logSqrtDetPosDefMat_diff
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getMatDetSqrtLogPosDefMat_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getMatDetSqrtLogPosDefMat_2() result(assertion)

        use pm_kind, only: IK, RK
        implicit none

        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: mat(nd,nd) = reshape(  [ 1._RK, 0._RK, 1._RK &
                                                                , 0._RK, 2._RK, 0._RK &
                                                                , 1._RK, 0._RK, 3._RK ], shape = shape(mat) )
        real(RK)    , parameter :: logSqrtDetPosDefMat_ref = log(2._RK)
        real(RK)                :: logSqrtDetPosDefMat
        real(RK), allocatable   :: Matrix(:,:)
        logical(LK)             :: failed

        Matrix = -mat

        call getMatDetSqrtLogPosDefMat(nd = nd, mat = Matrix, logSqrtDetPosDefMat = logSqrtDetPosDefMat, failed = failed)

        assertion = failed
        if (.not. assertion) return

    end function test_getMatDetSqrtLogPosDefMat_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isPosDef_1() result(assertion)
        use pm_kind, only: IK, RK
        implicit none
        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: mat(nd,nd) = reshape(  [ 1._RK, 0._RK, 1._RK &
                                                                , 0._RK, 2._RK, 0._RK &
                                                                , 1._RK, 0._RK, 3._RK ], shape = shape(mat) )

        assertion = isPosDef(nd = nd, Matrix = mat)
    end function test_isPosDef_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isPosDef_2() result(assertion)
        use pm_kind, only: IK, RK
        implicit none
        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: mat(nd,nd) = reshape(  [ -1._RK, -0._RK, -1._RK &
                                                                , -0._RK, -2._RK, -0._RK &
                                                                , -1._RK, -0._RK, -3._RK ], shape = shape(mat) )

        assertion = .not. isPosDef(nd = nd, Matrix = mat)
    end function test_isPosDef_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_sortPosDefMat_1() result(assertion)

        use pm_kind, only: RK, IK

        implicit none
        logical(LK)             :: assertion
        integer(IK), parameter  :: ColIndx(*) = [4]
        integer(IK), parameter  :: ColIndxMap(*) = [5]
        integer(IK)             :: rank
        real(RK), allocatable   :: mat(:,:), OutPosDefMat(:,:), OutPosDefMat_ref(:,:)
        integer(IK)             :: i, j

        assertion = .true._LK

        rank = 5
        allocate(mat(rank,rank))
        do j = 1,rank
            do i = 1,j
                mat(i,j) = i*10 + j
            end do
        end do

        ! switch the variable 4 with 5, such that the output remains a positive-definite matrix.

        OutPosDefMat = getSortedPosDefMat(rank, mat, 1_IK, ColIndx, ColIndxMap)

        OutPosDefMat_ref = reshape( [ 11._RK, 12._RK, 13._RK, 15._RK, 14._RK &
                                    ,  0._RK, 22._RK, 23._RK, 25._RK, 24._RK &
                                    ,  0._RK,  0._RK, 33._RK, 35._RK, 34._RK &
                                    ,  0._RK,  0._RK,  0._RK, 55._RK, 45._RK &
                                    ,  0._RK,  0._RK,  0._RK,  0._RK, 44._RK &
                                    ], shape = [rank,rank])
        OutPosDefMat_ref = transpose(OutPosDefMat_ref)
        do j = 1, rank
            do i = 1, j
                assertion = assertion .and. OutPosDefMat_ref(i,j) == OutPosDefMat(i,j)
            end do
        end do

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START

            write(test%disp%unit,"(*(g0))")
            write(test%disp%unit,"(*(g0))") "OutPosDefMat_ref:"
            do i = 1,rank
                write(test%disp%unit,"(*(F7.1))") (OutPosDefMat_ref(i,j),j=1,rank)
            end do
            write(test%disp%unit,"(*(g0))")

            write(test%disp%unit,"(*(g0))")
            write(test%disp%unit,"(*(g0))") "Output Positive-Definite Matrix with variables 4 and 5 swapped:"
            do i = 1,rank
                write(test%disp%unit,"(*(F7.1))") (OutPosDefMat(i,j),j=1,rank)
            end do
            write(test%disp%unit,"(*(g0))")

            write(test%disp%unit,"(*(g0))")
            write(test%disp%unit,"(*(g0))") "Original Matrix:"
            do i = 1,rank
                write(test%disp%unit,"(*(F7.1))") (mat(i,j),j=1,rank)
            end do

        end if
        ! LCOV_EXCL_STOP

    end function test_sortPosDefMat_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getRegresCoef_1() result(assertion)

        use pm_kind, only: RK, IK

        implicit none
        logical(LK)             ::  assertion
        real(RK)                ::  normalizedDifference
        integer , parameter     ::  rankPDM = 4, rankS11 = 3, rankS22 = 1
        real(RK), parameter     ::  mat(rankPDM,rankPDM) = reshape( &
                                    [ 4.414182620515998_RK, 1.173760167060120_RK, 0.757607629189287_RK, 5.075277296976230_RK &
                                    , 1.173760167060120_RK, 0.866956750091570_RK, 0.310654936099342_RK, 1.621274787164182_RK &
                                    , 0.757607629189287_RK, 0.310654936099342_RK, 0.955157221699132_RK, 1.254186231887444_RK &
                                    , 5.075277296976230_RK, 1.621274787164182_RK, 1.254186231887444_RK, 6.407791808961157_RK ] &
                                    , shape=shape(mat) )
        real(RK), parameter     ::  RegresCoefMatRef(rankS11,rankS22) = reshape( &
                                    [ 0.792047783119072 &
                                    , 0.253016145889269 &
                                    , 0.195728305363088 ] &
                                    , shape=shape(RegresCoefMatRef) )
        real(RK), parameter     ::  SchurComplementRef(rankS11,rankS11) = reshape( &
                                    [ 0.3943204887314210_RK, -0.110366933940115_RK, -0.235767795395625_RK &
                                    , -0.110366933940115_RK,  0.456748052015843_RK, -0.006674430520204_RK &
                                    , -0.235767795395625_RK, -0.006674430520204_RK,  0.709677475922085_RK ] &
                                    , shape=shape(SchurComplementRef) )
        real(RK)                ::  SchurComplement(rankS11,rankS11)
        real(RK)                ::  RegresCoefMat(rankS11,rankS22)
        integer(IK)             ::  i,j

        assertion = .true._LK

        call getRegresCoef  ( rankPDM           = rankPDM           &
                            , rankS11           = rankS11           &
                            , rankS22           = rankS22           &
                            , mat         = mat         &
                            , RegresCoefMat     = RegresCoefMat     &
                            , SchurComplement   = SchurComplement   &
                            , failed            = assertion         &
                            )
        assertion = .not. assertion
        call test%assert(assertion)

        do i = 1,rankS11
            !write(test%disp%unit,"(*(F22.15))") (RegresCoefMat(i,j),j=1,rankS22), (RegresCoefMatRef(i,j),j=1,rankS22)
            do j = 1,rankS22
                normalizedDifference = abs(RegresCoefMatRef(i,j) - RegresCoefMat(i,j)) / (RegresCoefMatRef(i,j) + RegresCoefMat(i,j))
                assertion = assertion .and. normalizedDifference < 1.e-5_RK
            end do
        end do

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0))")
            write(test%disp%unit,"(*(g0))") "Original Covariance Matrix:"
            do i = 1,rankPDM
                write(test%disp%unit,"(*(F22.15))") (mat(i,j),j=1,rankPDM)
            end do
            write(test%disp%unit,"(*(g0))")
            write(test%disp%unit,"(*(g0))") "RegresCoefMat, RegresCoefMatRef, difference, normalizedDifference:"
            write(test%disp%unit,"(*(g0))")
            write(test%disp%unit,"(*(g0))") "SchurComplement, SchurComplementRef, difference, normalizedDifference:"
            write(test%disp%unit,"(*(g0))")
            do i = 1,rankS11
                do j = 1,rankS22
                    write(test%disp%unit,"(*(F22.15))") RegresCoefMat(i,j) &
                                                        , RegresCoefMatRef(i,j) &
                                                        , abs(RegresCoefMat(i,j)-RegresCoefMatRef(i,j)) &
                                                        , normalizedDifference
                end do
            end do
            write(test%disp%unit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
        call test%assert(assertion)

        do i = 1,rankS11
            do j = 1,rankS11
                normalizedDifference = abs(SchurComplementRef(i,j) - SchurComplement(i,j)) / (SchurComplementRef(i,j) + SchurComplement(i,j))
                assertion = assertion .and. normalizedDifference < 1.e-6_RK
            end do
        end do

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            do i = 1,rankS11
                do j = 1,rankS11
                    write(test%disp%unit,"(*(F22.15))") SchurComplement(i,j) &
                                                        , SchurComplementRef(i,j) &
                                                        , abs(SchurComplement(i,j)-SchurComplementRef(i,j)) &
                                                        , normalizedDifference
                end do
            end do
        end if
        ! LCOV_EXCL_STOP

    end function test_getRegresCoef_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_matrix ! LCOV_EXCL_LINE