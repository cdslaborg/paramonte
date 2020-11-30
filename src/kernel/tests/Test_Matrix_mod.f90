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

!>  \brief This module contains tests of the module [Matrix_mod](@ref matrix_mod).
!>  @author Amir Shahmoradi

module Test_Matrix_mod

    use Matrix_mod
    use Err_mod, only: Err_type
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_Matrix

    type(Test_type) :: Test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_Matrix()

        implicit none

        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_isPosDef_1, "test_isPosDef_1")
        call Test%run(test_isPosDef_2, "test_isPosDef_2")
        call Test%run(test_getInvMat_1, "test_getInvMat_1")
        call Test%run(test_getOuterProd_1, "test_getOuterProd_1")
        call Test%run(test_getInvMatDet_1, "test_getInvMatDet_1")
        call Test%run(test_sortPosDefMat_1, "test_sortPosDefMat_1")
        call Test%run(test_getRegresCoef_1, "test_getRegresCoef_1")
        call Test%run(test_multiplyMatrix_1, "test_multiplyMatrix_1")
        call Test%run(test_getDeterminant_1, "test_getDeterminant_1")
        call Test%run(test_getInvPosDefMat_1, "test_getInvPosDefMat_1")
        call Test%run(test_getInvPosDefMat_2, "test_getInvPosDefMat_2")
        call Test%run(test_getCholeskyFactor_1, "test_getCholeskyFactor_1")
        call Test%run(test_getCholeskyFactor_2, "test_getCholeskyFactor_2")
        call Test%run(test_getInvMatFromCholFac, "test_getInvMatFromCholFac")
        call Test%run(test_getSqrtDetPosDefMat_1, "test_getSqrtDetPosDefMat_1")
        call Test%run(test_getSqrtDetPosDefMat_2, "test_getSqrtDetPosDefMat_2")
        call Test%run(test_getInvPosDefMatSqrtDet_1, "test_getInvPosDefMatSqrtDet_1")
        call Test%run(test_getInvPosDefMatSqrtDet_2, "test_getInvPosDefMatSqrtDet_2")
        call Test%run(test_getLogSqrtDetPosDefMat_1, "test_getLogSqrtDetPosDefMat_1")
        call Test%run(test_getLogSqrtDetPosDefMat_2, "test_getLogSqrtDetPosDefMat_2")
        call Test%run(test_symmetrizeUpperSquareMatrix_1, "test_symmetrizeUpperSquareMatrix_1")
        call Test%finalize()

    end subroutine test_Matrix

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getCholeskyFactor_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape(  [ 1._RK, 0._RK, 1._RK &
                                                                , 0._RK, 2._RK, 0._RK &
                                                                , 1._RK, 0._RK, 3._RK ], shape = shape(PosDefMat) )
        real(RK)    , parameter :: CholeskyLower_ref(nd,nd) = reshape(  [ 1.000000000000000_RK, 0.000000000000000_RK, 1.000000000000000_RK &
                                                                        , 0.000000000000000_RK, 2.000000000000000_RK, 0.000000000000000_RK &
                                                                        , 1.000000000000000_RK, 0.000000000000000_RK, 3.000000000000000_RK ] &
                                                                        , shape = shape(CholeskyLower_ref) )
        real(RK)    , parameter :: CholeskyDiagonal_ref(nd) = [ 1.000000000000000_RK, 1.414213562373095_RK, 1.414213562373095_RK ]
        real(RK)                :: CholeskyLower(nd,nd), CholeskyDiagonal(nd)
        real(RK), allocatable   :: CholeskyLower_diff(:,:), CholeskyDiagonal_diff(:)
        CholeskyLower = PosDefMat
        call getCholeskyFactor(nd = nd, PosDefMat = CholeskyLower, Diagonal = CholeskyDiagonal)

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(CholeskyLower_diff)) deallocate(CholeskyLower_diff); allocate(CholeskyLower_diff, mold = PosDefMat)
        CholeskyLower_diff = abs(PosDefMat - CholeskyLower_ref)

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(CholeskyDiagonal_diff)) deallocate(CholeskyDiagonal_diff); allocate(CholeskyDiagonal_diff, mold = CholeskyDiagonal)
        CholeskyDiagonal_diff = abs(CholeskyDiagonal - CholeskyDiagonal_ref)
        assertion = all(CholeskyLower_diff < tolerance) .and. all(CholeskyDiagonal_diff < tolerance)

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "CholeskyLower_ref  = ", CholeskyLower_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "CholeskyLower      = ", CholeskyLower
            write(Test%outputUnit,"(*(g0,:,', '))") "CholeskyLower_diff = ", CholeskyLower_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "CholeskyDiagonal_ref   = ", CholeskyDiagonal_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "CholeskyDiagonal       = ", CholeskyDiagonal
            write(Test%outputUnit,"(*(g0,:,', '))") "CholeskyDiagonal_diff  = ", CholeskyDiagonal_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getCholeskyFactor_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getCholeskyFactor_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape(  [ 1._RK, 0._RK, -1._RK &
                                                                , 0._RK, 2._RK, -0._RK &
                                                                , 1._RK, 0._RK, -3._RK ], shape = shape(PosDefMat) )
        real(RK)                :: CholeskyLower(nd,nd), CholeskyDiagonal(nd)
        CholeskyLower = PosDefMat
        call getCholeskyFactor(nd = nd, PosDefMat = CholeskyLower, Diagonal = CholeskyDiagonal)
        assertion = CholeskyDiagonal(1) < 0._RK
        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "CholeskyLower      = ", CholeskyLower
            write(Test%outputUnit,"(*(g0,:,', '))") "CholeskyDiagonal   = ", CholeskyDiagonal
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getCholeskyFactor_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getInvPosDefMatSqrtDet_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape(  [ 1._RK, 0._RK, 1._RK &
                                                                , 0._RK, 2._RK, 0._RK &
                                                                , 1._RK, 0._RK, 3._RK ], shape = shape(PosDefMat) )
        real(RK)    , parameter :: MatInvMat_ref(nd,nd) = reshape(  [ 1.500000000000000_RK, 0.000000000000000_RK, -0.50000000000000_RK &
                                                                    , 0.000000000000000_RK, 0.500000000000000_RK, 0.000000000000000_RK &
                                                                    , -0.50000000000000_RK, 0.000000000000000_RK, 0.500000000000000_RK ] &
                                                                    , shape = shape(MatInvMat_ref) )
        real(RK)    , parameter :: CholeskyDiagonal_ref(nd) = [ 1.000000000000000_RK, 1.414213562373095_RK, 1.414213562373095_RK ]
        real(RK)    , parameter :: sqrtDetInvPosDefMat_ref = 0.5_RK
        real(RK)                :: MatInvMat(nd,nd), sqrtDetInvPosDefMat
        real(RK), allocatable   :: MatInvMat_diff(:,:), sqrtDetInvPosDefMat_diff

        MatInvMat = PosDefMat

        call getInvPosDefMatSqrtDet(nd = nd, MatInvMat = MatInvMat, sqrtDetInvPosDefMat = sqrtDetInvPosDefMat)

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(MatInvMat_diff)) deallocate(MatInvMat_diff); allocate(MatInvMat_diff, mold = MatInvMat)

        MatInvMat_diff = abs(MatInvMat - MatInvMat_ref)
        sqrtDetInvPosDefMat_diff = abs(sqrtDetInvPosDefMat - sqrtDetInvPosDefMat_ref)

        assertion = all(MatInvMat_diff < tolerance) .and. sqrtDetInvPosDefMat_diff < tolerance

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "MatInvMat_ref  = ", MatInvMat_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "MatInvMat      = ", MatInvMat
            write(Test%outputUnit,"(*(g0,:,', '))") "MatInvMat_diff = ", MatInvMat_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "sqrtDetInvPosDefMat_ref   = ", sqrtDetInvPosDefMat_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "sqrtDetInvPosDefMat_diff  = ", sqrtDetInvPosDefMat
            write(Test%outputUnit,"(*(g0,:,', '))") "sqrtDetInvPosDefMat       = ", sqrtDetInvPosDefMat_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getInvPosDefMatSqrtDet_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getInvPosDefMatSqrtDet_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape(  [ 1._RK, 0._RK, -1._RK &
                                                                , 0._RK, 2._RK, -0._RK &
                                                                , 1._RK, 0._RK, -3._RK ], shape = shape(PosDefMat) )
        real(RK)                :: MatInvMat(nd,nd), sqrtDetInvPosDefMat

        MatInvMat = PosDefMat

        call getInvPosDefMatSqrtDet(nd = nd, MatInvMat = MatInvMat, sqrtDetInvPosDefMat = sqrtDetInvPosDefMat)

        assertion = sqrtDetInvPosDefMat < 0._RK

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "MatInvMat              = ", MatInvMat
            write(Test%outputUnit,"(*(g0,:,', '))") "sqrtDetInvPosDefMat    = ", sqrtDetInvPosDefMat
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getInvPosDefMatSqrtDet_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getInvMatFromCholFac() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: CholeskyLower(nd,nd) = reshape(  [ 1.000000000000000_RK, 0.000000000000000_RK, 1.000000000000000_RK &
                                                                    , 0.000000000000000_RK, 2.000000000000000_RK, 0.000000000000000_RK &
                                                                    , 1.000000000000000_RK, 0.000000000000000_RK, 3.000000000000000_RK ] &
                                                                    , shape = shape(CholeskyLower) )
        real(RK)    , parameter :: CholeskyDiagonal(nd) = [ 1.000000000000000_RK, 1.414213562373095_RK, 1.414213562373095_RK ]
        real(RK)    , parameter :: CholeskyDiagonal_ref(nd) = [ 1.000000000000000_RK, 1.414213562373095_RK, 1.414213562373095_RK ]
        real(RK)    , parameter :: InvMatFromCholFac_ref(nd,nd) = reshape(  [ 1.500000000000000_RK, 0.000000000000000_RK, -0.50000000000000_RK &
                                                                            , 0.000000000000000_RK, 0.500000000000000_RK, 0.000000000000000_RK &
                                                                            , -0.50000000000000_RK, 0.000000000000000_RK, 0.500000000000000_RK ] &
                                                                            , shape = shape(InvMatFromCholFac_ref) )
        real(RK)    , parameter :: sqrtDetInvPosDefMat_ref = 0.5_RK
        real(RK)                :: InvMatFromCholFac(nd,nd)
        real(RK), allocatable   :: InvMatFromCholFac_diff(:,:)

        InvMatFromCholFac = getInvMatFromCholFac(nd = nd, CholeskyLower = CholeskyLower, Diagonal = CholeskyDiagonal)

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(InvMatFromCholFac_diff)) deallocate(InvMatFromCholFac_diff); allocate(InvMatFromCholFac_diff, mold = InvMatFromCholFac)

        InvMatFromCholFac_diff = abs(InvMatFromCholFac - InvMatFromCholFac_ref)

        assertion = all(InvMatFromCholFac_diff < tolerance)

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "InvMatFromCholFac_ref  = ", InvMatFromCholFac_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "InvMatFromCholFac      = ", InvMatFromCholFac
            write(Test%outputUnit,"(*(g0,:,', '))") "InvMatFromCholFac_diff = ", InvMatFromCholFac_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getInvMatFromCholFac

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getInvPosDefMat_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
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

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "MatInvMat_ref  = ", MatInvMat_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "MatInvMat      = ", MatInvMat
            write(Test%outputUnit,"(*(g0,:,', '))") "MatInvMat_diff = ", MatInvMat_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getInvPosDefMat_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getInvPosDefMat_2() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape(  [ 1._RK, 0._RK, -1._RK &
                                                                , 0._RK, 2._RK, -0._RK &
                                                                , 1._RK, 0._RK, -3._RK ], shape = shape(PosDefMat) )
        real(RK)                :: MatInvMat(nd,nd)

        MatInvMat = getInvPosDefMat(nd = nd, PosDefMat = PosDefMat)

        assertion = MatInvMat(1,1) < 0._RK

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "MatInvMat      = ", MatInvMat
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getInvPosDefMat_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getInvMatDet_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape(  [ 1._RK, 0._RK, 1._RK &
                                                                , 0._RK, 2._RK, 0._RK &
                                                                , 1._RK, 0._RK, 3._RK ], shape = shape(PosDefMat) )
        real(RK)    , parameter :: MatInvMat_ref(nd,nd) = reshape(  [ 1.500000000000000_RK, 0.000000000000000_RK, -0.50000000000000_RK &
                                                                    , 0.000000000000000_RK, 0.500000000000000_RK, 0.000000000000000_RK &
                                                                    , -0.50000000000000_RK, 0.000000000000000_RK, 0.500000000000000_RK ] &
                                                                    , shape = shape(MatInvMat_ref) )
        real(RK)    , parameter :: MatrixLU_ref(nd,nd)  = reshape(  [ 1.000000000000000_RK, 0.000000000000000_RK, 1.000000000000000_RK &
                                                                    , 0.000000000000000_RK, 2.000000000000000_RK, 0.000000000000000_RK &
                                                                    , 1.000000000000000_RK, 0.000000000000000_RK, 2.000000000000000_RK ] &
                                                                    , shape = shape(MatInvMat_ref) )
        real(RK)    , parameter :: detInvMat_ref = 0.25_RK
        real(RK)                :: MatInvMat(nd,nd), detInvMat, detInvMat_diff
        real(RK), allocatable   :: MatrixLU(:,:), MatInvMat_diff(:,:), MatrixLU_diff(:,:)

        MatrixLU = PosDefMat

        call getInvMatDet(nd = nd, MatrixLU = MatrixLU, InverseMatrix = MatInvMat, detInvMat = detInvMat)

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(MatrixLU_diff)) deallocate(MatrixLU_diff); allocate(MatrixLU_diff, mold = MatrixLU)
        MatrixLU_diff = abs(MatrixLU - MatrixLU_ref)

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(MatInvMat_diff)) deallocate(MatInvMat_diff); allocate(MatInvMat_diff, mold = MatInvMat)
        MatInvMat_diff = abs(MatInvMat - MatInvMat_ref)

        detInvMat_diff = abs(detInvMat - detInvMat_ref)

        assertion = all(MatInvMat_diff < tolerance) .and. detInvMat_diff < tolerance

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "MatrixLU_ref  = ", MatrixLU_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "MatrixLU      = ", MatrixLU
            write(Test%outputUnit,"(*(g0,:,', '))") "MatrixLU_diff = ", MatrixLU_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "MatInvMat_ref  = ", MatInvMat_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "MatInvMat      = ", MatInvMat
            write(Test%outputUnit,"(*(g0,:,', '))") "MatInvMat_diff = ", MatInvMat_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "detInvMat_ref  = ", detInvMat_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "detInvMat      = ", detInvMat
            write(Test%outputUnit,"(*(g0,:,', '))") "detInvMat_diff = ", detInvMat_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getInvMatDet_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getInvMat_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
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

        InverseMatrix = getInvMat(nd = nd, Matrix = PosDefMat)

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(InverseMatrix_diff)) deallocate(InverseMatrix_diff); allocate(InverseMatrix_diff, mold = InverseMatrix)

        InverseMatrix_diff = abs(InverseMatrix - InverseMatrix_ref)

        assertion = all(InverseMatrix_diff < tolerance)

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "InverseMatrix_ref  = ", InverseMatrix_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "InverseMatrix      = ", InverseMatrix
            write(Test%outputUnit,"(*(g0,:,', '))") "InverseMatrix_diff = ", InverseMatrix_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getInvMat_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_multiplyMatrix_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape(  [ 1._RK, 0._RK, 1._RK &
                                                                , 0._RK, 2._RK, 0._RK &
                                                                , 1._RK, 0._RK, 3._RK ], shape = shape(PosDefMat) )
        real(RK)                :: MatrixProduct(nd,nd)
        real(RK), allocatable   :: MatrixProduct_diff(:,:), MatrixProduct_ref(:,:)

        MatrixProduct_ref = matmul(PosDefMat,PosDefMat)

        call multiplyMatrix(A = PosDefMat, rowsA = nd, colsA = nd, B = PosDefMat, rowsB = nd, colsB = nd, C = MatrixProduct)

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(MatrixProduct_diff)) deallocate(MatrixProduct_diff); allocate(MatrixProduct_diff, mold = MatrixProduct)

        MatrixProduct_diff = abs(MatrixProduct - MatrixProduct_ref)

        assertion = all(MatrixProduct_diff < tolerance)

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "MatrixProduct_ref  = ", MatrixProduct_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "MatrixProduct      = ", MatrixProduct
            write(Test%outputUnit,"(*(g0,:,', '))") "MatrixProduct_diff = ", MatrixProduct_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_multiplyMatrix_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getDeterminant_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape(  [ -1._RK, -0._RK, -1._RK &
                                                                , -0._RK, -2._RK, -0._RK &
                                                                , -1._RK, -0._RK, -3._RK ], shape = shape(PosDefMat) )
        real(RK)    , parameter :: determinant_ref = -4._RK
        real(RK)                :: determinant, determinant_diff

        determinant = getDeterminant(nd = nd, Matrix = PosDefMat)

        determinant_diff = abs(determinant - determinant_ref)

        assertion = determinant_diff < tolerance

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "determinant_ref  = ", determinant_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "determinant      = ", determinant
            write(Test%outputUnit,"(*(g0,:,', '))") "determinant_diff = ", determinant_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getDeterminant_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getSqrtDetPosDefMat_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape(  [ 1._RK, 0._RK, 1._RK &
                                                                , 0._RK, 2._RK, 0._RK &
                                                                , 1._RK, 0._RK, 3._RK ], shape = shape(PosDefMat) )
        real(RK)    , parameter :: sqrtDetPosDefMat_ref = 2._RK
        real(RK)                :: sqrtDetPosDefMat, sqrtDetPosDefMat_diff

        sqrtDetPosDefMat = getSqrtDetPosDefMat(nd = nd, PosDefMat = PosDefMat)

        sqrtDetPosDefMat_diff = abs(sqrtDetPosDefMat - sqrtDetPosDefMat_ref)

        assertion = sqrtDetPosDefMat_diff < tolerance

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "sqrtDetPosDefMat_ref  = ", sqrtDetPosDefMat_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "sqrtDetPosDefMat      = ", sqrtDetPosDefMat
            write(Test%outputUnit,"(*(g0,:,', '))") "sqrtDetPosDefMat_diff = ", sqrtDetPosDefMat_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getSqrtDetPosDefMat_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getSqrtDetPosDefMat_2() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape(  [ -1._RK, -0._RK, -1._RK &
                                                                , -0._RK, -2._RK, -0._RK &
                                                                , -1._RK, -0._RK, -3._RK ], shape = shape(PosDefMat) )
        real(RK)                :: sqrtDetPosDefMat

        sqrtDetPosDefMat = getSqrtDetPosDefMat(nd = nd, PosDefMat = PosDefMat)

        assertion = sqrtDetPosDefMat < 0._RK

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "sqrtDetPosDefMat      = ", sqrtDetPosDefMat
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getSqrtDetPosDefMat_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogSqrtDetPosDefMat_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape(  [ 1._RK, 0._RK, 1._RK &
                                                                , 0._RK, 2._RK, 0._RK &
                                                                , 1._RK, 0._RK, 3._RK ], shape = shape(PosDefMat) )
        real(RK)    , parameter :: logSqrtDetPosDefMat_ref = log(2._RK)
        real(RK)                :: logSqrtDetPosDefMat, logSqrtDetPosDefMat_diff
        real(RK), allocatable   :: Matrix(:,:)
        logical                 :: failed

        Matrix = PosDefMat

        call getLogSqrtDetPosDefMat(nd = nd, PosDefMat = Matrix, logSqrtDetPosDefMat = logSqrtDetPosDefMat, failed = failed)

        assertion = .not. failed
        if (.not. assertion) return

        logSqrtDetPosDefMat_diff = abs(logSqrtDetPosDefMat - logSqrtDetPosDefMat_ref)

        assertion = logSqrtDetPosDefMat_diff < tolerance

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "logSqrtDetPosDefMat_ref  = ", logSqrtDetPosDefMat_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "logSqrtDetPosDefMat      = ", logSqrtDetPosDefMat
            write(Test%outputUnit,"(*(g0,:,', '))") "logSqrtDetPosDefMat_diff = ", logSqrtDetPosDefMat_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogSqrtDetPosDefMat_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogSqrtDetPosDefMat_2() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape(  [ 1._RK, 0._RK, 1._RK &
                                                                , 0._RK, 2._RK, 0._RK &
                                                                , 1._RK, 0._RK, 3._RK ], shape = shape(PosDefMat) )
        real(RK)    , parameter :: logSqrtDetPosDefMat_ref = log(2._RK)
        real(RK)                :: logSqrtDetPosDefMat
        real(RK), allocatable   :: Matrix(:,:)
        logical                 :: failed

        Matrix = -PosDefMat

        call getLogSqrtDetPosDefMat(nd = nd, PosDefMat = Matrix, logSqrtDetPosDefMat = logSqrtDetPosDefMat, failed = failed)

        assertion = failed
        if (.not. assertion) return

    end function test_getLogSqrtDetPosDefMat_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isPosDef_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape(  [ 1._RK, 0._RK, 1._RK &
                                                                , 0._RK, 2._RK, 0._RK &
                                                                , 1._RK, 0._RK, 3._RK ], shape = shape(PosDefMat) )

        assertion = isPosDef(nd = nd, Matrix = PosDefMat)
    end function test_isPosDef_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isPosDef_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: PosDefMat(nd,nd) = reshape(  [ -1._RK, -0._RK, -1._RK &
                                                                , -0._RK, -2._RK, -0._RK &
                                                                , -1._RK, -0._RK, -3._RK ], shape = shape(PosDefMat) )

        assertion = .not. isPosDef(nd = nd, Matrix = PosDefMat)
    end function test_isPosDef_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_symmetrizeUpperSquareMatrix_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: UpperSquareMatrix_ref(nd,nd) = reshape(  [ -1._RK, -0._RK, -1._RK &
                                                                            , -0._RK, -2._RK, -0._RK &
                                                                            , -1._RK, -0._RK, -3._RK ], shape = shape(UpperSquareMatrix_ref) )
        real(RK)                :: UpperSquareMatrix(nd,nd)
        integer(IK)             :: i, j

        UpperSquareMatrix = 0._RK
        do j = 1, nd
            do i = 1, j
                UpperSquareMatrix(i,j) = UpperSquareMatrix_ref(i,j)
            end do
        end do

        call symmetrizeUpperSquareMatrix(nd,UpperSquareMatrix)
        assertion = all( abs(UpperSquareMatrix-UpperSquareMatrix_ref) < 1.e-14_RK )

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "UpperSquareMatrix_ref  = ", UpperSquareMatrix_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "UpperSquareMatrix      = ", UpperSquareMatrix
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_symmetrizeUpperSquareMatrix_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getOuterProd_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: Vector2(*) = [4._RK, 5._RK]
        real(RK)    , parameter :: Vector1(*) = [1._RK, 2._RK, 3._RK]
        real(RK)    , parameter :: OuterProduct_ref(3,2) = reshape( [ 4._RK,  8._RK, 12._RK &
                                                                    , 5._RK, 10._RK, 15._RK ], shape = shape(OuterProduct_ref) )
        real(RK), allocatable   :: OuterProduct(:,:), OuterProduct_diff(:,:)

        OuterProduct = getOuterProd(Vector1 = Vector1, Vector2 = Vector2)

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(OuterProduct_diff)) deallocate(OuterProduct_diff); allocate(OuterProduct_diff, mold = OuterProduct)

        OuterProduct_diff = abs(OuterProduct - OuterProduct_ref)

        assertion = all(OuterProduct_diff < tolerance)

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "OuterProduct_ref  = ", OuterProduct_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "OuterProduct      = ", OuterProduct
            write(Test%outputUnit,"(*(g0,:,', '))") "OuterProduct_diff = ", OuterProduct_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getOuterProd_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_sortPosDefMat_1() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none
        logical(IK)             :: assertion
        integer(IK), parameter  :: ColIndx(*) = [4]
        integer(IK), parameter  :: ColIndxMap(*) = [5]
        integer(IK)             :: rank
        real(RK), allocatable   :: PosDefMat(:,:), OutPosDefMat(:,:), OutPosDefMat_ref(:,:)
        integer(IK)             :: i, j

        assertion = .true.

        rank = 5
        allocate(PosDefMat(rank,rank))
        do j = 1,rank
            do i = 1,j
                PosDefMat(i,j) = i*10 + j
            end do
        end do

        ! switch the variable 4 with 5, such that the output remains a positive-definite matrix.

        OutPosDefMat = sortPosDefMat(rank, PosDefMat, 1_IK, ColIndx, ColIndxMap)

        OutPosDefMat_ref = reshape(  [ 11._RK, 12._RK, 13._RK, 15._RK, 14._RK &
                                    ,  0._RK, 22._RK, 23._RK, 25._RK, 24._RK &
                                    ,  0._RK,  0._RK, 33._RK, 35._RK, 34._RK &
                                    ,  0._RK,  0._RK,  0._RK, 55._RK, 45._RK &
                                    ,  0._RK,  0._RK,  0._RK,  0._RK, 44._RK ] &
                                    , shape = [rank,rank] )
        OutPosDefMat_ref = transpose(OutPosDefMat_ref)
        do j = 1, rank
            do i = 1, j
                assertion = assertion .and. OutPosDefMat_ref(i,j) == OutPosDefMat(i,j)
            end do
        end do

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START

            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "OutPosDefMat_ref:"
            do i = 1,rank
                write(Test%outputUnit,"(*(F7.1))") (OutPosDefMat_ref(i,j),j=1,rank)
            end do
            write(Test%outputUnit,"(*(g0))")

            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "Output Positive-Definite Matrix with variables 4 and 5 swapped:"
            do i = 1,rank
                write(Test%outputUnit,"(*(F7.1))") (OutPosDefMat(i,j),j=1,rank)
            end do
            write(Test%outputUnit,"(*(g0))")

            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "Original Matrix:"
            do i = 1,rank
                write(Test%outputUnit,"(*(F7.1))") (PosDefMat(i,j),j=1,rank)
            end do

        end if
        ! LCOV_EXCL_STOP

    end function test_sortPosDefMat_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getRegresCoef_1() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none
        logical                 ::  assertion
        real(RK)                ::  normalizedDifference
        integer , parameter     ::  rankPDM = 4, rankS11 = 3, rankS22 = 1
        real(RK), parameter     ::  PosDefMat(rankPDM,rankPDM) = reshape( &
                                    [ 4.414182620515998_RK, 1.173760167060120_RK, 0.757607629189287_RK, 5.075277296976230_RK &
                                    , 1.173760167060120_RK, 0.866956750091570_RK, 0.310654936099342_RK, 1.621274787164182_RK &
                                    , 0.757607629189287_RK, 0.310654936099342_RK, 0.955157221699132_RK, 1.254186231887444_RK &
                                    , 5.075277296976230_RK, 1.621274787164182_RK, 1.254186231887444_RK, 6.407791808961157_RK ] &
                                    , shape=shape(PosDefMat) )
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
        integer                 ::  i,j

        assertion = .true.

        call getRegresCoef  ( rankPDM           = rankPDM           &
                            , rankS11           = rankS11           &
                            , rankS22           = rankS22           &
                            , PosDefMat         = PosDefMat         &
                            , RegresCoefMat     = RegresCoefMat     &
                            , SchurComplement   = SchurComplement   &
                            )

        do i = 1,rankS11
            !write(Test%outputUnit,"(*(F22.15))") (RegresCoefMat(i,j),j=1,rankS22), (RegresCoefMatRef(i,j),j=1,rankS22)
            do j = 1,rankS22
                normalizedDifference = abs(RegresCoefMatRef(i,j) - RegresCoefMat(i,j)) / (RegresCoefMatRef(i,j) + RegresCoefMat(i,j))
                assertion = assertion .and. normalizedDifference < 1.e-5_RK
            end do
        end do

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "Original Covariance Matrix:"
            do i = 1,rankPDM
                write(Test%outputUnit,"(*(F22.15))") (PosDefMat(i,j),j=1,rankPDM)
            end do
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "RegresCoefMat, RegresCoefMatRef, difference, normalizedDifference:"
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "SchurComplement, SchurComplementRef, difference, normalizedDifference:"
            write(Test%outputUnit,"(*(g0))")
            do i = 1,rankS11
                do j = 1,rankS22
                    write(Test%outputUnit,"(*(F22.15))" ) RegresCoefMat(i,j) &
                                                        , RegresCoefMatRef(i,j) &
                                                        , abs(RegresCoefMat(i,j)-RegresCoefMatRef(i,j)) &
                                                        , normalizedDifference
                end do
            end do
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

        do i = 1,rankS11
            do j = 1,rankS11
                normalizedDifference = abs(SchurComplementRef(i,j) - SchurComplement(i,j)) / (SchurComplementRef(i,j) + SchurComplement(i,j))
                assertion = assertion .and. normalizedDifference < 1.e-6_RK
            end do
        end do

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            do i = 1,rankS11
                do j = 1,rankS11
                    write(Test%outputUnit,"(*(F22.15))" ) SchurComplement(i,j) &
                                                        , SchurComplementRef(i,j) &
                                                        , abs(SchurComplement(i,j)-SchurComplementRef(i,j)) &
                                                        , normalizedDifference
                end do
            end do
        end if
        ! LCOV_EXCL_STOP

    end function test_getRegresCoef_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_Matrix_mod ! LCOV_EXCL_LINE