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
!>  This include file contains procedure implementation of the generic interfaces of [pm_matrixMulTri](@ref pm_matrixMulTri).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the blas dispatch function.
#if     trmv_ENABLED
#define blasTRXX blasTRMV
#elif   trsv_ENABLED
#define blasTRXX blasTRSV
#elif   trmm_ENABLED
#define blasTRXX blasTRMM
#elif   trsm_ENABLED
#define blasTRXX blasTRSM
#else
#error  "Unrecognized interface."
#endif
        ! Define the blas dispatch transposition.
#if     (CGMB_ENABLED && ONOA_ENABLED) || (CGMA_ENABLED && ONOB_ENABLED)
#define BLAS_TRANSA "N"
#elif   (CGMB_ENABLED && INVA_ENABLED) || (CGMA_ENABLED && INVB_ENABLED)
#define BLAS_TRANSA "N"
#elif   (CGMB_ENABLED && (OTSA_ENABLED || OTOA_ENABLED)) || (CGMA_ENABLED && (OTSB_ENABLED || OTOB_ENABLED))
#define BLAS_TRANSA "T"
#elif   (CGMB_ENABLED && (OTHA_ENABLED || OTUA_ENABLED)) || (CGMA_ENABLED && (OTHB_ENABLED || OTUB_ENABLED))
#define BLAS_TRANSA "C"
#else
#error  "Unrecognized interface."
#endif
        ! Define the blas dispatch matrix class.
#if     (trmv_ENABLED || trsv_ENABLED || trmm_ENABLED || trsm_ENABLED) && ((CGMB_ENABLED && (CLDA_ENABLED || CLUA_ENABLED)) || (CGMA_ENABLED && (CLDB_ENABLED || CLUB_ENABLED)))
#define BLAS_UPLO "L"
#elif   (trmv_ENABLED || trsv_ENABLED || trmm_ENABLED || trsm_ENABLED) && ((CGMB_ENABLED && (CUDA_ENABLED || CUUA_ENABLED)) || (CGMA_ENABLED && (CUDB_ENABLED || CUUB_ENABLED)))
#define BLAS_UPLO "U"
#else
#error  "Unrecognized interface."
#endif
        ! Define the constants.
#if     CK_ENABLED
#define TYPE_KIND complex(CKG)
        complex(CKG), parameter :: ZERO = (0._CKG, 0._CKG), ONE = (1._CKG, 0._CKG)
#elif   RK_ENABLED
#define TYPE_KIND real(RKG)
        real(RKG), parameter :: ZERO = 0._RKG, ONE = 1._RKG
#else
#error  "Unrecognized interface."
#endif
        ! Define the input and output matrices.
#if     CGMA_ENABLED
#define BLAS_SIDE "R"
#define SOLMAT matA
#define TRIMAT matB
#elif   CGMB_ENABLED
#define BLAS_SIDE "L"
#define SOLMAT matB
#define TRIMAT matA
#else
#error  "Unrecognized interface."
#endif
        ! Define the conjugation of the triangular matrix.
#if     CK_ENABLED && ((OTUA_ENABLED || OTUB_ENABLED) || (OTHA_ENABLED || OTHB_ENABLED))
#define GET_CONJG(X)conjg(X)
#elif   (ONOA_ENABLED || ONOB_ENABLED) || (INVA_ENABLED || INVB_ENABLED) || (OTOA_ENABLED || OTOB_ENABLED) || (OTUA_ENABLED || OTUB_ENABLED)
#define GET_CONJG(X)X
#else
#error  "Unrecognized interface."
#endif
        ! Define the diag/unit triangular matrices.
#if     (CGMB_ENABLED && (CLDA_ENABLED || CUDA_ENABLED)) || (CGMA_ENABLED && (CLDB_ENABLED || CUDB_ENABLED))
#define BLAS_DIAG "N"
#define CXD_ENABLED 1
#elif   (CGMB_ENABLED && (CLUA_ENABLED || CUUA_ENABLED)) || (CGMA_ENABLED && (CLUB_ENABLED || CUUB_ENABLED))
#define BLAS_DIAG "U"
#define CXD_ENABLED 0
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     setMatMulTri_ENABLED && (trmv_ENABLED || trsv_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     !(BLAS_ENABLED && DISPATCH_ENABLED)
        integer(IK) :: irow, icol, bstart, ib, jb
        TYPE_KIND :: temp
#endif
        ! Define the assumed shape of the matrices.
#if     ASS_ENABLED
        integer(IK) :: nrow
        integer(IK) , parameter :: incB = 1
        nrow = size(TRIMAT, 1, IK)
        CHECK_ASSERTION(__LINE__, size(matA, 1, IK) == size(matA, 2, IK), SK_"@setMatMulTri(): The condition `size(matA, 1) == size(matA, 2)` must hold. shape(matA) = "//getStr(shape(matA, IK))) ! fpp
        CHECK_ASSERTION(__LINE__, size(matA, 1, IK) == size(matB, 1, IK), SK_"@setMatMulTri(): The condition `size(matA, 1) == size(matB, 1)` must hold. shape(matA), shape(matB) = "//getStr([shape(matA), shape(matB)])) ! fpp
#elif   EXP_ENABLED
        ! Ensure offsets and subset bounds are logical.
        CHECK_ASSERTION(__LINE__, 0_IK <= nrow, SK_"@setMatMulTri(): The condition `0 <= nrow` must hold. nrow = "//getStr(nrow)) ! fpp
        CHECK_ASSERTION(__LINE__, 0_IK /= incB, SK_"@setMatMulTri(): The condition `0 /= incB` must hold. nrow = "//getStr(incB)) ! fpp
        CHECK_ASSERTION(__LINE__, all(0_IK <= 1_IK - lbound(matA, kind = IK)), SK_"@setMatMulTri(): The condition `all(0 <= [roffA, coffA])` must hold. roffA, coffA = "//getStr(1_IK - lbound(matA, kind = IK))) ! fpp
        CHECK_ASSERTION(__LINE__, abs(incB) * (nrow - 1) < size(matB, 1, IK), SK_"@setMatMulTri(): The condition `size(matB(1::abs(incB))) <= nrow` must hold. incB, nrow, size(matB) = "//getStr([incB, nrow, size(matB, kind = IK)])) ! fpp
        CHECK_ASSERTION(__LINE__, all(nrow <= ubound(matA, kind = IK)), SK_"@setMatMulTri(): The condition `all(nrow == shape(matA) - [roffA, coffA])` must hold. nrow, roffA, coffA, shape(matA) = "//getStr([nrow, 1_IK - lbound(matA, kind = IK), shape(matA, IK)])) ! fpp
        ! Set up the start point in matB if the increment is not unity.
        ! This will be (nrow - 1) * incB too small for descending loops.
#else
#error  "Unrecognized interface."
#endif
        ! Quick return if possible.
        if (nrow == 0_IK) return
#if     BLAS_ENABLED && DISPATCH_ENABLED
        call blasTRXX(BLAS_UPLO, BLAS_TRANSA, BLAS_DIAG, nrow, TRIMAT(1,1), size(TRIMAT, 1, IK), SOLMAT(1), incB)
#else
        ! Start the operations.
        ! In this version the elements of matA are accessed sequentially with one pass through matA.
        if (incB < 0_IK) then
            bstart = 1 - (nrow - 1) * incB
        elseif (incB /= 1) then
            bstart = 1
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     CGMB_ENABLED && ONOB_ENABLED && ONOA_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! BLAS 2 TRMV - Form  matB := matA * matB.
#if     CUDA_ENABLED || CUUA_ENABLED
        if (incB == 1) then
            do icol = 1, nrow
                if (matB(icol) /= ZERO) then
                    temp = matB(icol)
                    do irow = 1, icol - 1
                        matB(irow) = matB(irow) + temp * matA(irow, icol)
                    end do
#if                 CXD_ENABLED
                    matB(icol) = matB(icol) * matA(icol, icol)
#endif
                end if
            end do
        else
            jb = bstart
            do icol = 1, nrow
                if (matB(jb) /= ZERO) then
                    temp = matB(jb)
                    ib = bstart
                    do irow = 1, icol - 1
                        matB(ib) = matB(ib) + temp * matA(irow, icol)
                        ib = ib + incB
                    end do
#if                 CXD_ENABLED
                    matB(jb) = matB(jb) * matA(icol, icol)
#endif
                end if
                jb = jb + incB
            end do
        end if
#elif   CLDA_ENABLED || CLUA_ENABLED
        if (incB == 1) then
            do icol = nrow, 1, -1
                if (matB(icol) /= ZERO) then
                    temp = matB(icol)
                    do irow = nrow, icol + 1, -1
                        matB(irow) = matB(irow) + temp * matA(irow, icol)
                    end do
#if                 CXD_ENABLED
                    matB(icol) = matB(icol) * matA(icol, icol)
#endif
                end if
            end do
        else
            bstart = bstart + (nrow - 1) * incB
            jb = bstart
            do icol = nrow, 1, -1
                if (matB(jb) /= ZERO) then
                    temp = matB(jb)
                    ib = bstart
                    do irow = nrow, icol + 1, -1
                        matB(ib) = matB(ib) + temp * matA(irow, icol)
                        ib = ib - incB
                    end do
#if                 CXD_ENABLED
                    matB(jb) = matB(jb) * matA(icol, icol)
#endif
                end if
                jb = jb - incB
            end do
        end if
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   CGMB_ENABLED && ONOB_ENABLED && INVA_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! BLAS TRSV - Form matB := inv(matA) * matB.

#if     CUDA_ENABLED || CUUA_ENABLED
        if (incB == 1) then
            do icol = nrow, 1, -1
                if (matB(icol) /= ZERO) then
#if                 CXD_ENABLED
                    matB(icol) = matB(icol) / matA(icol, icol)
#endif
                    temp = matB(icol)
                    do irow = icol - 1, 1, -1
                        matB(irow) = matB(irow) - temp * matA(irow,icol)
                    end do
                end if
            end do
        else
            jb = bstart + (nrow - 1) * incB
            do icol = nrow, 1, -1
                if (matB(jb) /= ZERO) then
#if                 CXD_ENABLED
                    matB(jb) = matB(jb) / matA(icol, icol)
#endif
                    temp = matB(jb)
                    ib = jb
                    do irow = icol - 1, 1, -1
                        ib = ib - incB
                        matB(ib) = matB(ib) - temp * matA(irow, icol)
                    end do
                end if
                jb = jb - incB
            end do
        end if
#elif   CLDA_ENABLED || CLUA_ENABLED
        if (incB == 1) then
            do icol = 1, nrow
                if (matB(icol) /= ZERO) then
#if                 CXD_ENABLED
                    matB(icol) = matB(icol) / matA(icol, icol)
#endif
                    temp = matB(icol)
                    do irow = icol + 1, nrow
                        matB(irow) = matB(irow) - temp * matA(irow, icol)
                    end do
                end if
            end do
        else
            jb = bstart
            do icol = 1, nrow
                if (matB(jb) /= ZERO) then
#if                 CXD_ENABLED
                    matB(jb) = matB(jb) / matA(icol, icol)
#endif
                    temp = matB(jb)
                    ib = jb
                    do irow = icol + 1, nrow
                        ib = ib + incB
                        matB(ib) = matB(ib) - temp * matA(irow, icol)
                    end do
                end if
                jb = jb + incB
            end do
        end if
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   CGMB_ENABLED && ONOB_ENABLED && (OTSA_ENABLED || OTHA_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! BLAS TRMV - Form matB := matA^T * matB or matB := matA^H * matB.

#if     CUDA_ENABLED || CUUA_ENABLED
        if (incB == 1) then
            do icol = nrow, 1, -1
                temp = matB(icol)
#if             CXD_ENABLED
                temp = temp * GET_CONJG(matA(icol, icol))
#endif
                do irow = icol - 1, 1, -1
                    temp = temp + GET_CONJG(matA(irow, icol)) * matB(irow)
                end do
                matB(icol) = temp
            end do
        else
            jb = bstart + (nrow - 1) * incB
            do icol = nrow, 1, -1
                temp = matB(jb)
                ib = jb
#if             CXD_ENABLED
                temp = temp*GET_CONJG(matA(icol,icol))
#endif
                do irow = icol - 1, 1, -1
                    ib = ib - incB
                    temp = temp + GET_CONJG(matA(irow, icol)) * matB(ib)
                end do
                matB(jb) = temp
                jb = jb - incB
            end do
        end if
#elif   CLDA_ENABLED || CLUA_ENABLED
        if (incB == 1) then
            do icol = 1, nrow
                temp = matB(icol)
#if             CXD_ENABLED
                temp = temp * GET_CONJG(matA(icol, icol))
#endif
                do irow = icol + 1, nrow
                    temp = temp + GET_CONJG(matA(irow, icol)) * matB(irow)
                end do
                matB(icol) = temp
            end do
        else
            jb = bstart
            do icol = 1, nrow
                temp = matB(jb)
                ib = jb
#if             CXD_ENABLED
                temp = temp * GET_CONJG(matA(icol, icol))
#endif
                do irow = icol + 1, nrow
                    ib = ib + incB
                    temp = temp + GET_CONJG(matA(irow,icol))*matB(ib)
                end do
                matB(jb) = temp
                jb = jb + incB
            end do
        end if
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   CGMB_ENABLED && ONOB_ENABLED && (OTOA_ENABLED || OTUA_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! BLAS 2 TRSV - Form matB := inv(matA^T) * matB  or  matB := inv(matA^H) * matB.
#if     CUDA_ENABLED || CUUA_ENABLED
        if (incB == 1) then
            do icol = 1, nrow
                temp = matB(icol)
                do irow = 1, icol - 1
                    temp = temp - GET_CONJG(matA(irow, icol)) * matB(irow)
                end do
#if             CXD_ENABLED
                temp = temp / GET_CONJG(matA(icol, icol))
#endif
                matB(icol) = temp
            end do
        else
            jb = bstart
            do icol = 1, nrow
                ib = bstart
                temp = matB(jb)
                do irow = 1, icol - 1
                    temp = temp - GET_CONJG(matA(irow, icol)) * matB(ib)
                    ib = ib + incB
                end do
#if             CXD_ENABLED
                temp = temp / GET_CONJG(matA(icol, icol))
#endif
                matB(jb) = temp
                jb = jb + incB
            end do
        end if
#elif   CLDA_ENABLED || CLUA_ENABLED
        if (incB == 1) then
            do icol = nrow, 1, -1
                temp = matB(icol)
                do irow = nrow,icol + 1, -1
                    temp = temp - GET_CONJG(matA(irow, icol)) * matB(irow)
                end do
#if             CXD_ENABLED
                temp = temp / GET_CONJG(matA(icol, icol))
#endif
                matB(icol) = temp
            end do
        else
            bstart = bstart + (nrow-1)*incB
            jb = bstart
            do icol = nrow, 1, -1
                ib = bstart
                temp = matB(jb)
                do irow = nrow, icol + 1, -1
                    temp = temp - GET_CONJG(matA(irow, icol)) * matB(ib)
                    ib = ib - incB
                end do
#if             CXD_ENABLED
                temp = temp / GET_CONJG(matA(icol, icol))
#endif
                matB(jb) = temp
                jb = jb - incB
            end do
        end if
#else
#error  "Unrecognized interface."
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatMulTri_ENABLED && (trmm_ENABLED || trsm_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     !(BLAS_ENABLED && DISPATCH_ENABLED)
        integer(IK) :: irow, icol, idum
#endif
        ! Define the assumed shape of the matrices.
#if     ASS_ENABLED
#define ALPHA alpha_def
        integer(IK) :: nrow, ncol
        TYPE_KIND :: alpha_def
        if (present(alpha)) then
            alpha_def = alpha
        else
            alpha_def = ONE
        end if
        nrow = size(SOLMAT, 1, IK)
        ncol = size(SOLMAT, 2, IK)
#elif   EXP_ENABLED
        ! Ensure offsets and subset bounds are logical.
        CHECK_ASSERTION(__LINE__, all(0_IK <= [nrow, ncol]), SK_"@setMatMulTri(): The condition `all(0_IK <= [nrow, ncol])` must hold. nrow, ncol = "//getStr([nrow, ncol])) ! fpp
        CHECK_ASSERTION(__LINE__, all(0_IK <= 1_IK - lbound(matA, kind = IK)), SK_"@setMatMulTri(): The condition `all(0 <= [roffA, coffA])` must hold. roffA, coffA = "//getStr(1_IK - lbound(matA, kind = IK))) ! fpp
        CHECK_ASSERTION(__LINE__, all(0_IK <= 1_IK - lbound(matB, kind = IK)), SK_"@setMatMulTri(): The condition `all(0 <= [roffB, coffB])` must hold. roffB, coffB = "//getStr(1_IK - lbound(matB, kind = IK))) ! fpp
#else
#error  "Unrecognized interface."
#endif
        ! Ensure subset of `TRIMAT` matches the corresponding subset of `SOLMAT`.
#if     CGMA_ENABLED
        CHECK_ASSERTION(__LINE__, all([nrow, ncol] <= ubound(matA, kind = IK)), SK_"@setMatMulTri(): The condition `all([nrow + roffA, ncol + coffA] <= shape(matA, kind = IK))` must hold. nrow, ncol, roffA, coffA, shape(matA) = "//getStr([nrow, ncol, 1_IK - lbound(matA, kind = IK), shape(matA, IK)])) ! fpp
        CHECK_ASSERTION(__LINE__, all(ncol <= ubound(matB, kind = IK)), SK_"@setMatMulTri(): The condition `all([roffB, coffB] + ncol <= shape(matB))` must hold. roffB, coffB, ncol, shape(matB) = "//getStr([1_IK - lbound(matB, kind = IK), ncol, shape(matB, IK)])) ! fpp
#elif   CGMB_ENABLED
        CHECK_ASSERTION(__LINE__, all(nrow <= ubound(matA, kind = IK)), \
        SK_"@setMatMulTri(): The condition `all([roffA, coffA] + nrow <= shape(matA))` must hold. roffA, coffA, nrow, shape(matA) = "//getStr([1_IK - lbound(matA, kind = IK), nrow, shape(matA, IK)])) ! fpp
        CHECK_ASSERTION(__LINE__, all([nrow, ncol] <= ubound(matB, kind = IK)), SK_"@setMatMulTri(): The condition `all([nrow + roffB, ncol + coffB] <= shape(matB, kind = IK))` must hold. nrow, ncol, roffB, coffB, shape(matB) = "//getStr([nrow, ncol, 1_IK - lbound(matB, kind = IK), shape(matB, IK), ubound(matB, kind = IK)])) ! fpp
#endif
        ! Set the iteration limits based on the storage format.
#if     (CUDA_ENABLED || CUUA_ENABLED) || (CUDB_ENABLED || CUUB_ENABLED)
#define GET_HALF_RANGE(i,j,k) i, j - 1
#elif   (CLDA_ENABLED || CLUA_ENABLED) || (CLDB_ENABLED || CLUB_ENABLED)
#define GET_HALF_RANGE(i,j,k) j + 1, k
#else
#error  "Unrecognized interface."
#endif
        ! Quick return if possible.
        if (nrow == 0_IK .or. ncol == 0_IK) return
        if (ALPHA == ZERO) then
            SOLMAT(1 : nrow, 1 : ncol) = ZERO
            return
        end if
#if     BLAS_ENABLED && DISPATCH_ENABLED
        call blasTRXX(BLAS_SIDE, BLAS_UPLO, BLAS_TRANSA, BLAS_DIAG, nrow, ncol, ALPHA, TRIMAT(1,1), size(TRIMAT, 1, IK), SOLMAT(1,1), size(SOLMAT, 1, IK))
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     trmm_ENABLED && CGMB_ENABLED && ONOB_ENABLED && ONOA_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! BLAS 3 TRMM - Form SOLMAT := alpha * TRIMAT * SOLMAT.
        block
            TYPE_KIND :: temp
#if         CUDA_ENABLED || CUUA_ENABLED
            do icol = 1, ncol
                do idum = 1, nrow
                    if (SOLMAT(idum, icol) /= ZERO) then
                        temp = ALPHA * SOLMAT(idum, icol)
                        do irow = 1, idum - 1
                            SOLMAT(irow, icol) = SOLMAT(irow, icol) + temp * TRIMAT(irow, idum)
                        end do
#if                     CXD_ENABLED
                        temp = temp * TRIMAT(idum, idum)
#endif
                        SOLMAT(idum, icol) = temp
                    end if
                end do
            end do
#elif       CLDA_ENABLED || CLUA_ENABLED
            do icol = 1, ncol
                do idum = nrow, 1, -1
                    if (SOLMAT(idum, icol) /= ZERO) then
                        temp = ALPHA * SOLMAT(idum, icol)
                        SOLMAT(idum, icol) = temp
#if                     CXD_ENABLED
                        SOLMAT(idum, icol) = SOLMAT(idum, icol) * TRIMAT(idum, idum)
#endif
                        do irow = idum + 1, nrow
                            SOLMAT(irow, icol) = SOLMAT(irow, icol) + temp * TRIMAT(irow, idum)
                        end do
                    end if
                end do
            end do
#else
#error      "Unrecognized interface."
#endif
        end block

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   trsm_ENABLED && CGMB_ENABLED && ONOB_ENABLED && INVA_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! BLAS 3 TRSM - Form SOLMAT := alpha * inv(TRIMAT) * SOLMAT.
#if     CUDA_ENABLED || CUUA_ENABLED
#define GET_FULL_RANGE(i,j) j, i, -1
#elif   CLDA_ENABLED || CLUA_ENABLED
#define GET_FULL_RANGE(i,j) i, j
#else
#error  "Unrecognized interface."
#endif
        do icol = 1, ncol
            if (ALPHA /= ONE) then
                do  irow = 1, nrow
                    SOLMAT(irow, icol) = ALPHA * SOLMAT(irow, icol)
                end do
            end if
            do idum = GET_FULL_RANGE(1, nrow)
                if (SOLMAT(idum, icol) /= ZERO) then
#if                 CXD_ENABLED
                    SOLMAT(idum, icol) = SOLMAT(idum, icol) / TRIMAT(idum, idum)
#endif
                    do irow = GET_HALF_RANGE(1, idum, nrow)
                        SOLMAT(irow, icol) = SOLMAT(irow, icol) - SOLMAT(idum, icol) * TRIMAT(irow, idum)
                    end do
                end if
            end do
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   trmm_ENABLED && CGMB_ENABLED && ONOB_ENABLED && (OTSA_ENABLED || OTHA_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! BLAS 3 TRMM - Form SOLMAT := alpha * TRIMAT^T * SOLMAT or SOLMAT := alpha * TRIMAT^H * SOLMAT

        block
            TYPE_KIND :: temp
            do icol = 1, ncol
#if             CUDA_ENABLED || CUUA_ENABLED
                do irow = nrow, 1, -1
                    temp = SOLMAT(irow, icol)
#if                 CXD_ENABLED
                    temp = temp * GET_CONJG(TRIMAT(irow, irow))
#endif
                    do idum = 1, irow - 1
                        temp = temp + GET_CONJG(TRIMAT(idum, irow)) * SOLMAT(idum, icol)
                    end do
                    SOLMAT(irow, icol) = ALPHA * temp
                end do
#elif           CLDA_ENABLED || CLUA_ENABLED
                do irow = 1, nrow
                    temp = SOLMAT(irow, icol)
#if                 CXD_ENABLED
                    temp = temp * GET_CONJG(TRIMAT(irow, irow))
#endif
                    do idum = irow + 1, nrow
                        temp = temp + GET_CONJG(TRIMAT(idum,irow))*SOLMAT(idum, icol)
                    end do
                    SOLMAT(irow, icol) = ALPHA*temp
                end do
#else
#error          "Unrecognized interface."
#endif
            end do
        end block

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   trsm_ENABLED && CGMB_ENABLED && ONOB_ENABLED && (OTOA_ENABLED || OTUA_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! BLAS 3 TRSM - Form SOLMAT := alpha * inv(transpose(TRIMAT)) * SOLMAT
        ! BLAS 3 TRSM - Form SOLMAT := alpha * inv(conjg(transpose(TRIMAT))) * SOLMAT
#if     CUDA_ENABLED || CUUA_ENABLED
#define GET_FULL_RANGE(i,j) i, j
#elif   CLDA_ENABLED || CLUA_ENABLED
#define GET_FULL_RANGE(i,j) j, i, -1
#else
#error  "Unrecognized interface."
#endif
        !error stop "SOLMAT := alpha * inv(transpose(TRIMAT)) * SOLMAT"
        block
            TYPE_KIND :: temp
            do icol = 1, ncol
                do irow = GET_FULL_RANGE(1, nrow)
                    temp = ALPHA * SOLMAT(irow, icol)
                    do idum = GET_HALF_RANGE(1, irow, nrow)
                        temp = temp - GET_CONJG(TRIMAT(idum, irow)) * SOLMAT(idum, icol)
                    end do
#if                 CXD_ENABLED
                    temp = temp / GET_CONJG(TRIMAT(irow, irow))
#endif
                    SOLMAT(irow, icol) = temp
                end do
            end do
        end block

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   trmm_ENABLED && CGMA_ENABLED && ONOA_ENABLED && ONOB_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! BLAS 3 TRMM - Form  SOLMAT := alpha * SOLMAT * TRIMAT.

        block
            TYPE_KIND :: temp
#if         CUDB_ENABLED || CUUB_ENABLED
            do icol = ncol, 1, -1
                temp = ALPHA
#if             CXD_ENABLED
                temp = temp * TRIMAT(icol, icol)
#endif
                do irow = 1, nrow
                    SOLMAT(irow, icol) = temp * SOLMAT(irow, icol)
                end do
                do idum = 1, icol - 1
                    if (TRIMAT(idum, icol) /= ZERO) then
                        temp = ALPHA * TRIMAT(idum, icol)
                        do irow = 1, nrow
                            SOLMAT(irow, icol) = SOLMAT(irow, icol) + temp * SOLMAT(irow, idum)
                        end do
                    end if
                end do
            end do
#elif       CLDB_ENABLED || CLUB_ENABLED
            do  icol = 1, ncol
                temp = ALPHA
#if             CXD_ENABLED
                temp = temp * TRIMAT(icol, icol)
#endif
                do irow = 1, nrow
                    SOLMAT(irow, icol) = temp * SOLMAT(irow, icol)
                end do
                do idum = icol + 1, ncol
                    if (TRIMAT(idum, icol) /= ZERO) then
                        temp = ALPHA * TRIMAT(idum, icol)
                        do irow = 1, nrow
                            SOLMAT(irow, icol) = SOLMAT(irow, icol) + temp * SOLMAT(irow, idum)
                        end do
                    end if
                end do
            end do
#else
#error      "Unrecognized interface."
#endif
        end block

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   trsm_ENABLED && CGMA_ENABLED && ONOA_ENABLED && INVB_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! BLAS 3 TRSM - Form SOLMAT := alpha * SOLMAT * inv(TRIMAT).
#if     CUDB_ENABLED || CUUB_ENABLED
#define GET_FULL_RANGE(i,j) i, j
#elif   CLDB_ENABLED || CLUB_ENABLED
#define GET_FULL_RANGE(i,j) j, i, -1
#else
#error  "Unrecognized interface."
#endif
        block
#if         CXD_ENABLED
            TYPE_KIND :: temp
#endif
            do icol = GET_FULL_RANGE(1, ncol)
                if (ALPHA /= ONE) then
                    do irow = 1, nrow
                        SOLMAT(irow,icol) = ALPHA * SOLMAT(irow, icol)
                    end do
                end if
                do idum = GET_HALF_RANGE(1, icol, ncol)
                    if (TRIMAT(idum, icol) /= ZERO) then
                        do irow = 1, nrow
                            SOLMAT(irow, icol) = SOLMAT(irow, icol) - TRIMAT(idum, icol) * SOLMAT(irow, idum)
                        end do
                    end if
                end do
#if             CXD_ENABLED
                temp = ONE / TRIMAT(icol, icol)
                do irow = 1, nrow
                    SOLMAT(irow, icol) = temp * SOLMAT(irow, icol)
                end do
#endif
            end do
        end block

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   trmm_ENABLED && CGMA_ENABLED && ONOA_ENABLED && (OTSB_ENABLED || OTHB_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! BLAS 3 TRMM - Form SOLMAT := alpha * SOLMAT * TRIMAT**T or SOLMAT := alpha * SOLMAT * TRIMAT^H.

        block
            TYPE_KIND :: temp
#if         CUDB_ENABLED || CUUB_ENABLED
            do idum = 1, ncol
                do icol = 1, idum - 1
                    if (TRIMAT(icol, idum) /= ZERO) then
                        temp = ALPHA * GET_CONJG(TRIMAT(icol, idum))
                        do irow = 1, nrow
                            SOLMAT(irow, icol) = SOLMAT(irow, icol) + temp * SOLMAT(irow, idum)
                        end do
                    end if
                end do
                temp = ALPHA
#if             CXD_ENABLED
                temp = temp * GET_CONJG(TRIMAT(idum, idum))
#endif
                if (temp /= ONE) then
                    do irow = 1, nrow
                        SOLMAT(irow, idum) = temp * SOLMAT(irow, idum)
                    end do
                end if
            end do
#elif       CLDB_ENABLED || CLUB_ENABLED
            do idum = ncol, 1, -1
                do icol = idum + 1, ncol
                    if (TRIMAT(icol, idum) /= ZERO) then
                        temp = ALPHA * GET_CONJG(TRIMAT(icol, idum))
                        do irow = 1, nrow
                            SOLMAT(irow, icol) = SOLMAT(irow, icol) + temp * SOLMAT(irow, idum)
                        end do
                    end if
                end do
                temp = ALPHA
#if             CXD_ENABLED
                temp = temp * GET_CONJG(TRIMAT(idum, idum))
#endif
                if (temp /= ONE) then
                    do irow = 1, nrow
                        SOLMAT(irow, idum) = temp * SOLMAT(irow, idum)
                    end do
                end if
            end do
#else
#error      "Unrecognized interface."
#endif
        end block

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   trsm_ENABLED && CGMA_ENABLED && ONOA_ENABLED && (OTOB_ENABLED || OTUB_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! BLAS 3 TRSM - Form SOLMAT := alpha * SOLMAT * inv(TRIMAT**T).
        ! BLAS 3 TRSM - Form SOLMAT := alpha * SOLMAT * inv(TRIMAT**H).
#if     CUDB_ENABLED || CUUB_ENABLED
#define GET_FULL_RANGE(i,j) j, i, -1
#elif   CLDB_ENABLED || CLUB_ENABLED
#define GET_FULL_RANGE(i,j) i, j
#else
#error  "Unrecognized interface."
#endif
        block
            TYPE_KIND :: temp
            do idum = GET_FULL_RANGE(1, ncol)
#if             CXD_ENABLED
                temp = ONE / GET_CONJG(TRIMAT(idum, idum))
                do irow = 1, nrow
                    SOLMAT(irow, idum) = temp * SOLMAT(irow, idum)
                end do
#endif
                do icol = GET_HALF_RANGE(1, idum, ncol)
                    if (TRIMAT(icol, idum) /= ZERO) then
                        temp = GET_CONJG(TRIMAT(icol, idum))
                        do irow = 1, nrow
                            SOLMAT(irow, icol) = SOLMAT(irow, icol) - temp * SOLMAT(irow, idum)
                        end do
                    end if
                end do
                if (ALPHA /= ONE) then
                    do irow = 1, nrow
                        SOLMAT(irow, idum) = ALPHA * SOLMAT(irow, idum)
                    end do
                end if
            end do
        end block
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  GET_FULL_RANGE
#undef  GET_HALF_RANGE
#undef  CXD_ENABLED
#undef  BLAS_TRANSA
#undef  BLAS_SIDE
#undef  BLAS_UPLO
#undef  BLAS_DIAG
#undef  TYPE_KIND
#undef  GET_CONJG
#undef  blasTRXX
#undef  SOLMAT
#undef  TRIMAT
#undef  ALPHA