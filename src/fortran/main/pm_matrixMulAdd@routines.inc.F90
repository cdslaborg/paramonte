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
!>  This include file contains procedure implementation of the generic interfaces of [pm_matrixMulAdd](@ref pm_matrixMulAdd).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the default optional constants.
#if     ASS_ENABLED && (gemv_ENABLED || spmv_ENABLED || hpmv_ENABLED || symm_ENABLED || symv_ENABLED || hemv_ENABLED || hemm_ENABLED || gemm_ENABLED)
#define BETA beta_def
#define ALPHA alpha_def
#define DECLARE_SET_DEFAULT_ALPHA_BETA \
TYPE_KIND :: alpha_def, beta_def; \
beta_def = ONE; alpha_def = ONE; \
if (present(beta)) beta_def = beta; \
if (present(alpha)) alpha_def = alpha;
#elif   !EXP_ENABLED
#error  "Unrecognized interface."
#endif
        ! Define the transposition of `matA` for BLAS gemv/gemm call routines.
#if     TNA_ENABLED
#define TMA "N"
#elif   TSA_ENABLED
#define TMA "T"
#elif   THA_ENABLED
#define TMA "C"
#endif
        ! Define the transposition of `matB` for BLAS gemm call routines.
#if     TNB_ENABLED
#define TMB "N"
#elif   TSB_ENABLED
#define TMB "T"
#elif   THB_ENABLED
#define TMB "C"
#endif
        ! Define the call interface for BLAS symm/hemm call routines.
#if     CHA_ENABLED && CHB_ENABLED
#define SHMM blasHEMM
#elif   CSA_ENABLED || CNA_ENABLED || CSB_ENABLED || CNB_ENABLED
#define SHMM blasSYMM
#endif
        ! Define blas gemm call interface.
#define GEMM(transA, tranB, alpha, beta) \
blasGEMM(transA, tranB, nrow, ncol, ndum, alpha, matA(1, 1), size(matA, 1, IK), matB(1, 1), size(matB, 1, IK), beta, matC(1, 1), size(matC, 1, IK))
        ! Define matA transpose form.
#if     CK_ENABLED && (CHA_ENABLED || THA_ENABLED)
#define CONJG_A(X) conjg(X)
#elif   THA_ENABLED || TSA_ENABLED || TNA_ENABLED || CHA_ENABLED || CSA_ENABLED || CNA_ENABLED
#define CONJG_A(X) X
#else
#error  "Unrecognized interface."
#endif
        ! Define matB transpose form.
#if     CK_ENABLED && (CHB_ENABLED || THB_ENABLED)
#define CONJG_B(X) conjg(X)
#elif   THB_ENABLED || TSB_ENABLED || TNB_ENABLED || CHB_ENABLED || CSB_ENABLED || CNB_ENABLED
#define CONJG_B(X) X
#else
#error  "Unrecognized interface."
#endif
        ! Define the `real` component function.
#if     CK_ENABLED && (THA_ENABLED || THB_ENABLED || CHA_ENABLED || CHB_ENABLED)
#define GET_RE(X) X%re
#elif   THA_ENABLED || TSA_ENABLED || TNA_ENABLED || THB_ENABLED || TSB_ENABLED || TNB_ENABLED || \
        CHA_ENABLED || CSA_ENABLED || CNA_ENABLED || CHB_ENABLED || CSB_ENABLED || CNB_ENABLED
#define GET_RE(X) X
#else
#error  "Unrecognized interface."
#endif
        ! Define the additional temporary variable for triangular `matA`.
!#if     SLA_ENABLED || SUA_ENABLED || D2_D1_ENABLED || D1_D1_ENABLED
!#define ndum nrow
!#elif   SLB_ENABLED || SUB_ENABLED
!#define ndum ncol
!#elif   SFA_ENABLED && SFB_ENABLED
!#else
!#error  "Unrecognized interface."
!#endif
        ! Declare constants and temporary variables.
#if     IK_ENABLED
#define TYPE_KIND integer(IKG)
        integer(IKG), parameter :: ZERO = 0_IKG, ONE = 1_IKG
#elif   CK_ENABLED
#define TYPE_KIND complex(CKG)
        complex(CKG), parameter :: ZERO = (0._CKG, 0._CKG), ONE = (1._CKG, 0._CKG)
#elif   RK_ENABLED
#define TYPE_KIND real(RKG)
        real(RKG)   , parameter :: ZERO = 0._RKG, ONE = 1._RKG
#else
#error  "Unrecognized interface."
#endif

!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!#if     setMatMulAdd_ENABLED && axpby_ENABLED
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!        integer(IK) :: nell, iell, bell, cell
!        if (matA == ZERO) return
!        nell =
!#if     ASS_ENABLED
!        ! Code for both increments equal to 1
!        do iell = 1, nell
!            matC(iell) = matA * matB(iell) + matC(iell)
!        end do
!#elif   EXP_ENABLED
!        if (incB == 1 .and. incC == 1) then
!        else ! code for unequal increments or equal increments not equal to 1.
!            bell = 1
!            cell = 1
!            if (incB < 0) bell = (1 - nell) * incB + 1
!            if (incC < 0) cell = (1 - nell) * incC + 1
!            do iell = 1, nell
!                matC(cell) = matC(cell) + matA * matB(bell)
!                bell = bell + incB
!                cell = cell + incC
!            end do
!        end if
!#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     setMatMulAdd_ENABLED && gemv_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the implied shapes and perform bound checks.
#if     !(BLAS_ENABLED && DISPATCH_ENABLED)
        integer(IK) :: irow, icol, iell, bstart, cstart
#endif
#if     ASS_ENABLED
        integer(IK) :: nrow, ncol
        integer(IK) , parameter :: incB = 1_IK, incC = 1_IK
        DECLARE_SET_DEFAULT_ALPHA_BETA
        nrow = size(matA, 1, IK)
        ncol = size(matA, 2, IK)
#elif   EXP_ENABLED
        ! Ensure offsets and subset bounds are logical.
        CHECK_ASSERTION(__LINE__, all(0_IK <= [nrow, ncol]), \
        SK_"@setMatMulAdd(): The condition `all(0_IK <= [nrow, ncol])` must hold. nrow, ncol = "//getStr([nrow, ncol])) ! fpp
        CHECK_ASSERTION(__LINE__, all(0_IK <= 1_IK - lbound(matA, kind = IK)), \
        SK_"@setMatMulAdd(): The condition `all(0 <= [roffA, coffA])` must hold. roffA, coffA = "//getStr(1_IK - lbound(matA))) ! fpp
        ! Ensure subset of `matA` does not overflow the matrix bounds.
        CHECK_ASSERTION(__LINE__, all([roffA + nrow, coffA + ncol] <= shape(matA, kind = IK)), \
        SK_"@setMatMulAdd(): The condition `all([roffA + nrow, coffA + ncol] <= shape(matA))` must hold. roffA, coffA, nrow, ncol, shape(matA) = "\
        //getStr([roffA, coffA, nrow, ncol, shape(matA, IK)])) ! fpp
#else
#error  "Unrecognized interface."
#endif
        ! Check bounds.
#if     TNA_ENABLED
#define LENB ncol
#define LENC nrow
        CHECK_ASSERTION(__LINE__, all([ncol, nrow] == [size(matB(::abs(incB)), 1, IK), size(matC(::abs(incC)), 1, IK)]), \
        SK_"@setMatMulAdd(): The condition `all([nrow, ncol] == [size(matC(::abs(incC))), size(matB(::abs(incB)))])` must hold. [nrow, ncol], size(matC(::abs(incC))), size(matB(::abs(incB))) = "\
        //getStr([[nrow, ncol], size(matC(::abs(incC)), 1, IK), size(matB(::abs(incB)), 1, IK)])) ! fpp
#elif   TSA_ENABLED || THA_ENABLED
#define LENB nrow
#define LENC ncol
        CHECK_ASSERTION(__LINE__, all([nrow, ncol] == [size(matB(::abs(incB)), 1, IK), size(matC(::abs(incC)), 1, IK)]), \
        SK_"@setMatMulAdd(): The condition `all([nrow, ncol] == [size(matB(::abs(incB))), size(matC(::abs(incC)))])` must hold. [nrow, ncol], size(matB(::abs(incB))), size(matC(::abs(incC))) = "\
        //getStr([shape(matA, IK), [nrow, ncol], size(matB(::abs(incB)), 1, IK), size(matC(::abs(incC)), 1, IK)])) ! fpp
#endif
        ! Quick return if possible.
        if (nrow == 0_IK .or. ncol == 0_IK .or. (ALPHA == ZERO .and. BETA == ONE)) return
#if     BLAS_ENABLED && DISPATCH_ENABLED
        call blasGEMV(TMA, nrow, ncol, ALPHA, matA(1,1), size(matA, 1, IK), matB(1), incB, BETA, matC(1), incC)
#else
        ! Set `lenB` and `lenC`, the lengths of the vectors `matB` and `matC`, and set up the start points in `matB` and `matC`.
        if (0_IK < incB) then
            bstart = 1
        else
            bstart = 1 - (LENB - 1) * incB
        end if
        if (0_IK < incC) then
            cstart = 1
        else
            cstart = 1 - (LENC - 1) * incC
        end if
        ! Start the operations. in this version the elements of matA are accessed sequentially with one pass through matA.
        ! First form  matC := beta * matC.
        if (BETA /= ONE) then
            if (incC == 1_IK) then
                if (BETA == ZERO) then
                    matC(1 : LENC) = ZERO
                else
                    do concurrent(irow = 1 : LENC)
                        matC(irow) = BETA * matC(irow)
                    end do
                end if
            else
                iell = cstart
                if (BETA == ZERO) then
                    do irow = 1, LENC
                        matC(iell) = ZERO
                        iell = iell + incC
                    end do
                else
                    do irow = 1, LENC
                        matC(iell) = BETA * matC(iell)
                        iell = iell + incC
                    end do
                end if
            end if
        end if
        if (ALPHA == ZERO) return
#if     TNA_ENABLED
        ! Form matC := alpha * matA * matB + matC.
        block
            integer(IK) :: cell
            TYPE_KIND :: temp
            iell = bstart
            if (incC == 1_IK) then
                do icol = 1, ncol
                    temp = ALPHA * matB(iell)
                    do irow = 1, nrow
                        matC(irow) = matC(irow) + temp * matA(irow, icol)
                    end do
                    iell = iell + incB
                end do
            else
                do icol = 1, ncol
                    temp = ALPHA * matB(iell)
                    cell = cstart
                    do irow = 1, nrow
                        matC(cell) = matC(cell) + temp * matA(irow, icol)
                        cell = cell + incC
                    end do
                    iell = iell + incB
                end do
            end if
        end block
#elif   TSA_ENABLED || THA_ENABLED
        ! Form matC := alpha * matA**T * matB + matC or matC := alpha * matA**H * matB + matC.
        block
            TYPE_KIND :: temp
            iell = cstart
            if (incB == 1_IK) then
                do icol = 1, ncol
                    temp = ZERO
                    do irow = 1, nrow
                        temp = temp + CONJG_A(matA(irow, icol)) * matB(irow)
                    end do
                    matC(iell) = matC(iell) + ALPHA * temp
                    iell = iell + incC
                end do
            else
                do icol = 1, ncol
                    temp = ZERO
                    cstart = bstart
                    do irow = 1, nrow
                        temp = temp + CONJG_A(matA(irow, icol)) * matB(cstart)
                        cstart = cstart + incB
                    end do
                    matC(iell) = matC(iell) + ALPHA * temp
                    iell = iell + incC
                end do
            end if
        end block
#else
#error  "Unrecognized interface."
#endif
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatMulAdd_ENABLED && (spmv_ENABLED || hpmv_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     spmv_ENABLED
#define PMV blasSPMV
#elif   hpmv_ENABLED
#define PMV blasHPMV
#else
#error  "Unrecognized interface."
#endif
#if     !(BLAS_ENABLED && DISPATCH_ENABLED)
        integer(IK) :: ib, ic, jb, jc, bstart, cstart, irow, icol, idum, jdum
#endif
        ! Define the implied shapes and perform bound checks.
#if     ASS_ENABLED
        integer(IK) , parameter :: incB = 1_IK, incC = 1_IK
        integer(IK) :: ndim
        DECLARE_SET_DEFAULT_ALPHA_BETA
        ndim = size(matB, 1, IK)
        CHECK_ASSERTION(__LINE__, ndim == size(matC, kind = IK), \
        SK_"@setMatMulAdd(): The condition `size(matB) == size(matC)` must hold. size(matB), size(matC) = "\
        //getStr([ndim, size(matC, kind = IK)])) ! fpp
#elif   EXP_ENABLED
        CHECK_ASSERTION(__LINE__, 0_IK <= ndim, \
        SK_"@setMatMulAdd(): The condition `0 <= ndim` must hold. ndim = "//getStr(ndim)) ! fpp
        CHECK_ASSERTION(__LINE__, all(0_IK /= [incB, incC]), \
        SK_"@setMatMulAdd(): The condition `all(0 /= [incB, incC])` must hold. ndim, roffA, coffA, incB, incC = "//\
        getStr([incB, incC])) ! fpp
        CHECK_ASSERTION(__LINE__, abs(incB) * (ndim - 1) < size(matB, kind = IK), \
        SK_"@setMatMulAdd(): The condition `abs(incB) * (ndim - 1) < size(matB(:))` must hold. incB, size(matB) = "//\
        getStr([incB, size(matB, kind = IK)])) ! fpp
        CHECK_ASSERTION(__LINE__, abs(incC) * (ndim - 1) < size(matC, kind = IK), \
        SK_"@setMatMulAdd(): The condition `abs(incC) * (ndim - 1) < size(matC(:))` must hold. incC, size(matC) = "//\
        getStr([incC, size(matC, kind = IK)])) ! fpp
#else
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, ndim * (ndim + 1_IK) / 2_IK == size(matA, kind = IK), \
        SK_"@setMatMulAdd(): The condition `ndim * (ndim + 1) / 2 == size(matA)` must hold. ndim, size(matA) = "\
        //getStr([ndim, size(matA, kind = IK)])) ! fpp
        ! Quick return if possible.
        if (ndim == 0_IK .or. (ALPHA == ZERO .and. BETA == ONE)) return
#if     BLAS_ENABLED && DISPATCH_ENABLED && ((spmv_ENABLED && RK_ENABLED) || (hpmv_ENABLED && CK_ENABLED)) && SUA_ENABLED
        call PMV("U", ndim, ALPHA, matA(1), matB(1), incB, BETA, matC(1), incC)
#elif   BLAS_ENABLED && DISPATCH_ENABLED && ((spmv_ENABLED && RK_ENABLED) || (hpmv_ENABLED && CK_ENABLED)) && SLA_ENABLED
        call PMV("L", ndim, ALPHA, matA(1), matB(1), incB, BETA, matC(1), incC)
#else
        ! Set up the start points in matC.
        if (incC > 0_IK) then
            cstart = 1_IK
        else
            cstart = 1_IK - (ndim - 1_IK) * incC
        end if
        ! Start the operations. In this version the elements of the array matA
        ! are accessed sequentially with ONE pass through matA.
        ! First Form  matC := BETA * matC.
        if (BETA /= ONE) then
            if (incC == 1_IK) then
                if (BETA == ZERO) then
                    do irow = 1_IK, ndim
                        matC(irow) = ZERO
                    end do
                else
                    do irow = 1_IK, ndim
                        matC(irow) = BETA * matC(irow)
                    end do
                end if
            else
                ic = cstart
                if (BETA == ZERO) then
                    do irow = 1_IK, ndim
                        matC(ic) = ZERO
                        ic = ic + incC
                    end do
                else
                    do irow = 1_IK, ndim
                        matC(ic) = BETA * matC(ic)
                        ic = ic + incC
                    end do
                end if
            end if
        end if
        if (ALPHA == ZERO) return
        ! Set up the start points in matB.
        if (incB > 0_IK) then
            bstart = 1_IK
        else
            bstart = 1_IK - (ndim - 1_IK) * incB
        end if
        jdum = 1_IK
#if     SUA_ENABLED
        ! Form matC when matA contains the upper triangle.
        block
            TYPE_KIND :: temp, dumm
            if (incB == 1_IK .and. incC == 1_IK) then
                do icol = 1_IK, ndim
                    temp = ALPHA * matB(icol)
                    dumm = ZERO
                    idum = jdum
                    do irow = 1_IK,icol - 1_IK
                        matC(irow) = matC(irow) + temp * matA(idum)
                        dumm = dumm + CONJG_A(matA(idum)) * matB(irow)
                        idum = idum + 1_IK
                    end do
                    matC(icol) = matC(icol) + temp * GET_RE(matA(jdum + icol - 1)) + ALPHA * dumm
                    jdum = jdum + icol
                end do
            else
                jb = bstart
                jc = cstart
                do icol = 1_IK, ndim
                    temp = ALPHA * matB(jb)
                    dumm = ZERO
                    ib = bstart
                    ic = cstart
                    do idum = jdum, jdum + icol - 2
                        matC(ic) = matC(ic) + temp * matA(idum)
                        dumm = dumm + CONJG_A(matA(idum)) * matB(ib)
                        ib = ib + incB
                        ic = ic + incC
                    end do
                    matC(jc) = matC(jc) + temp * GET_RE(matA(jdum + icol - 1)) + ALPHA * dumm
                    jb = jb + incB
                    jc = jc + incC
                    jdum = jdum + icol
                end do
            end if
        end block
#elif   SLA_ENABLED
        ! Form matC when matA contains the lower triangle.
        block
            TYPE_KIND :: temp, dumm
            if (incB == 1_IK .and. incC == 1_IK) then
                do icol = 1_IK, ndim
                    temp = ALPHA * matB(icol)
                    dumm = ZERO
                    matC(icol) = matC(icol) + temp * GET_RE(matA(jdum))
                    idum = jdum + 1_IK
                    do irow = icol + 1_IK, ndim
                        matC(irow) = matC(irow) + temp * matA(idum)
                        dumm = dumm + CONJG_A(matA(idum)) * matB(irow)
                        idum = idum + 1_IK
                    end do
                    matC(icol) = matC(icol) + ALPHA * dumm
                    jdum = jdum + (ndim - icol + 1)
                end do
            else
                jb = bstart
                jc = cstart
                do icol = 1_IK, ndim
                    temp = ALPHA * matB(jb)
                    dumm = ZERO
                    matC(jc) = matC(jc) + temp * GET_RE(matA(jdum))
                    ib = jb
                    ic = jc
                    do idum = jdum + 1_IK, jdum + ndim - icol
                        ib = ib + incB
                        ic = ic + incC
                        matC(ic) = matC(ic) + temp * matA(idum)
                        dumm = dumm + CONJG_A(matA(idum)) * matB(ib)
                    end do
                    matC(jc) = matC(jc) + ALPHA * dumm
                    jb = jb + incB
                    jc = jc + incC
                    jdum = jdum + (ndim - icol + 1)
                end do
            end if
        end block
#else
#error  "Unrecognized interface."
#endif
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatMulAdd_ENABLED && (symv_ENABLED || hemv_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     symv_ENABLED
#define PMV blasSYMV
#elif   hemv_ENABLED
#define PMV blasHEMV
#else
#error  "Unrecognized interface."
#endif
#if     !(BLAS_ENABLED && DISPATCH_ENABLED && ((spmv_ENABLED && RK_ENABLED) || (hpmv_ENABLED && CK_ENABLED)) && (SLA_ENABLED || SUA_ENABLED))
        integer(IK) :: ib, ic, jb, jc, bstart, cstart, irow, icol
#endif
        ! Define the implied shapes and perform bound checks.
#if     ASS_ENABLED
        integer(IK) , parameter :: incB = 1_IK, incC = 1_IK
        integer(IK) :: ndim
        DECLARE_SET_DEFAULT_ALPHA_BETA
        ndim = size(matA, 1, IK)
        CHECK_ASSERTION(__LINE__, all(ndim == [size(matA, 2, IK), shape(matB, IK), shape(matC, IK)]), \
        SK_"@setMatMulAdd(): The condition `all(ndim == [shape(matA(:,:)), shape(matB(:)), shape(matC(:))])` must hold. ndim, shape(matA), shape(matB), shape(matC) = "\
        //getStr([ndim, shape(matA, IK), shape(matB, IK), shape(matC, IK)])) ! fpp
#elif   EXP_ENABLED
        CHECK_ASSERTION(__LINE__, all(0_IK <= [ndim, roffA, coffA]), \
        SK_"@setMatMulAdd(): The condition `all(0 <= [ndim, roffA, coffA])` must hold. ndim, roffA, coffA = "//\
        getStr([ndim, roffA, coffA])) ! fpp
        CHECK_ASSERTION(__LINE__, all(0_IK /= [incB, incC]), \
        SK_"@setMatMulAdd(): The condition `all(0 /= [incB, incC])` must hold. ndim, roffA, coffA, incB, incC = "//\
        getStr([incB, incC])) ! fpp
        CHECK_ASSERTION(__LINE__, all(ndim <= ubound(matA, kind = IK)), \
        SK_"@setMatMulAdd(): The condition `all(ndim + [roffA, coffA] <= shape(matA))` must hold. ndim, ubound(matA) = "//\
        getStr([ndim, ubound(matA, kind = IK)])) ! fpp
        CHECK_ASSERTION(__LINE__, abs(incB) * (ndim - 1) < size(matB, kind = IK), \
        SK_"@setMatMulAdd(): The condition `abs(incB) * (ndim - 1) < size(matB(:))` must hold. incB, size(matB) = "//\
        getStr([incB, size(matB, kind = IK)])) ! fpp
        CHECK_ASSERTION(__LINE__, abs(incC) * (ndim - 1) < size(matC, kind = IK), \
        SK_"@setMatMulAdd(): The condition `abs(incC) * (ndim - 1) < size(matC(:))` must hold. incC, size(matC) = "//\
        getStr([incC, size(matC, kind = IK)])) ! fpp
#else
#error  "Unrecognized interface."
#endif
        if (ndim == 0_IK .or. (ALPHA == ZERO .and. BETA == ONE)) return ! Quick return if possible.
#if     BLAS_ENABLED && DISPATCH_ENABLED && ((spmv_ENABLED && RK_ENABLED) || (hpmv_ENABLED && CK_ENABLED)) && SUA_ENABLED
        call PMV("U", ndim, ALPHA, matA(1,1), size(matA, 1, IK), matB(1), incB, BETA, matC(1), incC)
#elif   BLAS_ENABLED && DISPATCH_ENABLED && ((spmv_ENABLED && RK_ENABLED) || (hpmv_ENABLED && CK_ENABLED)) && SLA_ENABLED
        call PMV("L", ndim, ALPHA, matA(1,1), size(matA, 1, IK), matB(1), incB, BETA, matC(1), incC)
#else
        if(incC > 0_IK) then
            cstart = 1_IK
        else
            cstart = 1_IK - (ndim - 1_IK) * incC
        end if
        ! Start the operations. In this version the elements of matA are
        ! accessed sequentially with ONE pass through the triangular part of matA.
        ! First form  matC := BETA * matC.
        if(BETA /= ONE) then
            if(incC == 1_IK) then
                if (BETA == ZERO) then
                    do irow = 1_IK, ndim
                        matC(irow) = ZERO
                    end do
                else
                    do irow = 1_IK, ndim
                        matC(irow) = BETA * matC(irow)
                    end do
                end if
            else
                ic = cstart
                if(BETA == ZERO) then
                    do irow = 1_IK, ndim
                        matC(ic) = ZERO
                        ic = ic + incC
                    end do
                else
                    do irow = 1_IK, ndim
                        matC(ic) = BETA * matC(ic)
                        ic = ic + incC
                    end do
                end if
            end if
        end if
        if(ALPHA == ZERO) return
        ! set up the start points in matB and matC.
        if(incB > 0_IK) then
            bstart = 1_IK
        else
            bstart = 1_IK - (ndim - 1_IK) * incB
        end if
#if     SUA_ENABLED
        ! Form  matC when matA is stored in upper triangle.
        block
            TYPE_KIND :: temp, dumm
            if(incB == 1_IK .and. incC == 1_IK) then
                do icol = 1_IK, ndim
                    temp = ALPHA * matB(icol)
                    dumm = ZERO
                    do irow = 1_IK, icol - 1_IK
                        matC(irow) = matC(irow) + temp * matA(irow, icol)
                        dumm = dumm + CONJG_A(matA(irow, icol)) * matB(irow)
                    end do
                    matC(icol) = matC(icol) + temp * GET_RE(matA(icol, icol)) + ALPHA * dumm
                end do
            else
                jb = bstart
                jc = cstart
                do icol = 1_IK, ndim
                    temp = ALPHA * matB(jb)
                    dumm = ZERO
                    ib = bstart
                    ic = cstart
                    do irow = 1_IK, icol - 1_IK
                        matC(ic) = matC(ic) + temp * matA(irow, icol)
                        dumm = dumm + CONJG_A(matA(irow, icol)) * matB(ib)
                        ib = ib + incB
                        ic = ic + incC
                    end do
                    matC(jc) = matC(jc) + temp * GET_RE(matA(icol, icol)) + ALPHA * dumm
                    jb = jb + incB
                    jc = jc + incC
                end do
            end if
        end block
#elif   SLA_ENABLED
        ! Form matC when matA is stored in lower triangle.
        block
            TYPE_KIND :: temp, dumm
            if (incB == 1_IK .and. incC == 1_IK) then
                do icol = 1_IK, ndim
                    temp = ALPHA * matB(icol)
                    dumm = ZERO
                    matC(icol) = matC(icol) + temp * GET_RE(matA(icol, icol))
                    do irow = icol + 1_IK, ndim
                        matC(irow) = matC(irow) + temp * matA(irow, icol)
                        dumm = dumm + CONJG_A(matA(irow, icol)) * matB(irow)
                    end do
                    matC(icol) = matC(icol) + ALPHA * dumm
                end do
            else
                jb = bstart
                jc = cstart
                do icol = 1_IK, ndim
                    temp = ALPHA * matB(jb)
                    dumm = ZERO
                    matC(jc) = matC(jc) + temp * GET_RE(matA(icol, icol))
                    ib = jb
                    ic = jc
                    do irow = icol + 1_IK, ndim
                        ib = ib + incB
                        ic = ic + incC
                        matC(ic) = matC(ic) + temp * matA(irow, icol)
                        dumm = dumm + CONJG_A(matA(irow, icol)) * matB(ib)
                    end do
                    matC(jc) = matC(jc) + ALPHA * dumm
                    jb = jb + incB
                    jc = jc + incC
                end do
            end if
        end block
#else
#error  "Unrecognized interface."
#endif
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatMulAdd_ENABLED && (symm_ENABLED || hemm_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     !(BLAS_ENABLED && DISPATCH_ENABLED)
        integer(IK) :: idum, irow, icol
#endif
        ! Define the implied shapes and perform bound checks.
#if     ASS_ENABLED
        integer(IK) :: nrow, ncol
        DECLARE_SET_DEFAULT_ALPHA_BETA
        nrow = size(matC, 1, IK)
        ncol = size(matC, 2, IK)
        ! Check bounds.
#if     CNA_ENABLED
        CHECK_ASSERTION(__LINE__, all(shape(matA, IK) == [nrow, ncol]), \
        SK_"@setMatMulAdd(): The condition `all(shape(matA) == shape(matC))` must hold. shape(matA), nrow, ncol = "\
        //getStr([shape(matA, IK), nrow, ncol])) ! fpp
#else
        CHECK_ASSERTION(__LINE__, all(shape(matA, IK) == [nrow, nrow]), \
        SK_"@setMatMulAdd(): The condition `all(shape(matA) == size(matC, 1))` must hold. shape(matA), size(matC, 1) = "\
        //getStr([shape(matA, IK), nrow])) ! fpp
#endif
#if     CNB_ENABLED
        CHECK_ASSERTION(__LINE__, all(shape(matB, IK) == [nrow, ncol]), \
        SK_"@setMatMulAdd(): The condition `all(shape(matB) == shape(matC))` must hold. shape(matB), shape(matC) = "\
        //getStr([shape(matB, IK), nrow, ncol])) ! fpp
#else
        CHECK_ASSERTION(__LINE__, all(shape(matB, IK) == [ncol, ncol]), \
        SK_"@setMatMulAdd(): The condition `all(shape(matB) == size(matC, 2))` must hold. shape(matB), size(matC, 2) = "\
        //getStr([shape(matB, IK), ncol])) ! fpp
#endif
#elif   EXP_ENABLED
        ! Ensure offsets and subset bounds are logical.
        CHECK_ASSERTION(__LINE__, all(0_IK <= [nrow, ncol]), \
        SK_"@setMatMulAdd(): The condition `all(0_IK <= shape(matC))` must hold. shape(matC) = "//getStr([nrow, ncol])) ! fpp
        CHECK_ASSERTION(__LINE__, all(0_IK <= 1_IK - lbound(matA, kind = IK)), \
        SK_"@setMatMulAdd(): The condition `all(0 <= [roffA, coffA])` must hold. roffA, coffA = "//getStr(1_IK - lbound(matA))) ! fpp
        CHECK_ASSERTION(__LINE__, all(0_IK <= 1_IK - lbound(matB, kind = IK)), \
        SK_"@setMatMulAdd(): The condition `all(0 <= [roffB, coffB])` must hold. roffB, coffB = "//getStr(1_IK - lbound(matB))) ! fpp
        CHECK_ASSERTION(__LINE__, all(0_IK <= 1_IK - lbound(matC, kind = IK)), \
        SK_"@setMatMulAdd(): The condition `all(0 <= [roffC, coffC])` must hold. roffC, coffC = "//getStr(1_IK - lbound(matC))) ! fpp
        ! Ensure subset of `matA` does not overflow the matrix bounds.
#if     CNA_ENABLED
        CHECK_ASSERTION(__LINE__, all([nrow, ncol] <= ubound(matA, kind = IK)), \
        SK_"@setMatMulAdd(): The condition `all([nrow + roffA, ncol + coffA] <= shape(matA))` must hold. nrow, ncol, roffA, coffA, shape(matA) = "\
        //getStr([nrow, ncol, roffA, coffA, shape(matA, IK)])) ! fpp
#else
        CHECK_ASSERTION(__LINE__, all([nrow, nrow] <= ubound(matA, kind = IK)), \
        SK_"@setMatMulAdd(): The condition `all([nrow + roffA, nrow + coffA] <= shape(matA))` must hold. nrow, roffA, coffA, shape(matA) = "\
        //getStr([nrow, roffA, coffA, shape(matA, IK)])) ! fpp
#endif
        ! Ensure subset of `matB` does not overflow the matrix bounds.
#if     CNB_ENABLED
        CHECK_ASSERTION(__LINE__, all([nrow, ncol] <= ubound(matB, kind = IK)), \
        SK_"@setMatMulAdd(): The condition `all([nrow + roffB, ncol + coffB] <= shape(matB))` must hold. nrow, ncol, roffB, coffB, shape(matB) = "\
        //getStr([nrow, ncol, roffB, coffB, shape(matB, IK)])) ! fpp
#else
        CHECK_ASSERTION(__LINE__, all([ncol, ncol] <= ubound(matB, kind = IK)), \
        SK_"@setMatMulAdd(): The condition `all([ncol + roffB, ncol + coffB] <= shape(matB))` must hold. ncol, roffB, coffB, shape(matB) = "\
        //getStr([ncol, roffB, coffB, shape(matB, IK)])) ! fpp
#endif
#else
#error  "Unrecognized interface."
#endif
        ! Quick return if possible.
        if (nrow == 0_IK .or. ncol == 0_IK .or. (ALPHA == ZERO .and. BETA == ONE)) return
        if (ALPHA == ZERO) then
            if (BETA == ZERO) then
                do concurrent(icol = 1 : ncol, irow = 1 : nrow)
                    matC(irow, icol) = ZERO
                end do
            else
                do concurrent(icol = 1 : ncol, irow = 1 : nrow)
                    matC(irow, icol) = BETA * matC(irow, icol)
                end do
            end if
            return
        end if
        ! BLAS 3 SYMM/HEMM: Form matC := alpha * matA(upper-triangular) * matB + beta * matC.
#if     SUA_ENABLED && SFB_ENABLED && (CSA_ENABLED || CHA_ENABLED) && CNB_ENABLED && BLAS_ENABLED && DISPATCH_ENABLED
        call SHMM("L", "U", nrow, ncol, ALPHA, matA, size(matA, 1, IK), matB, size(matB, 1, IK), BETA, matC, size(matC, 1, IK))
#elif   SUA_ENABLED && SFB_ENABLED && (CSA_ENABLED || CHA_ENABLED) && CNB_ENABLED
        block
            TYPE_KIND :: temp, dumm
            do icol = 1, ncol
                do irow = 1, nrow
                        temp = ALPHA * matB(irow, icol)
                        dumm = ZERO
                        do idum = 1, irow - 1
                            matC(idum, icol) = matC(idum, icol) + temp * matA(idum, irow)
                            dumm = dumm + matB(idum, icol) * CONJG_A(matA(idum, irow))
                        end do
                    if (BETA == ZERO) then
                        matC(irow, icol) = temp * GET_RE(matA(irow, irow)) + ALPHA * dumm
                    else
                        matC(irow, icol) = BETA * matC(irow, icol) + temp * GET_RE(matA(irow, irow)) + ALPHA * dumm
                    end if
                end do
            end do
        end block
        ! BLAS 3 SYMM/HEMM: Form matC := alpha * matA(lower-triangular) * matB + beta * matC.
#elif   SLA_ENABLED && SFB_ENABLED && (CSA_ENABLED || CHA_ENABLED) && CNB_ENABLED && BLAS_ENABLED && DISPATCH_ENABLED
        call SHMM("L", "L", nrow, ncol, ALPHA, matA, size(matA, 1, IK), matB, size(matB, 1, IK), BETA, matC, size(matC, 1, IK))
#elif   SLA_ENABLED && SFB_ENABLED && (CSA_ENABLED || CHA_ENABLED) && CNB_ENABLED
        block
            TYPE_KIND :: temp, dumm
            do icol = 1, ncol
                do irow = nrow, 1, -1
                    temp = ALPHA * matB(irow, icol)
                    dumm = ZERO
                    do idum = irow + 1, nrow
                        matC(idum, icol) = matC(idum, icol) + temp * matA(idum, irow)
                        dumm = dumm + matB(idum, icol) * CONJG_A(matA(idum, irow))
                    end do
                    if (BETA == ZERO) then
                        matC(irow, icol) = temp * GET_RE(matA(irow, irow)) + ALPHA * dumm
                    else
                        matC(irow, icol) = BETA * matC(irow, icol) + temp * GET_RE(matA(irow, irow)) + ALPHA * dumm
                    end if
                end do
            end do
        end block
        ! BLAS 3 SYMM/HEMM: Form matC := alpha * matA * matB(triangular) + beta * matC.
#elif   SFA_ENABLED && SUB_ENABLED && CNA_ENABLED && (CSB_ENABLED || CHB_ENABLED) && BLAS_ENABLED && DISPATCH_ENABLED
        call SHMM("R", "U", ncol, nrow, ALPHA, matB, size(matB, 1, IK), matA, size(matA, 1, IK), BETA, matC, size(matC, 1, IK))
#elif   SFA_ENABLED && SLB_ENABLED && CNA_ENABLED && (CSB_ENABLED || CHB_ENABLED) && BLAS_ENABLED && DISPATCH_ENABLED
        call SHMM("R", "L", ncol, nrow, ALPHA, matB, size(matB, 1, IK), matA, size(matA, 1, IK), BETA, matC, size(matC, 1, IK))
#elif   SFA_ENABLED && (SUB_ENABLED || SLB_ENABLED) && CNA_ENABLED && (CSB_ENABLED || CHB_ENABLED)
        block
            TYPE_KIND :: temp
            do icol = 1, ncol
                temp = ALPHA * GET_RE(matB(icol, icol))
                if (BETA == ZERO) then
                    do irow = 1, nrow
                        matC(irow, icol) = temp * matA(irow, icol)
                    end do
                else
                    do irow = 1, nrow
                        matC(irow, icol) = BETA * matC(irow, icol) + temp * matA(irow, icol)
                    end do
                end if
                do idum = 1, icol - 1
#if                 SUB_ENABLED
                    temp = ALPHA * matB(idum, icol)
#elif               SLB_ENABLED
                    temp = ALPHA * CONJG_B(matB(icol, idum))
#else
#error              "Unrecognized interface."
#endif
                    do irow = 1, nrow
                        matC(irow, icol) = matC(irow, icol) + temp * matA(irow, idum)
                    end do
                end do
                do idum = icol + 1, ncol
#if                 SUB_ENABLED
                    temp = ALPHA * CONJG_B(matB(icol, idum))
#elif               SLB_ENABLED
                    temp = ALPHA * matB(idum, icol)
#else
#error              "Unrecognized interface."
#endif
                    do irow = 1, nrow
                        matC(irow, icol) = matC(irow, icol) + temp * matA(irow, idum)
                    end do
                end do
            end do
        end block
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatMulAdd_ENABLED && gemm_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: irow, icol
#if     !(BLAS_ENABLED && DISPATCH_ENABLED)
        integer(IK) :: idum
#endif
        ! Define the implied shapes and perform bound checks.
#if     ASS_ENABLED
        integer(IK) :: nrow, ncol, ndum
        DECLARE_SET_DEFAULT_ALPHA_BETA
        nrow = size(matC, 1, IK)
        ncol = size(matC, 2, IK)
        ! Set the length of the dummy axis of `matA` but only if it is a general asymmetric matrix.
#if     SFA_ENABLED && TNA_ENABLED
        ndum = size(matA, 2, IK)
#elif   SFA_ENABLED && (TSA_ENABLED || THA_ENABLED)
        ndum = size(matA, 1, IK)
#else
#error  "Unrecognized interface."
#endif
        ! Check bounds.
#if     TNA_ENABLED
        CHECK_ASSERTION(__LINE__, all(shape(matA, IK) == [nrow, ndum]), \
        SK_"@setMatMulAdd(): The condition `all(shape(matA) == [nrow, ndum])` must hold. shape(matA), nrow, ndum = "\
        //getStr([shape(matA, IK), nrow, ndum])) ! fpp
#else
        CHECK_ASSERTION(__LINE__, all(shape(matA, IK) == [ndum, nrow]), \
        SK_"@setMatMulAdd(): The condition `all(shape(matA) == [ndum, nrow])` must hold. shape(matA), ndum, nrow = "\
        //getStr([shape(matA, IK), ndum, nrow])) ! fpp
#endif
#if     TNB_ENABLED
        CHECK_ASSERTION(__LINE__, all(shape(matB, IK) == [ndum, ncol]), \
        SK_"@setMatMulAdd(): The condition `all(shape(matB) == [ndum, ncol])` must hold. shape(matB), ndum, ncol = "\
        //getStr([shape(matB, IK), ndum, ncol])) ! fpp
#else
        CHECK_ASSERTION(__LINE__, all(shape(matB, IK) == [ncol, ndum]), \
        SK_"@setMatMulAdd(): The condition `all(shape(matB) == [ncol, ndum])` must hold. shape(matB), ncol, ndum = "\
        //getStr([shape(matB, IK), ncol, ndum])) ! fpp
#endif
#elif   EXP_ENABLED
        ! Ensure offsets and subset bounds are logical.
        CHECK_ASSERTION(__LINE__, all(0_IK <= [nrow, ncol, ndum]), \
        SK_"@setMatMulAdd(): The condition `all(0_IK <= [nrow, ncol, ndum])` must hold. nrow, ncol, ndum = "//getStr([nrow, ncol, ndum])) ! fpp
        CHECK_ASSERTION(__LINE__, all(0_IK <= 1_IK - lbound(matA, kind = IK)), \
        SK_"@setMatMulAdd(): The condition `all(0 <= [roffA, coffA])` must hold. roffA, coffA = "//getStr(1_IK - lbound(matA))) ! fpp
        CHECK_ASSERTION(__LINE__, all(0_IK <= 1_IK - lbound(matB, kind = IK)), \
        SK_"@setMatMulAdd(): The condition `all(0 <= [roffB, coffB])` must hold. roffB, coffB = "//getStr(1_IK - lbound(matB))) ! fpp
        CHECK_ASSERTION(__LINE__, all(0_IK <= 1_IK - lbound(matC, kind = IK)), \
        SK_"@setMatMulAdd(): The condition `all(0 <= [roffC, coffC])` must hold. roffC, coffC = "//getStr(1_IK - lbound(matC))) ! fpp
        ! Ensure subset of `matA` does not overflow the matrix bounds.
#if     TNA_ENABLED
        CHECK_ASSERTION(__LINE__, all([roffA + nrow, coffA + ndum] <= shape(matA, kind = IK)), \
        SK_"@setMatMulAdd(): The condition `all([roffA + nrow, coffA + ndum] <= shape(matA))` must hold. roffA, coffA, nrow, ndum, shape(matA) = "\
        //getStr([roffA, coffA, nrow, ndum, shape(matA, IK)])) ! fpp
#elif   TSA_ENABLED || THA_ENABLED
        CHECK_ASSERTION(__LINE__, all([roffA + ndum, coffA + nrow] <= shape(matA, IK)), \
        SK_"@setMatMulAdd(): The condition `all([roffA + ndum, coffA + nrow] <= shape(matA))` must hold. roffA, coffA, nrow, ndum, shape(matA) = "\
        //getStr([roffA, coffA, nrow, ndum, shape(matA, IK)])) ! fpp
#else
#error  "Unrecognized interface."
#endif
        ! Ensure subset of `matB` does not overflow the matrix bounds.
#if     TNB_ENABLED
        CHECK_ASSERTION(__LINE__, all([roffB + ndum, coffB + ncol] <= shape(matB, IK)), \
        SK_"@setMatMulAdd(): The condition `all([roffB + ndum, coffB + ncol] <= shape(matB))` must hold. roffB, coffB, ncol, ndum, shape(matB) = "\
        //getStr([roffB, coffB, ncol, ndum, shape(matB, IK)])) ! fpp
#elif   TSB_ENABLED || THB_ENABLED
        CHECK_ASSERTION(__LINE__, all([roffB + ncol, coffB + ndum] <= shape(matB, IK)), \
        SK_"@setMatMulAdd(): The condition `all([roffB + ncol, coffB + ndum] <= shape(matB))` must hold. roffB, coffB, ncol, ndum, shape(matB) = "\
        //getStr([roffB, coffB, ncol, ndum, shape(matB, IK)])) ! fpp
#else
#error  "Unrecognized interface."
#endif
#else
#error  "Unrecognized interface."
#endif
        ! Quick return if possible.
        if (nrow == 0_IK .or. ncol == 0_IK .or. ((ALPHA == ZERO .or. ndum == 0_IK) .and. BETA == ONE)) return
        if (ALPHA == ZERO .or. ndum == 0_IK) then ! The condition `.or. ndum == 0_IK` is extra to BLAS here, to bypass the implicit-interface limitation of BLAS.
            if (BETA == ZERO) then
                do concurrent(icol = 1 : ncol, irow = 1 : nrow)
                    matC(irow, icol) = ZERO
                end do
            else
                do concurrent(icol = 1 : ncol, irow = 1 : nrow)
                    matC(irow, icol) = BETA * matC(irow, icol)
                end do
            end if
            return
        end if
        ! BLAS 3 GEMM: Form matC := alpha * matA * matB + beta * matC.
        ! \todo There is a lot of room for improvement here, particularly when order of matC is large.
        ! For example, the BETA value check should be completely taken out of the icol loop.
#if     SFA_ENABLED && SFB_ENABLED && TNA_ENABLED && TNB_ENABLED && BLAS_ENABLED && DISPATCH_ENABLED
        call GEMM(TMA, TMB, ALPHA, BETA)
#elif   SFA_ENABLED && SFB_ENABLED && TNA_ENABLED && TNB_ENABLED
        block
            TYPE_KIND :: temp
#if         0
            ! assume BETA /= ZERO .and. BETA /= ONE
            do icol = 1, ncol
                idum = 1
                temp = ALPHA * matB(idum, icol)
                do irow = 1, nrow
                    matC(irow, icol) = BETA * matC(irow, icol) + temp * matA(irow, idum)
                end do
                do idum = 2, ndum
                    temp = ALPHA * matB(idum, icol)
                    do irow = 1, nrow
                        matC(irow, icol) = matC(irow, icol) + temp * matA(irow, idum)
                    end do
                end do
            end do
#else
            do icol = 1, ncol
                if (BETA == ZERO) then
                    matC(1 : nrow, icol) = ZERO
                elseif (BETA /= ONE) then
                    do irow = 1, nrow
                        matC(irow, icol) = BETA * matC(irow, icol)
                    end do
                end if
                do idum = 1, ndum
                    temp = ALPHA * matB(idum, icol)
                    do irow = 1, nrow
                        matC(irow, icol) = matC(irow, icol) + temp * matA(irow, idum)
                    end do
                end do
            end do
        end block
#endif
        ! BLAS 3 GEMM: Form matC := alpha * transpose(matA) * matB + beta * matC
#elif   SFA_ENABLED && SFB_ENABLED && (TSA_ENABLED || THA_ENABLED) && TNB_ENABLED && BLAS_ENABLED && DISPATCH_ENABLED
        call GEMM(TMA, TMB, ALPHA, BETA)
#elif   SFA_ENABLED && SFB_ENABLED && (TSA_ENABLED || THA_ENABLED) && TNB_ENABLED
        block
            TYPE_KIND :: temp
            do icol = 1, ncol
                do irow = 1, nrow
                    temp = ZERO
                    do idum = 1, ndum
                        temp = temp + CONJG_A(matA(idum, irow)) * matB(idum, icol)
                    end do
                    if (BETA == ZERO) then
                        matC(irow, icol) = ALPHA * temp
                    else
                        matC(irow, icol) = ALPHA * temp + BETA * matC(irow, icol)
                    end if
                end do
            end do
        end block
        ! BLAS 3 GEMM: Form matC := alpha * matA * transpose(matB) + beta * matC
#elif   SFA_ENABLED && SFB_ENABLED && TNA_ENABLED && (TSB_ENABLED || THB_ENABLED) && BLAS_ENABLED && DISPATCH_ENABLED
        call GEMM(TMA, TMB, ALPHA, BETA)
#elif   SFA_ENABLED && SFB_ENABLED && TNA_ENABLED && (TSB_ENABLED || THB_ENABLED)
        block
            TYPE_KIND :: temp
            do icol = 1, ncol
                if (BETA == ZERO) then
                    matC(1 : nrow, icol) = ZERO
                else if (BETA /= ONE) then
                    do irow = 1, nrow
                        matC(irow, icol) = BETA * matC(irow, icol)
                    end do
                end if
                do idum = 1, ndum
                    temp = ALPHA * CONJG_B(matB(icol, idum))
                    do irow = 1, nrow
                        matC(irow, icol) = matC(irow, icol) + temp * matA(irow, idum)
                    end do
                end do
            end do
        end block
        ! BLAS 3 GEMM: Form matC := alpha * transpose(matA) * transpose(matB) + beta * matC
#elif   SFA_ENABLED && SFB_ENABLED && (TSA_ENABLED || THA_ENABLED) && (TSB_ENABLED || THB_ENABLED) && BLAS_ENABLED && DISPATCH_ENABLED
        call GEMM(TMA, TMB, ALPHA, BETA)
#elif   SFA_ENABLED && SFB_ENABLED && (TSA_ENABLED || THA_ENABLED) && (TSB_ENABLED || THB_ENABLED)
        block
            TYPE_KIND :: temp
            do icol = 1, ncol
                do irow = 1, nrow
                    temp = ZERO
                    do idum = 1, ndum
                        temp = temp + CONJG_A(matA(idum, irow)) * CONJG_B(matB(icol, idum))
                    end do
                    if (BETA == ZERO) then
                        matC(irow, icol) = ALPHA * temp
                    else
                        matC(irow, icol) = ALPHA * temp + BETA * matC(irow, icol)
                    end if
                end do
            end do
        end block
#else
#error  "Unrecognized interface."
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  DECLARE_SET_DEFAULT_ALPHA_BETA
#undef  TYPE_KIND
#undef  DECLARE
#undef  CONJG_A
#undef  CONJG_B
#undef  GET_RE
#undef  ALPHA
#undef  BETA
#undef  ndum
#undef  SHMM
#undef  GEMM
#undef  LENB
#undef  LENC
#undef  TMA
#undef  TMB
#undef  PMV