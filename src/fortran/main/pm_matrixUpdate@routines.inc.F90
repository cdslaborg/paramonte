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
!>  This include file contains procedure implementation of the generic interfaces of [pm_matrixUpdate](@ref pm_matrixUpdate).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the BLAS storage triangle.
#if     shrk_ENABLED && SLD_ENABLED
#define BLAS_UPLO "L"
#elif   shrk_ENABLED && SUD_ENABLED
#define BLAS_UPLO "U"
#elif   !(setMatUpdateR1_ENABLED || setMatUpdateTriang_ENABLED)
#error  "Unrecognized interface."
#endif
        ! Define the BLAS storage triangle.
#if     shrk_ENABLED && OTP_ENABLED && CHM_ENABLED && CK_ENABLED
#define blasXXRK blasHERK
#define BLAS_TRANS "C"
#elif   shrk_ENABLED && ONO_ENABLED && CHM_ENABLED && CK_ENABLED
#define blasXXRK blasHERK
#define BLAS_TRANS "N"
#elif   shrk_ENABLED && OTP_ENABLED && (CHM_ENABLED || CSM_ENABLED) && (CK_ENABLED || RK_ENABLED)
#define blasXXRK blasSYRK
#define BLAS_TRANS "T"
#elif   shrk_ENABLED && ONO_ENABLED && (CHM_ENABLED || CSM_ENABLED) && (CK_ENABLED || RK_ENABLED)
#define blasXXRK blasSYRK
#define BLAS_TRANS "N"
#elif   !(IK_ENABLED || setMatUpdateR1_ENABLED || setMatUpdateTriang_ENABLED)
#error  "Unrecognized interface."
#endif
        ! Define the constants.
#if     IK_ENABLED
#define TYPE_KIND integer(TKG)
        integer(TKG), parameter :: ZERO = 0_TKG, ONE = 1_TKG
#elif   CK_ENABLED
#define TYPE_KIND complex(TKG)
        complex(TKG), parameter :: ZERO = (0._TKG, 0._TKG), ONE = (1._TKG, 0._TKG)
#elif   RK_ENABLED
#define TYPE_KIND real(TKG)
        real(TKG)   , parameter :: ZERO = 0._TKG, ONE = 1._TKG
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     setMatUpdate_ENABLED && shrk_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     !(BLAS_ENABLED && DISPATCH_ENABLED)
        TYPE_KIND :: temp
        integer(IK) :: irow, icol, idum
#endif
#if     ASS_ENABLED
#define BETA beta_def
#define ALPHA alpha_def
        integer(IK) :: ndim, ndum
        TYPE_KIND :: alpha_def, beta_def
        if (present(beta)) then; beta_def = beta; else; beta_def = ONE; end if
        if (present(alpha)) then; alpha_def = alpha; else; alpha_def = ONE; end if
        ndim = size(mat, 1, IK)
        CHECK_ASSERTION(__LINE__, size(mat, 1, IK) == size(mat, 2, IK), SK_"@setMatUpdate(): The condition `size(mat, 1) == size(mat, 2)` must hold. shape(mat) = "//getStr(shape(mat, IK)))
#if     ONO_ENABLED
        ndum = size(matA, 2, IK)
        CHECK_ASSERTION(__LINE__, size(mat, 2, IK) == size(matA, 1, IK), SK_"@setMatUpdate(): The condition `size(mat, 2) == size(matA, 1)` must hold. shape(mat), shape(matA) = "//getStr([shape(mat, IK), shape(matA, IK)]))
#elif   OTP_ENABLED
        ndum = size(matA, 1, IK)
        CHECK_ASSERTION(__LINE__, size(mat, 2, IK) == size(matA, 2, IK), SK_"@setMatUpdate(): The condition `size(mat, 2) == size(matA, 2)` must hold. shape(mat), shape(matA) = "//getStr([shape(mat, IK), shape(matA, IK)]))
#else
#error  "Unrecognized interface."
#endif
#elif   EXP_ENABLED
        ! Ensure roffs and subset bounds are logical.
        CHECK_ASSERTION(__LINE__, all(0_IK <= [ndim, ndum]), SK_"@setMatUpdate(): The condition `all(0_IK <= [ndim, ndum])` must hold. ndim, ndum = "//getStr([ndim, ndum]))
        CHECK_ASSERTION(__LINE__, all(0_IK <= 1_IK - lbound(mat, kind = IK)), SK_"@setMatUpdate(): The condition `all(0 <= [roff, coff])` must hold. roff, coff = "//getStr([1_IK - lbound(mat)]))
        CHECK_ASSERTION(__LINE__, all(0_IK <= 1_IK - lbound(matA, kind = IK)), SK_"@setMatUpdate(): The condition `all(0 <= [roffA, coffA])` must hold. roffA, coffA = "//getStr([1_IK - lbound(matA)]))
#if     ONO_ENABLED
        CHECK_ASSERTION(__LINE__, all([ndim, ndum] <= shape(matA, IK)), SK_"@setMatUpdate(): The condition `all([roffA + ndim, coffA + ndum] <= shape(matA))` must hold. roffA, coffA, ndim, ndum, shape(matA) = "//getStr([1_IK - lbound(matA, kind = IK), ndim, ndum, shape(matA, IK)]))
#elif   OTP_ENABLED
        CHECK_ASSERTION(__LINE__, all([ndum, ndim] <= shape(matA, IK)), SK_"@setMatUpdate(): The condition `all([roffA + ndum, coffA + ndim] <= shape(matA))` must hold. roffA, coffA, ndim, ndum, shape(matA) = "//getStr([1_IK - lbound(matA, kind = IK), ndim, ndum, shape(matA, IK)]))
#else
#error  "Unrecognized interface."
#endif
#else
#error  "Unrecognized interface."
#endif
        ! Set the updating components.
#if     CHM_ENABLED && CK_ENABLED
#define GET_CONJG(X) conjg(X)
#define GET_RE(X) X%re
#define GET_REAL(X)real(X, TKG)
#elif   CHM_ENABLED || (CSM_ENABLED && (ONO_ENABLED || OTP_ENABLED))
#define GET_CONJG(X) X
#define GET_REAL(X)X
#define GET_RE(X) X
#else
#error  "Unrecognized interface."
#endif
        ! Set loop bounds.
#if     SLD_ENABLED
#define POFFSET + 1
#define GET_SLICE(i,j,k) j, k
#define SLD_SET(X,Y) X = Y
#define SUD_SET(X,Y)
#elif   SUD_ENABLED
#define POFFSET - 1
#define SLD_SET(X,Y)
#define SUD_SET(X,Y) X = Y
#define GET_SLICE(i,j,k) i, j
#else
#error  "Unrecognized interface."
#endif
#if     BLAS_ENABLED && DISPATCH_ENABLED
        call blasXXRK(BLAS_UPLO, BLAS_TRANS, ndim, ndum, GET_RE(ALPHA), matA(1,1), size(matA, 1, IK), GET_RE(BETA), mat(1,1), size(mat, 1, IK))
#else
        ! Quick return if possible.
        if (ndim == 0_IK .or. ((GET_RE(ALPHA) == GET_RE(ZERO) .or. ndum == 0_IK) .and. GET_RE(BETA) == GET_RE(ONE))) return
        if (GET_RE(ALPHA) == GET_RE(ZERO)) then
            if (GET_RE(BETA) == GET_RE(ZERO)) then
                do icol = 1, ndim
                    do irow = GET_SLICE(1, icol, ndim)
                        mat(irow, icol) = ZERO
                    end do
                end do
            else
                do icol = 1, ndim
                    SLD_SET(mat(icol, icol), GET_RE(BETA) * GET_RE(mat(icol, icol)))
                    do irow = GET_SLICE(1, icol POFFSET, ndim)
                        mat(irow, icol) = GET_RE(BETA) * mat(irow, icol)
                    end do
                    SUD_SET(mat(icol, icol), GET_RE(BETA) * GET_RE(mat(icol, icol)))
                end do
            end if
            return
        end if
        ! Perform Update.
#if     ONO_ENABLED
        ! Form  mat := GET_RE(ALPHA) * matA * matA ** T + GET_RE(BETA) * mat.
        ! Form  mat := GET_RE(ALPHA) * matA * matA ** H + GET_RE(BETA) * mat.
        do icol = 1, ndim
            if (GET_RE(BETA) == GET_RE(ZERO)) then
                do irow = GET_SLICE(1, icol, ndim)
                    mat(irow, icol) = ZERO
                end do
            else if (GET_RE(BETA) /= GET_RE(ONE)) then
                SLD_SET(mat(icol, icol), GET_RE(BETA) * GET_RE(mat(icol, icol)))
                do irow = GET_SLICE(1, icol POFFSET, ndim)
                    mat(irow, icol) = GET_RE(BETA) * mat(irow, icol)
                end do
                SUD_SET(mat(icol, icol), GET_RE(BETA) * GET_RE(mat(icol, icol)))
#if         CHM_ENABLED && CK_ENABLED
            else
                mat(icol, icol)%im = ZERO%im
#endif
            end if
            do idum = 1, ndum
                if (matA(icol, idum) /= ZERO) then
                    temp = GET_RE(ALPHA) * GET_CONJG(matA(icol, idum))
                    SLD_SET(mat(icol, icol), GET_RE(mat(icol, icol)) + GET_REAL(temp * matA(icol, idum)))
                    do irow = GET_SLICE(1, icol POFFSET, ndim)
                        mat(irow, icol) = mat(irow, icol) + temp * matA(irow, idum)
                    end do
                    SUD_SET(mat(icol, icol), GET_RE(mat(icol, icol)) + GET_REAL(temp * matA(icol, idum)))
                end if
            end do
        end do
#elif   OTP_ENABLED && CSM_ENABLED
        ! Form mat := GET_RE(ALPHA) * matA ** T * matA + GET_RE(BETA) * mat.
        ! This block is indeed the same as the block for HA_ENABLED.
        ! However, it is kept separately for now for its cleaner implementation.
        do icol = 1, ndim
            do irow = GET_SLICE(1, icol, ndim)
                temp = ZERO
                do idum = 1, ndum
                    temp = temp + matA(idum, irow) * matA(idum, icol)
                end do
                if (GET_RE(BETA) == ZERO) then
                    mat(irow, icol) = GET_RE(ALPHA) * temp
                else
                    mat(irow, icol) = GET_RE(ALPHA) * temp + GET_RE(BETA) * mat(irow, icol)
                end if
            end do
        end do
#elif   OTP_ENABLED && CHM_ENABLED
        ! Form mat := GET_RE(ALPHA) * matA ** H * matA + GET_RE(BETA) * mat.
#define SET_DIAGONAL \
GET_RE(temp) = GET_RE(ZERO); \
do idum = 1, ndum; \
    GET_RE(temp) = GET_RE(temp) + GET_REAL(GET_CONJG(matA(idum, icol)) * matA(idum, icol)); \
end do; \
if (GET_RE(BETA) == GET_RE(ZERO)) then; \
    mat(icol, icol) = GET_RE(ALPHA) * GET_RE(temp); \
else; \
    mat(icol, icol) = GET_RE(ALPHA) * GET_RE(temp) + GET_RE(BETA) * GET_RE(mat(icol, icol)); \
end if;
#define SET_OFFDIAG \
do irow = GET_SLICE(1, icol POFFSET, ndim); \
    temp = ZERO; \
    do idum = 1, ndum; \
        temp = temp + GET_CONJG(matA(idum, irow)) * matA(idum, icol); \
    end do; \
    if (GET_RE(BETA) == GET_RE(ZERO)) then; \
        mat(irow, icol) = GET_RE(ALPHA) * temp; \
    else; \
        mat(irow, icol) = GET_RE(ALPHA) * temp + GET_RE(BETA) * mat(irow, icol); \
    end if; \
end do;
        do icol = 1, ndim
#if         SLD_ENABLED
            SET_DIAGONAL
            SET_OFFDIAG
#elif       SUD_ENABLED
            SET_OFFDIAG
            SET_DIAGONAL
#endif
        end do
#else
#error  "Unrecognized interface."
#endif
#endif

        !%%%%%%%%%%%%%%%%%%%%%
#elif   setMatUpdateR1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, j, jy, ix, kx, effLenX, srow, erow
        TYPE_KIND :: temp
        ! Define transpose type for complex arguments.
#if     setMatUpdateR1H_CK_ENABLED || setMatUpdateR1AH_CK_ENABLED
#define GET_CONJG(X) conjg(X)
#else
#define GET_CONJG(X) X
#endif
        ! Define the `alpha` factor of the transformation.
#if     Fixed_ENABLED
#define ALPHA_TIMES
#elif   Alpha_ENABLED
#define ALPHA_TIMES alpha *
        if (alpha == ZERO) return ! Quick return if possible.
#else
#error  "Unrecognized interface."
#endif
        kx = 1_IK
        if (incA /= 1_IK) then
            CHECK_ASSERTION(__LINE__, incA /= 0_IK, SK_"@setMatUpdateR1(): The condition `incA /= 0_IK` must hold. incA = "//getStr(incA)) ! fpp
            effLenX = (size(vecA, 1, IK) - 1_IK) / abs(incA) + 1_IK ! the effective length of `vecA`.
            CHECK_ASSERTION(__LINE__, incA * effLenX /= size(vecA, 1, IK), \
            SK_"@setMatUpdateR1(): The condition `size(vecA) == incA * ((size(vecA, 1) - 1) / abs(incA) + 1)` must hold. incA, size(vecA) = "//\
            getStr([incA, size(vecA, 1, IK)])) ! fpp
            if (incA < 0_IK) kx = size(vecA, 1, IK)
        else
            effLenX = size(vecA, 1, IK)
        end if

        jy = 1_IK
        if (incB /= 1_IK) then
            CHECK_ASSERTION(__LINE__, incB /= 0_IK, SK_"@setMatUpdateR1(): The condition `incB /= 0_IK` must hold. incB = "//getStr(incB)) ! fpp
            if (incB < 0_IK) jy = size(vecB, 1, IK)
        end if

        CHECK_ASSERTION(__LINE__, 0_IK <= roff, SK_"@setMatUpdateR1(): The conditions `0_IK <= roff` must hold. roff = "//getStr(roff)) ! fpp

        CHECK_ASSERTION(__LINE__, size(mat, 2, IK) == (size(vecB, 1, IK) - 1_IK) / abs(incB) + 1_IK, \
        SK_"@setMatUpdateR1(): The condition `size(mat, 2) == (size(vecB) - 1) / abs(incB) + 1` must hold. size(mat, 2), (size(vecB) - 1) / abs(incB) + 1 = "//\
        getStr([size(mat, 2, IK), (size(vecB, 1, IK) - 1_IK) / abs(incB) + 1_IK])) ! fpp
        CHECK_ASSERTION(__LINE__, size(mat, 1, IK) >= effLenX + roff, \
        SK_"@setMatUpdateR1(): The conditions `size(mat, 1) >= (size(vecA, 1) - 1) / abs(incA) + 1 + roff` must hold. size(mat, 1), size(vecA), incA, roff = "//\
        getStr([size(mat, 1, IK), size(vecA, 1, IK), incA, roff])) ! fpp

        !   Start the operations. In this version the elements of `mat` are accessed sequentially with one pass through `mat`.

        srow = roff + 1_IK
        erow = roff + effLenX
        if (incA == 1_IK) then
            do j = 1_IK, size(mat, 2, IK)
                if (vecB(jy) /= ZERO) then
                    temp = ALPHA_TIMES GET_CONJG(vecB(jy))
                    do i = srow, erow
                        mat(i,j) = mat(i,j) + vecA(i) * temp
                    end do
                end if
                jy = jy + incB
            end do
        else
            do j = 1, size(mat, 2, IK)
                if (vecB(jy) /= ZERO) then
                    temp = ALPHA_TIMES GET_CONJG(vecB(jy))
                    ix = kx
                    do i = srow, erow
                        mat(i,j) = mat(i,j) + vecA(ix) * temp
                        ix = ix + incA
                    end do
                end if
                jy = jy + incB
            end do
        end if
#undef  ALPHA_TIMES

        !%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatUpdateTriang_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: irow, icol, idum
        TYPE_KIND :: temp

        ! Ensure roffs and subset bounds are logical.
        CHECK_ASSERTION(__LINE__, all(0_IK <= [ndim, ndum]), \
        SK_"@setMatUpdateTriang(): The condition `all(0_IK <= [ndim, ndum])` must hold. ndim, ndum = "//getStr([ndim, ndum])) ! fpp
        CHECK_ASSERTION(__LINE__, all(0_IK <= 1_IK - lbound(matA, kind = IK)), \
        SK_"@setMatUpdateTriang(): The condition `all(0 <= [roffA, coffA])` must hold. roffA, coffA = "//getStr([1_IK - lbound(matA)])) ! fpp
        CHECK_ASSERTION(__LINE__, all(0_IK <= 1_IK - lbound(mat, kind = IK)), \
        SK_"@setMatUpdateTriang(): The condition `all(0 <= [roff, coff])` must hold. roff, coff = "//getStr([1_IK - lbound(mat)])) ! fpp
        ! Ensure subset of `matA` does not overflow the matrix bounds.
#if     AS_ENABLED || AH_ENABLED
        CHECK_ASSERTION(__LINE__, all([ndim, ndum] <= shape(matA, kind = IK)), \
        SK_"@setMatUpdateTriang(): The condition `all([roffA + ndim, coffA + ndum] <= shape(matA))` must hold. roffA, coffA, ndim, ndum, shape(matA) = "\
        //getStr([1_IK - lbound(matA, kind = IK), ndim, ndum, shape(matA, IK)])) ! fpp
#elif   SA_ENABLED || HA_ENABLED
        CHECK_ASSERTION(__LINE__, all([ndum, ndim] <= shape(matA,IK)), \
        SK_"@setMatUpdateTriang(): The condition `all([roffA + ndum, coffA + ndim] <= shape(matA))` must hold. roffA, coffA, ndim, ndum, shape(matA) = "\
        //getStr([1_IK - lbound(matA, kind = IK), ndim, ndum, shape(matA, IK)])) ! fpp
#else
#error  "Unrecognized interface."
#endif
        ! Set the updating components.
#if     AS_ENABLED || SA_ENABLED
#define GET_CONJG(X) X
#define GET_REAL(X)X
#define GET_RE(X) X
#elif   AH_ENABLED || HA_ENABLED
#define GET_CONJG(X) conjg(X)
#define GET_RE(X) X%re
#define GET_REAL(X)real(X, TKG)
#else
#error  "Unrecognized interface."
#endif
        ! Set loop bounds.
#if     lowDiaC_ENABLED
#define POFFSET + 1
#define GET_SLICE(i,j,k) j, k
#define lowDiaC_SET(X,Y) X = Y
#define uppDiaC_SET(X,Y)
#elif   uppDiaC_ENABLED
#define POFFSET - 1
#define lowDiaC_SET(X,Y)
#define uppDiaC_SET(X,Y) X = Y
#define GET_SLICE(i,j,k) i, j
#else
#error  "Unrecognized interface."
#endif
        ! Quick return if possible.
        if (ndim == 0_IK .or. ((alpha == ZERO .or. ndum == 0_IK) .and. beta == ONE)) return
        if (alpha == GET_RE(ZERO)) then
            if (beta == GET_RE(ZERO)) then
                do icol = 1, ndim
                    do irow = GET_SLICE(1, icol, ndim)
                        mat(irow, icol) = ZERO
                    end do
                end do
            else
                do icol = 1, ndim
                    lowDiaC_SET(mat(icol, icol), beta * GET_RE(mat(icol, icol)))
                    do irow = GET_SLICE(1, icol POFFSET, ndim)
                        mat(irow, icol) = beta * mat(irow, icol)
                    end do
                    uppDiaC_SET(mat(icol, icol), beta * GET_RE(mat(icol, icol)))
                end do
            end if
            return
        end if
        ! Perform Update.
#if     AS_ENABLED || AH_ENABLED
        ! Form  mat := alpha * matA * matA ** T + beta * mat.
        ! Form  mat := alpha * matA * matA ** H + beta * mat.
        do icol = 1, ndim
            if (beta == GET_RE(ZERO)) then
                do irow = GET_SLICE(1, icol, ndim)
                    mat(irow, icol) = ZERO
                end do
            else if (beta /= GET_RE(ONE)) then
                lowDiaC_SET(mat(icol, icol), beta * GET_RE(mat(icol, icol)))
                do irow = GET_SLICE(1, icol POFFSET, ndim)
                    mat(irow, icol) = beta * mat(irow, icol)
                end do
                uppDiaC_SET(mat(icol, icol), beta * GET_RE(mat(icol, icol)))
#if         AH_ENABLED
            else
                mat(icol, icol)%im = ZERO%im
#endif
            end if
            do idum = 1, ndum
                if (matA(icol, idum) /= ZERO) then
                    temp = alpha * GET_CONJG(matA(icol, idum))
                    lowDiaC_SET(mat(icol, icol), GET_RE(mat(icol, icol)) + GET_REAL(temp * matA(icol, idum)))
                    do irow = GET_SLICE(1, icol POFFSET, ndim)
                        mat(irow, icol) = mat(irow, icol) + temp * matA(irow, idum)
                    end do
                    uppDiaC_SET(mat(icol, icol), GET_RE(mat(icol, icol)) + GET_REAL(temp * matA(icol, idum)))
                end if
            end do
        end do
#elif   SA_ENABLED
        ! Form mat := alpha * matA ** T * matA + beta * mat.
        ! This block is indeed the same as the block for HA_ENABLED.
        ! However, it is kept separately for now for its cleaner implementation.
        do icol = 1, ndim
            do irow = GET_SLICE(1, icol, ndim)
                temp = ZERO
                do idum = 1, ndum
                    temp = temp + matA(idum, irow) * matA(idum, icol)
                end do
                if (beta == ZERO) then
                    mat(irow, icol) = alpha * temp
                else
                    mat(irow, icol) = alpha * temp + beta * mat(irow, icol)
                end if
            end do
        end do
#elif   HA_ENABLED
        ! Form mat := alpha * matA ** H * matA + beta * mat.
#define SET_DIAGONAL \
GET_RE(temp) = GET_RE(ZERO); \
do idum = 1, ndum; \
    GET_RE(temp) = GET_RE(temp) + GET_REAL(GET_CONJG(matA(idum, icol)) * matA(idum, icol)); \
end do; \
if (beta == GET_RE(ZERO)) then; \
    mat(icol, icol) = alpha * GET_RE(temp); \
else; \
    mat(icol, icol) = alpha * GET_RE(temp) + beta * GET_RE(mat(icol, icol)); \
end if;
#define SET_OFFDIAG \
do irow = GET_SLICE(1, icol POFFSET, ndim); \
    temp = ZERO; \
    do idum = 1, ndum; \
        temp = temp + GET_CONJG(matA(idum, irow)) * matA(idum, icol); \
    end do; \
    if (beta == GET_RE(ZERO)) then; \
        mat(irow, icol) = alpha * temp; \
    else; \
        mat(irow, icol) = alpha * temp + beta * mat(irow, icol); \
    end if; \
end do;
        do icol = 1, ndim
#if         lowDiaC_ENABLED
            SET_DIAGONAL
            SET_OFFDIAG
#elif       uppDiaC_ENABLED
            SET_OFFDIAG
            SET_DIAGONAL
#endif
        end do
#else
#error  "Unrecognized interface."
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  SET_DIAGONAL
#undef  SET_OFFDIAG
#undef  lowDiaC_SET
#undef  uppDiaC_SET
#undef  BLAS_TRANS
#undef  BLAS_UPLO
#undef  TYPE_KIND
#undef  GET_CONJG
#undef  GET_SLICE
#undef  blasXXRK
#undef  GET_REAL
#undef  SLD_SET
#undef  SUD_SET
#undef  POFFSET
#undef  GET_RE
#undef  ALPHA
#undef  BETA