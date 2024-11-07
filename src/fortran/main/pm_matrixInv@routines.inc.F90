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
!>  This file contains procedure implementations of [pm_matrixInv](@ref pm_matrixInv).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Set the type and kind.
#if     CK_ENABLED
#define GET_CONJG(X) conjg(X)
        complex(TKG), parameter :: ZERO = (0._TKG, 0._TKG), ONE = (1._TKG, 0._TKG)
#define TYPE_KIND complex(TKG)
#elif   RK_ENABLED
#define GET_CONJG(X) X
        real(TKG), parameter :: ZERO = 0._TKG, ONE = 1._TKG
#define TYPE_KIND real(TKG)
#else
#error  "Unrecognized interface."
#endif
        ! Define the runtime checks.
#define CHECK_SHAPE_MAT \
CHECK_ASSERTION(__LINE__, size(mat, 1, IK) == size(mat, 2, IK), SK_"@setMatInv(): The condition `size(mat, 1) == size(mat, 2)` must hold. shape(mat) = "//getStr(shape(mat, IK)))
#define CHECK_SHAPE_INV \
CHECK_ASSERTION(__LINE__, all(shape(inv, IK) == shape(mat, IK)), SK_"@setMatInv(): The condition `all(shape(inv) == shape(mat))` must hold. shape(inv), shape(mat) = "//getStr([shape(inv, IK), shape(mat, IK)]))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getMatInv_ENABLED && (Def_ENABLED || Det_ENABLED || Inf_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: rperm(size(mat, 1, IK))
        TYPE_KIND :: lup(size(mat, 1, IK), size(mat, 2, IK))
#if     Def_ENABLED || Det_ENABLED
#define AUXIL info
        integer(IK) :: info
#elif   !Inf_ENABLED
#error  "Unrecognized interface."
#endif
        lup = mat
        call setMatLUP(lup, rperm, AUXIL)
        if (AUXIL /= 0_IK) then
#if         Inf_ENABLED
            return
#else
            error stop MODULE_NAME//SK_"@getMatInv(): LUP factorization of the input `mat` failed."
#endif
        end if
        call setMatInv(inv, lup, rperm)
#if     Det_ENABLED
        block
            integer(IK) :: idim
            auxil = ONE
            do idim = 1, size(rperm, 1, IK)
                auxil = auxil * lup(idim, idim)
                if (rperm(idim) /= idim) auxil = -auxil
            end do
        end block
#endif
#undef  AUXIL

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getMatInv_ENABLED && (CUD_ENABLED || CUU_ENABLED || CLD_ENABLED || CLU_ENABLED || LUP_ENABLED || CCU_ENABLED || CCL_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setMatInv(inv, mat, auxil)
!#if     CCU_ENABLED
!        call setMatCopy(inv(1:,2:), rdpack, inv(2:,1:), rdpack, lowDia, transHerm)
!#elif   CCL_ENABLED
!        call setMatCopy(inv(2:,1:), rdpack, inv(1:,2:), rdpack, uppDia, transHerm)
!#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatInv_ENABLED && (CUD_ENABLED || CLD_ENABLED || CUU_ENABLED || CLU_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Set the indexing rules based on the subset type.
#if     CUD_ENABLED || CUU_ENABLED
#define GET_INDEX(I,J) J,I
#elif   CLD_ENABLED || CLU_ENABLED
#define GET_INDEX(I,J) I,J
#else
#error  "Unrecognized interface."
#endif
        TYPE_KIND   :: summ
        integer(IK) :: irow, icol, idim, ndim
        ndim = size(mat, 1, IK)
        CHECK_SHAPE_MAT
        CHECK_SHAPE_INV
        do irow = 1_IK, ndim
#if         CLD_ENABLED || CLU_ENABLED
            inv(1 : irow - 1, irow) = ZERO
#endif
#if         CUD_ENABLED || CLD_ENABLED
            CHECK_ASSERTION(__LINE__, mat(irow, irow) /= ZERO, SK_"@setMatInv(): The condition `mat(irow, irow) /= 0` must hold. irow = "//getStr(irow))
            inv(irow, irow) = ONE / mat(irow, irow)
#elif       CUU_ENABLED || CLU_ENABLED
            CHECK_ASSERTION(__LINE__, mat(irow, irow) == ONE, SK_"@setMatInv(): The condition `mat(irow, irow) == 1` must hold. irow = "//getStr(irow))
            inv(irow, irow) = ONE
#endif
#if         CUD_ENABLED || CUU_ENABLED
            inv(irow + 1 : ndim, irow) = ZERO
#endif
            do icol = 1_IK, irow - 1_IK
                summ = ZERO
                do idim = icol, irow - 1_IK
                    summ = summ + mat(GET_INDEX(irow, idim)) * inv(GET_INDEX(idim, icol))
                end do
#if             CUD_ENABLED || CLD_ENABLED
                inv(GET_INDEX(irow, icol)) = -summ * inv(irow, irow)
#elif           CUU_ENABLED || CLU_ENABLED
                inv(GET_INDEX(irow, icol)) = -summ
#endif
            end do
        end do
#undef  GET_INDEX

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatInv_ENABLED && (CCU_ENABLED || CCL_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Set the indexing rules based on the subset type.
#if     CCL_ENABLED && (FUL_ENABLED || UXD_ENABLED)
#define INDMAT(I,J)I,J
#define INDINV(I,J)I,J
#elif   CCU_ENABLED && (FUL_ENABLED || XLD_ENABLED)
#define INDMAT(I,J)J,I
#define INDINV(I,J)J,I
#elif   CCL_ENABLED && XLD_ENABLED
#define INDMAT(I,J)I,J
#define INDINV(I,J)J,I
#elif   CCU_ENABLED && UXD_ENABLED
#define INDMAT(I,J)J,I
#define INDINV(I,J)I,J
#else
#error  "Unrecognized interface."
#endif
       !TYPE_KIND :: summ
        integer(IK) :: irow, icol, ndim !, idim
        ndim = size(mat, 1, IK)
        CHECK_SHAPE_MAT
        CHECK_SHAPE_INV
        if (1_IK < ndim) then
            CHECK_ASSERTION(__LINE__, all(real(getMatCopy(lfpack, mat, rdpack, dia), TKG) /= 0._TKG), \
            SK_"@setMatInv(): The The diagonal elements of the Cholesky factorization must be non-zero. getMatCopy(lfpack, mat, rdpack, dia) = "//getStr(getMatCopy(lfpack, mat, rdpack, dia)))
            do irow = 1_IK, ndim - 1_IK
                inv(irow, irow) = ONE / mat(irow, irow)
                do icol = irow + 1_IK, ndim
#if                 FUL_ENABLED || (CCL_ENABLED && UXD_ENABLED) || (CCU_ENABLED && XLD_ENABLED)
                    inv(INDINV(irow, icol)) = -dot_product(mat(INDMAT(icol, irow : icol - 1)), inv(INDINV(irow, irow : icol - 1))) / mat(icol, icol)
#elif               (CCU_ENABLED && UXD_ENABLED) || (CCL_ENABLED && XLD_ENABLED)
                    inv(INDINV(irow, icol)) = -GET_CONJG(dot_product(mat(INDMAT(icol, irow : icol - 1)), GET_CONJG(inv(INDINV(irow, irow : icol - 1))))) / mat(icol, icol)
#else
#error              "Unrecognized interface."
#endif
                end do
            end do
            inv(irow, irow) = ONE / mat(irow, irow)
            do irow = 1_IK, ndim
                do icol = irow, ndim
                    inv(INDINV(irow, icol)) = dot_product(inv(INDINV(icol, icol : ndim)), inv(INDINV(irow, icol : ndim)))
#if                 FUL_ENABLED
                    inv(INDINV(icol, irow)) = GET_CONJG(inv(INDINV(irow, icol)))
#endif
                end do
            end do
        elseif (ndim == 1_IK) then
            inv(1,1) = ONE / mat(1,1)**2
        end if
#undef  INDINV
#undef  INDMAT

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatInv_ENABLED && LUP_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_KIND   :: summ
        integer(IK) :: irow, icol, idim, ndim
        ndim = size(mat, 1, IK)
        CHECK_SHAPE_MAT
        CHECK_SHAPE_INV
        CHECK_ASSERTION(__LINE__, size(auxil, 1, IK) == ndim, SK_"@setMatInv(): The condition `size(auxil) == size(inv, 1)` must hold. size(auxil), shape(inv) = "//getStr([size(auxil, 1, IK), shape(inv, IK)]))
        do icol = 1, ndim
            inv(1 : icol - 1, icol) = ZERO
            inv(icol, icol) = ONE
            inv(icol + 1 : ndim, icol) = ZERO
            idim = 0_IK
            do irow = 1_IK, ndim
                summ = inv(auxil(irow), icol)
                inv(auxil(irow), icol) = inv(irow, icol)
                if (idim /= 0_IK) then
                    summ = summ - dot_product(GET_CONJG(mat(irow, idim : irow - 1)), inv(idim : irow - 1, icol))
                elseif (summ /= ZERO) then
                    idim = irow
                end if
                inv(irow, icol) = summ
            end do
            do irow = ndim, 1_IK, -1_IK
                CHECK_ASSERTION(__LINE__, mat(irow, irow) /= ZERO, SK_"@setMatInv(): The condition `mat(irow, irow) /= ZERO` must hold. irow = "//getStr(irow))
                inv(irow, icol) = (inv(irow, icol) - dot_product(GET_CONJG(mat(irow, irow + 1 : ndim)), inv(irow + 1 : ndim, icol))) / mat(irow, irow)
            end do
        end do

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  CHECK_SHAPE_INV
#undef  CHECK_SHAPE_MAT
#undef  TYPE_KIND
#undef  GET_CONJG