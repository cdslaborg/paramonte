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
!>  This include file contains procedure implementation of the generic interface [pm_matrixDet](@ref pm_matrixDet).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK_ENABLED
#define TYPE_KIND complex(TKG)
#define GET_CONJG(X) conjg(X)
#elif   RK_ENABLED
#define GET_CONJG(X) X
#define TYPE_KIND real(TKG)
#else
#error  "Unrecognized interface."
#endif
        ! Define the output variable.
#if     getMatDetSqrt_ENABLED || setMatDetSqrt_ENABLED
#define GETLOG(X) X
#define OUTPUT detSqrt
#define INCREMENT(X,Y) X = X * Y
#elif   getMatDetSqrtLog_ENABLED || setMatDetSqrtLog_ENABLED
#define GETLOG(X) log(X)
#define OUTPUT detSqrtLog
#define INCREMENT(X,Y) X = X + log(Y)
#elif   !(getMatDet_ENABLED || setMatDet_ENABLED)
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getMatDet_ENABLED || setMatDet_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getMatDet_ENABLED
        integer(IK) :: info
        TYPE_KIND :: matLUP(size(mat, 1, IK), size(mat, 2, IK))
#elif   setMatDet_ENABLED
#define matLUP mat
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: idim, ndim
        integer(IK) :: rperm(size(mat,1))
        ndim = size(mat, 1, IK)
        CHECK_ASSERTION(__LINE__, size(mat, 2, IK) == ndim, SK_"The input `mat` must be square. shape(mat) = "//getStr(shape(mat)))
#if     getMatDet_ENABLED
        matLUP = mat
#endif
        call setMatLUP(matLUP, rperm, info) ! This returns det as +-1.
        if (info /= 0_IK) then
#if         getMatDet_ENABLED
            error stop "LU factorization failed." ! LCOV_EXCL_LINE
#elif       setMatDet_ENABLED
            return
#endif
        end if
        det = 1 ! parity
        do idim = 1_IK, ndim
            det = det * matLUP(idim, idim)
            if (rperm(idim) /= idim) det = -det
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getMatDetSqrt_ENABLED || getMatDetSqrtLog_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getMatDetSqrt_ENABLED
#define SETDET(SUBSET) \
call setMatDetSqrt(mat, SUBSET, detSqrt, info, chol, nothing)
#elif   getMatDetSqrtLog_ENABLED
#define SETDET(SUBSET) \
call setMatDetSqrtLog(mat, SUBSET, detSqrtLog, info, chol, nothing)
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: info
        TYPE_KIND :: chol(size(mat, 1, IK), size(mat, 2, IK))
        if (size(mat, 1, IK) < 64_IK) then
            ! Use the unblocked algorithm from this module.
            if (present(subset)) then
                select type (subset)
                type is (uppDia_type)
                    SETDET(uppDia)
                type is (lowDia_type)
                    SETDET(lowDia)
                class default
                    error stop "Unrecognized `subset`. Only objects of `uppDia_type()` and `lowDia_type()` are recognized." ! LCOV_EXCL_LINE
                end select
            else
                SETDET(uppDia)
            end if
            if (info /= 0_IK) error stop "The Cholesky factorization of small matrix failed."
        else
            ! Use the blocked algorithm potentially dispatched to LAPACK if available.
            if (present(subset)) then
                select type (subset)
                type is (uppDia_type)
                    call setMatCopy(chol, rdpack, mat, rdpack, uppDia)
                    call setMatChol(chol, uppDia, info, iteration)
                type is (lowDia_type)
                    call setMatCopy(chol, rdpack, mat, rdpack, lowDia)
                    call setMatChol(chol, lowDia, info, iteration)
                class default
                    error stop "Unrecognized `subset`. Only objects of `uppDia_type()` and `lowDia_type()` are recognized." ! LCOV_EXCL_LINE
                end select
            else
                call setMatCopy(chol, rdpack, mat, rdpack, uppDia)
                call setMatChol(chol, uppDia, info, iteration)
            end if
            if (info /= 0_IK) error stop "The Cholesky factorization of large matrix failed."
            OUTPUT = GETLOG(real(chol(1, 1), TKG))
            do info = 2, size(chol, 1, IK)
                INCREMENT(OUTPUT,real(chol(info, info), TKG))
            end do
            info = 0_IK
        end if
#undef  SETDET

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatDetSqrt_ENABLED || setMatDetSqrtLog_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKG) :: summ
        integer(IK) :: irow
!#if     IMP_ENABLED
        integer(IK) :: ndim
        ndim = size(mat, 1, IK)
        CHECK_ASSERTION(__LINE__, size(mat, 1, IK) == size(mat, 2, IK), SK_"@setMatChol(): The condition `size(mat, 1) == size(mat, 2)` must hold. ndim = "//getStr([size(mat, 1, IK), size(mat, 2, IK)]))
        CHECK_ASSERTION(__LINE__, all(shape(mat, IK) == shape(chol, IK)), SK_"@setMatChol(): The condition `all(shape(mat) == shape(chol))` must hold. shape(mat), shape(chol) = "//getStr([shape(mat, IK), shape(chol, IK)]))
!#elif   EXP_ENABLED
!        ! Ensure offsets and subset bounds are logical.
!        check_assertion(__LINE__, 0_IK <= ndim, SK_"@setMatChol(): The condition `0 <= ndim` must hold. ndim = "//getStr(ndim))
!        check_assertion(__LINE__, all(0_IK <= 1_IK - lbound(mat, kind = IK)), SK_"@setMatChol(): The condition `all(0 <= [roff, coff])` must hold. roff, coff = "//getStr(1_IK - lbound(mat)))
!        check_assertion(__LINE__, all(ndim <= ubound(mat, kind = IK)), SK_"@setMatChol(): The condition `ndim + [roff, coff] <= shape(mat)` must hold. ndim, roff, coff, shape(mat) = "//getStr([ndim, 1 - lbound(mat, kind = IK), shape(mat, kind = IK)]))
!        check_assertion(__LINE__, all(ndim <= ubound(chol, kind = IK)), SK_"@setMatChol(): The condition `ndim + [roffc, coffc] <= shape(chol)` must hold. ndim, roffc, coffc, shape(chol) = "//getStr([ndim, 1 - lbound(chol, kind = IK), shape(chol, kind = IK)]))
!#else
!#error  "Unrecognized interface."
!#endif
        CHECK_ASSERTION(__LINE__, all(ndim <= ubound(chol, kind = IK)), SK_"@setMatChol(): The condition `all(shape(mat) == shape(chol))` must hold. shape(mat), ubound(chol) = "//getStr([shape(mat, IK), ubound(chol, kind = IK)]))
        ! Define operations based on the specified triangular side.
#if     UXD_ENABLED
#define GET_MATMUL(I,J) matmul(J,I)
#define GETMAT(MAT,I,J) MAT(J,I)
#elif   XLD_ENABLED
#define GET_MATMUL(I,J) matmul(I,J)
#define GETMAT(MAT,I,J) MAT(I,J)
#else
#error  "Unrecognized interface."
#endif
        if (1_IK < ndim) then
            info = 1_IK
            summ = real(mat(info, info), TKG)
            if (0._TKG < summ) then
                summ = sqrt(summ)
                chol(info, info) = summ
                OUTPUT = GETLOG(summ)
                summ = 1._TKG / summ
#if             ONO_ENABLED
                do irow = info + 1, ndim
                    GETMAT(chol, irow, info) = (GETMAT(mat, irow, info) - dot_product(GETMAT(chol, info, 1 : info - 1), GETMAT(chol, irow, 1 : info - 1))) * summ
                end do
#elif           OTH_ENABLED
                do irow = info + 1, ndim
                    GETMAT(chol, info, irow) = GET_CONJG(GETMAT(mat, irow, info) - dot_product(GETMAT(chol, 1 : info - 1, irow), GETMAT(chol, 1 : info - 1, info))) * summ
                end do
#endif
                do info = 2_IK, ndim
#if                 ONO_ENABLED
                    summ = real(mat(info, info), TKG) - real(dot_product(GETMAT(chol, info, 1 : info - 1), GETMAT(chol, info, 1 : info - 1)), TKG)
#elif               OTH_ENABLED
                    summ = real(mat(info, info), TKG) - real(dot_product(GETMAT(chol, 1 : info - 1, info), GETMAT(chol, 1 : info - 1, info)), TKG)
#else
#error              "Unrecognized interface."
#endif
                    if (0._TKG < summ) then
                        summ = sqrt(summ)
                        chol(info, info) = summ
                        INCREMENT(OUTPUT,summ)
                        summ = 1._TKG / summ
#if                     ONO_ENABLED
                        !GETMAT(chol, info + 1 : ndim, info) = (GETMAT(mat, info + 1 : ndim, info) - GET_MATMUL(GET_CONJG(GETMAT(chol, info + 1 : ndim, 1 : info - 1)),GETMAT(chol, info, 1 : info - 1))) * summ
                        do irow = info + 1, ndim
                            GETMAT(chol, irow, info) = (GETMAT(mat, irow, info) - dot_product(GETMAT(chol, info, 1 : info - 1), GETMAT(chol, irow, 1 : info - 1))) * summ
                        end do
#elif                   OTH_ENABLED
                        do irow = info + 1, ndim
                            GETMAT(chol, info, irow) = (GET_CONJG(GETMAT(mat, irow, info) - dot_product(GETMAT(chol, 1 : info - 1, irow), GETMAT(chol, 1 : info - 1, info)))) * summ
                        end do
#endif
                        cycle
                    end if
                    return
                end do
                info = 0_IK
            end if
        elseif (1_IK == ndim) then
            summ = real(mat(1, 1), TKG)
            if (0._TKG < real(summ, TKG)) then
                summ = sqrt(summ)
                OUTPUT = GETLOG(summ)
                chol(1,1) = summ
                info = 0_IK
            else
                info = 1_IK
            end if
        else
            info = 0_IK
        end if
#undef  GET_MATMUL
#undef  GETMAT

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  INCREMENT
#undef  TYPE_KIND
#undef  GET_CONJG
#undef  OUTPUT
#undef  GETLOG
#undef  matLUP