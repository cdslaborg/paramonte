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
!>  This file contains the implementation details of the 2D routines under the generic interface [pm_sampleCor](@ref pm_sampleCor).
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Saturday 4:40 PM, August 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Set the conjugation rule.
#if     CK_ENABLED
#define TYPE_OF_SAMPLE complex(TKC)
#define GET_CONJG(X)conjg(X)
#define GET_RE(X)X%re
#elif   RK_ENABLED
#define TYPE_OF_SAMPLE real(TKC)
#define GET_CONJG(X)X
#define GET_RE(X)X
#elif   getCor_ENABLED || setCor_ENABLED
!elif   !(SK_ENABLED || IK_ENABLED || LK_ENABLED || PSSK_ENABLED || BSSK_ENABLED)
#error  "Unrecognized interface."
#endif
        ! Define weight type and kind and ZEROW.
#if     setRho_ENABLED && WTI_ENABLED
#define TYPE_OF_WEIGHT integer(IK)
#elif   setRho_ENABLED && (WTR_ENABLED || WNO_ENABLED)
#define TYPE_OF_WEIGHT real(RK)
#endif
        ! Define the weight arguments.
#if     WNO_ENABLED
#define WEIGHT_ARGS
#elif   WTI_ENABLED || WTR_ENABLED
#define WEIGHT_ARGS, weight, weisum
#elif   (Prs_ENABLED || getRho_ENABLED || setRho_ENABLED)
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%
#if     setCordance_ENABLED
        !%%%%%%%%%%%%%%%%%%

        ! Define the indexing rules.
#if     D0_ENABLED && SK_ENABLED
#define GET_SIZE(ARRAY)len(ARRAY, IK)
#define GETELL(ARRAY,I)ARRAY(I:I)
#elif   D1_ENABLED && (PSSK_ENABLED || BSSK_ENABLED)
#define GET_SIZE(ARRAY)size(ARRAY, 1, IK)
#define GETELL(ARRAY,I)ARRAY(I)%val
#elif   D1_ENABLED
#define GET_SIZE(ARRAY)size(ARRAY, 1, IK)
#define GETELL(ARRAY,I)ARRAY(I)
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: isam, jsam
        CHECK_ASSERTION(__LINE__, GET_SIZE(x) == GET_SIZE(y), SK_"@setCor(): The condition `size(x) == size(y)` must hold. size(x), size(y) = "//getStr([GET_SIZE(x), GET_SIZE(y)]))
        tiedx = 0_IK
        tiedy = 0_IK
        ! Set the output arguments.
#if     Sum_ENABLED
        cordance = 0_IK
#define INCREMENT_CORDANCE cordance = cordance + 1_IK
#define DECREMENT_CORDANCE cordance = cordance - 1_IK
#elif   All_ENABLED
        concordance = 0_IK
        discordance = 0_IK
#define INCREMENT_CORDANCE concordance = concordance + 1_IK
#define DECREMENT_CORDANCE discordance = discordance + 1_IK
#else
#error  "Unrecognized interface."
#endif
        do isam = 2, GET_SIZE(x)
            do jsam = 1, isam - 1
                if (GETELL(x,isam) < GETELL(x,jsam)) then
                    if (GETELL(y,isam) < GETELL(y,jsam)) then
                        INCREMENT_CORDANCE
                    elseif (GETELL(y,jsam) < GETELL(y,isam)) then
                        DECREMENT_CORDANCE
                    else
                        tiedy = tiedy + 1_IK
                    end if
                elseif (GETELL(x,jsam) < GETELL(x,isam)) then
                    if (GETELL(y,isam) < GETELL(y,jsam)) then
                        INCREMENT_CORDANCE
                    elseif (GETELL(y,jsam) < GETELL(y,isam)) then
                        DECREMENT_CORDANCE
                    else
                        tiedy = tiedy + 1_IK
                    end if
                else
                    tiedx = tiedx + 1_IK
                    if (GETELL(y,isam) == GETELL(y,jsam)) tiedy = tiedy + 1_IK
                end if
            end do
        end do
#undef  INCREMENT_CORDANCE
#undef  DECREMENT_CORDANCE
#undef  GET_SIZE
#undef  GETELL

!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!#elif   getCorMerged_ENABLED && New_ENABLED
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!        integer(IK) :: idim
!        call setCorMerged(corMerged, corB, corA, meanDiff, fracA, uppDia)
!        ! The `do concurrent` version is problematic for OpenMP builds of the ParaMonte library with the Intel ifort 2021.8 on heap, yielding segmentation fault.
!        !do concurrent(idim = 2 : size(corMerged, 1, IK))
!        do idim = 2, size(corMerged, 1, IK)
!            corMerged(idim, 1 : idim - 1) = corMerged(1 : idim - 1, idim)
!        end do
!
!        !%%%%%%%%%%%%%%%%%%%
!#elif   setCorMerged_ENABLED
!        !%%%%%%%%%%%%%%%%%%%
!
!        integer(IK) :: idim, jdim, ndim
!        real(TKC) :: fracB, fracAB, stdinv(size(meanDiff, 1, IK))
!        ! Define the loop ranges for different subsets.
!#if     UXD_ENABLED || UXX_ENABLED
!#define GET_RANGE(I,J,K)I, J - I
!#elif   XLD_ENABLED || XLX_ENABLED
!#define GET_RANGE(I,J,K)J + I, K
!#else
!#error  "Unrecognized interface."
!#endif
!        ! Define the output value for setCorMerged.
!#if     Old_ENABLED
!#define corMerged corB
!#elif   New_ENABLED
!        CHECK_ASSERTION(__LINE__, all(shape(corMerged, IK) == shape(corA, IK)), SK_"@setCorMerged(): The condition `all(shape(corMerged) == shape(corA))` must hold. shape(corMerged), shape(corA) = "//getStr([shape(corMerged, IK), shape(corA, IK)]))
!#else
!#error  "Unrecognized interface."
!#endif
!        CHECK_ASSERTION(__LINE__, 0._TKC < fracA .and. fracA < 1._TKC, SK_"@setCorMerged(): The condition `0 < fracA .and. fracA < 1` must hold. fracA = "//getStr(fracA))
!        CHECK_ASSERTION(__LINE__, all(shape(corB, IK) == shape(corA, IK)), SK_"@setCorMerged(): The condition `all(shape(corB, IK) == shape(corA, IK))` must hold. shape(corB), shape(corA) = "//getStr([shape(corB, IK), shape(corA, IK)]))
!        CHECK_ASSERTION(__LINE__, all(size(meanDiff, 1, IK) == shape(corA, IK)), SK_"@setCorMerged(): The condition `all(size(meanDiff) == shape(corA))` must hold. size(meanDiff), shape(corA) = "//getStr([size(meanDiff, 1, IK), shape(corA, IK)]))
!        fracB = 1._TKC - fracA
!        fracAB = fracA * fracB
!        ndim = size(corMerged, 1, IK)
!        do jdim = 1, size(corA, 1, IK)
!            stdinv(jdim) = 1._TKC / sqrt(1._TKC + fracAB * meanDiff(jdim)**2)
!#if         XLD_ENABLED
!            corMerged(jdim, jdim) = 1._TKC
!#endif
!            do idim = GET_RANGE(1,jdim,ndim)
!                corMerged(idim, jdim) = fracB * corB(idim, jdim) + fracA * corA(idim, jdim) + fracAB * meanDiff(idim) * meanDiff(jdim)
!            end do
!#if         UXD_ENABLED
!            corMerged(jdim, jdim) = 1._TKC
!#endif
!        end do
!        do jdim = 1, ndim
!            do idim = GET_RANGE(1,jdim,ndim)
!                corMerged(idim, jdim) = corMerged(idim, jdim) * stdinv(idim) * stdinv(jdim)
!            end do
!        end do
!#undef  GET_RANGE
!#undef  corMerged
!
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getCor_ENABLED && CFC_ENABLED && RULD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: idim
#if     VUXX_ENABLED || VXLX_ENABLED
        call setCor(cor, uppDia, cov, subsetv, stdinv)
#elif   VUXD_ENABLED || VXLD_ENABLED
        call setCor(cor, uppDia, cov, subsetv)
#else
#error  "Unrecognized interface."
#endif
        ! The `do concurrent` version is problematic for OpenMP builds of the ParaMonte library with the Intel ifort 2021.8 on heap, yielding segmentation fault.
        !do concurrent(idim = 2 : size(cor, 1, IK))
        do idim = 2, size(cor, 1, IK)
            cor(idim, 1 : idim - 1) = GET_CONJG(cor(1 : idim - 1, idim))
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getCor_ENABLED && (XY_ENABLED || ULD_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     WTI_ENABLED
        integer(IK) :: weisum
#elif   WTR_ENABLED
        real(TKC) :: weisum
#elif   !WNO_ENABLED
#error  "Unrecognized interface."
#endif
#if     XY_ENABLED
        TYPE_OF_SAMPLE :: mean(2)
        call setMean(mean, x, y WEIGHT_ARGS)
        call setCor(cor, mean, x, y WEIGHT_ARGS)
#elif   ULD_ENABLED
        TYPE_OF_SAMPLE :: mean(size(sample, 3 - dim, IK))
        call setMean(mean, sample, dim WEIGHT_ARGS)
        call setCor(cor, uppDia, mean, sample, dim WEIGHT_ARGS)
#else
#error  "Unrecognized interface."
#endif
#if     ULD_ENABLED
        block
            integer(IK) :: ndim, idim
            ndim = size(cor, 1, IK)
            ! The `do concurrent` version is problematic for OpenMP builds of the ParaMonte library with the Intel ifort 2021.8 on heap, yielding segmentation fault.
            !do concurrent(idim = 1 : ndim)
            do idim = 1, ndim
                cor(idim, 1 : idim - 1) = GET_CONJG(cor(1 : idim - 1, idim))
            end do
        end block
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setCor_ENABLED && CFC_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: idim, jdim, ndim
        ! Define the indexing rules for `cor`.
#if     RUXX_ENABLED || RUXD_ENABLED
#define RINDEX(I,J) I,J
#elif   RXLX_ENABLED || RXLD_ENABLED
#define RINDEX(I,J) J,I
#else
#error  "Unrecognized interface."
#endif
        ! Define the indexing rules for `cov`.
#if     VUXX_ENABLED || VUXD_ENABLED
#define VINDEX(I,J) I,J
#elif   VXLX_ENABLED || VXLD_ENABLED
#define VINDEX(I,J) J,I
#else
#error  "Unrecognized interface."
#endif
        ! Define the conjugation rule.
#if     (VUXD_ENABLED && (RUXD_ENABLED || RUXX_ENABLED)) || (VUXX_ENABLED && (RUXD_ENABLED || RUXX_ENABLED)) || \
        (VXLD_ENABLED && (RXLD_ENABLED || RXLX_ENABLED)) || (VXLX_ENABLED && (RXLD_ENABLED || RXLX_ENABLED))
#define CONJUGATE(X)X
#elif   (VUXD_ENABLED && (RXLD_ENABLED || RXLX_ENABLED)) || (VUXX_ENABLED && (RXLD_ENABLED || RXLX_ENABLED)) || \
        (VXLD_ENABLED && (RUXD_ENABLED || RUXX_ENABLED)) || (VXLX_ENABLED && (RUXD_ENABLED || RUXX_ENABLED))
#define CONJUGATE(X)GET_CONJG(X)
#else
#error  "Unrecognized interface."
#endif
        ! Define the diagonal elements of `cor` if needed.
#if     (RUXD_ENABLED || RXLD_ENABLED) && (VUXD_ENABLED || VXLD_ENABLED)
#define SET_CORDIA_BEG(DIM) cor(DIM, DIM) = 1._TKC
#define SET_CORDIA_END(DIM)
#elif   (RUXD_ENABLED || RXLD_ENABLED) && (VUXX_ENABLED || VXLX_ENABLED)
#define SET_CORDIA_BEG(DIM)
#define SET_CORDIA_END(DIM) cor(DIM, DIM) = 1._TKC
#elif   RUXX_ENABLED || RXLX_ENABLED
#define SET_CORDIA_BEG(DIM)
#define SET_CORDIA_END(DIM)
#else
#error  "Unrecognized interface."
#endif
        ! Define the inverse standard deviations if needed.
#if     VUXD_ENABLED || VXLD_ENABLED
#define SET_STDINV(DIM) stdinv(DIM) = 1._TKC / sqrt(GET_RE(cov(DIM, DIM)));
#define SET_RANGE(START,STOP) STOP,START, -1
        real(TKC) :: stdinv(size(cor, 1, IK))
        !do concurrent(idim = 1 : size(stdinv, 1, IK); stdinv(idim) = 1._TKC / sqrt(cov(idim, idim)); end do
#elif   VUXX_ENABLED || VXLX_ENABLED
#define SET_RANGE(START,STOP) START,STOP
#define SET_STDINV(DIM)
        CHECK_ASSERTION(__LINE__, size(cov, 1, IK) == size(stdinv, 1, IK), SK_"@setCor(): The condition `size(cov, 1) == size(stdinv)` must hold. size(cov, 1), size(stdinv) = "//getStr([size(cov, 1, IK), size(stdinv, 1, IK)]))
        CHECK_ASSERTION(__LINE__, all(0._TKC < stdinv), SK_"@setCor(): The condition `all(0. < stdinv)` must hold. stdinv = "//getStr(stdinv))
#else
#error  "Unrecognized interface."
#endif
        ndim = size(cov, 1, IK)
        CHECK_ASSERTION(__LINE__, ndim == size(cov, 2, IK), SK_"@setCor(): The condition `size(cov, 1) == size(cov, 2)` must hold. shape(cov) = "//getStr(shape(cov, IK)))
        CHECK_ASSERTION(__LINE__, all(ndim == shape(cor, IK)), SK_"@setCor(): The condition `all(size(cov, 1) == shape(cor))` must hold. size(cov, 1), shape(cor) = "//getStr([ndim, shape(cor, IK)]))
        do jdim = 1, ndim
            SET_STDINV(jdim)
            SET_CORDIA_BEG(jdim)
            do idim = SET_RANGE(1, jdim - 1)
                cor(RINDEX(idim, jdim)) = stdinv(idim) * stdinv(jdim) * CONJUGATE(cov(VINDEX(idim, jdim)))
            end do
            SET_CORDIA_END(jdim)
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setCor_ENABLED && Prs_ENABLED && XY_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_OF_SAMPLE :: cov(2,2)
#if     Org_ENABLED
        call setCov(cov, x, y WEIGHT_ARGS)
#elif   Avg_ENABLED
        call setCov(cov, mean, x, y WEIGHT_ARGS)
#else
#error  "Unrecognized interface."
#endif
        cor = cov(1,2) / sqrt(GET_RE(cov(1,1)) * GET_RE(cov(2,2)))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setCor_ENABLED && Prs_ENABLED && (UXD_ENABLED || XLD_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: idim, jdim, ndim
        real(TKC) :: stdinv(size(cor, 1, IK))
#if     Org_ENABLED
        call setCov(cor, subset, sample, dim WEIGHT_ARGS)
#elif   Avg_ENABLED
        call setCov(cor, subset, mean, sample, dim WEIGHT_ARGS)
#else
#error  "Unrecognized interface."
#endif
        ndim = size(cor, 1, IK)
        do idim = 1, ndim
            stdinv(idim) = 1._TKC / sqrt(GET_RE(cor(idim, idim)))
            cor(idim, idim) = 1._TKC
        end do
#if     UXD_ENABLED
#define THIS_RANGE 1, jdim - 1
#elif   XLD_ENABLED
#define THIS_RANGE jdim + 1, ndim
#else
#error  "Unrecognized interface."
#endif
        do jdim = 1, ndim
            do idim = THIS_RANGE
                cor(idim, jdim) = cor(idim, jdim) * stdinv(idim) * stdinv(jdim)
            end do
        end do
#undef  THIS_RANGE

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getRho_ENABLED && (XY_D0_ENABLED || XY_D1_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     XY_D0_ENABLED
        real(RK) :: frankx(len(x, IK)), franky(len(x, IK))
#elif   XY_D1_ENABLED
        real(RK) :: frankx(size(x, 1, IK)), franky(size(x, 1, IK))
#else
#error  "Unrecognized interface."
#endif
#if     WNO_ENABLED
        call setRho(rho, frankx, franky, x, y)
#elif   WTI_ENABLED || WTR_ENABLED
        call setRho(rho, frankx, franky, x, y, weight)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getRho_ENABLED && ULD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: idim, ndim
        real(RK) :: frank(size(sample, 1, IK), size(sample, 2, IK))
#if     WNO_ENABLED
        call setRho(rho, uppDia, frank, sample, dim)
#elif   WTI_ENABLED || WTR_ENABLED
        call setRho(rho, uppDia, frank, sample, dim, weight)
#else
#error  "Unrecognized interface."
#endif
        ndim = size(rho, 1, IK)
        ! The `do concurrent` version is problematic for OpenMP builds of the ParaMonte library with the Intel ifort 2021.8 on heap, yielding segmentation fault.
        !do concurrent(idim = 1 : ndim)
        do idim = 1, ndim
            rho(idim, 1 : idim - 1) = rho(1 : idim - 1, idim)
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setRho_ENABLED && (XY_D0_ENABLED || XY_D1_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RK) :: mean(2)
#if     WTI_ENABLED || WTR_ENABLED
        TYPE_OF_WEIGHT :: weisum
#elif   !WNO_ENABLED
#error  "Unrecognized interface."
#endif
#if     XY_D0_ENABLED
        CHECK_ASSERTION(__LINE__, size(frankx, 1, IK) == len(x, IK), SK_"@setRho(): The condition `size(frankx) == len(x)` must hold. size(frankx), size(x) = "//getStr([size(frankx, 1, IK), len(x, IK)]))
        CHECK_ASSERTION(__LINE__, size(franky, 1, IK) == len(y, IK), SK_"@setRho(): The condition `size(franky) == len(y)` must hold. size(franky), size(y) = "//getStr([size(franky, 1, IK), len(y, IK)]))
#elif   XY_D1_ENABLED
        CHECK_ASSERTION(__LINE__, size(frankx, 1, IK) == size(x, 1, IK), SK_"@setRho(): The condition `size(frankx) == size(x)` must hold. size(frankx), size(x) = "//getStr([size(frankx, 1, IK), size(x, 1, IK)]))
        CHECK_ASSERTION(__LINE__, size(franky, 1, IK) == size(y, 1, IK), SK_"@setRho(): The condition `size(franky) == size(y)` must hold. size(franky), size(y) = "//getStr([size(franky, 1, IK), size(y, 1, IK)]))
#else
#error  "Unrecognized interface."
#endif
        call setRankFractional(frankx, x)
        call setRankFractional(franky, y)
        call setMean(mean, frankx, franky WEIGHT_ARGS)
        call setCor(rho, mean, frankx, franky WEIGHT_ARGS)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setRho_ENABLED && (UXD_ENABLED || XLD_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     WTI_ENABLED || WTR_ENABLED
        TYPE_OF_WEIGHT :: weisum
#elif   !WNO_ENABLED
#error  "Unrecognized interface."
#endif
        integer(IK) :: idim, ndim, nsam
        real(RK) :: mean(size(sample, 3 - dim, IK))
        CHECK_ASSERTION(__LINE__, all(shape(frank, IK) == shape(sample, IK)), SK_"@setRho(): The condition `all(shape(frank) == shape(sample))` must hold. shape(frank), shape(sample) = "//getStr([shape(frank, IK), shape(sample, IK)]))
        ndim = size(rho, 3 - dim, IK)
        nsam = size(sample, dim, IK)
        if (dim == 1_IK) then
            do idim = 1, ndim
                call setRankFractional(frank(1:nsam, idim), sample(1:nsam, idim))
            end do
        else
            do idim = 1, ndim
                call setRankFractional(frank(idim, 1:nsam), sample(idim, 1:nsam))
            end do
        end if
        call setMean(mean, frank, dim WEIGHT_ARGS)
        call setCor(rho, subset, mean, frank, dim WEIGHT_ARGS)

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef TYPE_OF_SAMPLE
#undef TYPE_OF_WEIGHT
#undef SET_CORDIA_BEG
#undef SET_CORDIA_END
#undef WEIGHT_ARGS
#undef SET_STDINV
#undef CONJUGATE
#undef GET_CONJG
#undef SET_RANGE
#undef GET_RE
#undef RINDEX
#undef VINDEX