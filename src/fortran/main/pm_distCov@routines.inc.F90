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
!>  This include file contains the implementation of procedures in [pm_distCov](@ref pm_distCov).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getCovRand_ENABLED && GRAM_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        type(rngf_type) :: rng
#if     S0_ENABLED
#if     CK_ENABLED
        complex(TKG), parameter :: ONE = (1._TKG, 0._TKG)
#elif   RK_ENABLED
        real(TKG), parameter :: ONE = 1._TKG
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: jdim
        if (present(scale)) then
            call setCovRand(rng, rand, scale)
        else
            call setCovRand(rng, rand)
            ! Ensure all diagonals are 1.
            do jdim = 1, sizE(rand, 1, IK)
                rand(jdim, jdim) = ONE
            end do
        end if
#elif   S1_ENABLED
        call setCovRand(rng, rand, scale)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getCovRand_ENABLED && (DVINE_ENABLED || ONION_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        type(rngf_type) :: rng
        ! The onion method appears unstable, use the slower but more stable method for now.
        !type(onion_type) :: onion
#if     S0_ENABLED
        if (present(scale)) then
            call setCovRand(rng, rand, method, eta, scale)
        else
            call setCovRand(rng, rand, method, eta)
        end if
#elif   S1_ENABLED
        call setCovRand(rng, rand, method, eta, scale)
#else
#error  "Unrecognized interface."
#endif
#if     ONION_ENABLED
        if (onion%info /= 0_IK) error stop "@getCovRand(): The Cholesky factorization of the onion method failed."
#endif
        !call setMatCopy(rand(2:, 1:), rdpack, rand(1:, 2:), rdpack, uppDia)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getCovRand_ENABLED && GRAM_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#error  "Unrecognized interface."
!        type(rngf_type) :: rng
!#if     S0_ENABLED
!        if (present(scale)) then
!            call setCovRand(rng, rand, scale)
!        else
!            call setCovRand(rng, rand)
!        end if
!#elif   S1_ENABLED
!        call setCovRand(rng, rand, scale)
!#else
!#error  "Unrecognized interface."
!#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setCovRand_ENABLED && GRAM_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKG) :: normfac
        !real(TKG), parameter :: EPS = 2 * sqrt(epsilon(0._TKG))
        ! Define the output kind.
#if     CK_ENABLED
#define TYPE_OF_RAND complex(TKG)
#define GET_CONJG(X)conjg(X)
#define GET_RE(X)X%re
        complex(TKG), parameter :: LB = -cmplx(1._TKG, 1._TKG, TKG), UB = cmplx(1._TKG, 1._TKG, TKG), ZERO = (0._TKG, 0._TKG), ONE = (1._TKG, 0._TKG)
#elif   RK_ENABLED
#define GET_RE(X)X
#define GET_CONJG(X)X
#define TYPE_OF_RAND real(TKG)
        real(TKG), parameter :: LB = -1._TKG, UB = 1._TKG, ZERO = 0._TKG, ONE = 1._TKG
#else
#error  "Unrecognized interface."
#endif
        TYPE_OF_RAND :: upper(size(rand, 1, IK), size(rand, 2, IK))
        integer(IK) :: idim, jdim, ndim
#if     SD_ENABLED
#define GET_SCALED(X)X
#elif   S0_ENABLED
#define GET_SCALED(X)X * scaledSq
        real(TKG) :: scaledSq
        scaledSq = scale**2
        CHECK_ASSERTION(__LINE__, 0._TKG < scale, SK_"@setCovRand(): The condition `0. < scale` must hold. scale = "//getStr(scale))
#elif   S1_ENABLED
#define GET_SCALED(X)X * scale(idim) * scale(jdim)
        CHECK_ASSERTION(__LINE__, all([0._TKG < scale]), SK_"@setCovRand(): The condition `all([0. < scale])` must hold. scale = "//getStr(scale))
        CHECK_ASSERTION(__LINE__, size(scale, 1, IK) == size(rand, 1, IK), SK_"@setCovRand(): The condition `size(scale) == size(rand, 1, IK)` must hold. size(scale), shape(rand) = "//getStr([size(scale, 1, IK), shape(rand, IK)]))
#else
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, size(rand, 1, IK) == size(rand, 2, IK), SK_"@setCovRand(): The condition `size(rand, 1) == size(rand, 2)` must hold. shape(rand) = "//getStr(shape(rand, IK)))
        ndim = size(rand, 1, IK)
        if (ndim < 1_IK) return
        do idim = 1, ndim
            do
                call setUnifRand(rng, upper(1 : idim, idim), LB, UB)
                normfac = real(sqrt(dot_product(upper(1 : idim, idim), upper(1 : idim, idim))), TKG) ! \todo: The performance of this expression can be improved by replacing `dot_product` with `absq()`.
                if (ZERO == normfac) cycle
                exit
            end do
            normfac = 1._TKG / normfac
            upper(1 : idim, idim) = upper(1 : idim, idim) * normfac
        end do
        rand = matmul(transpose(GET_CONJG(upper)), upper)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Ensure the diagonals are all pure 1, free from numerical round-off errors.
        ! Perhaps not really essential as it is guaranteed to be 1, at least theoretically.
        ! This is now done only within `getCovRand` above.
        !do jdim = 1, ndim
        !    rand(jdim, jdim) = ONE
        !end do
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !block
        !   use pm_io
        !   call disp%show("rand")
        !   call disp%show( rand )
        !end block
        !call setCor(rand, uppDia, rand, uppDia)

        ! Rescale.
#if     S0_ENABLED || S1_ENABLED
        do jdim = 1, ndim
            do idim = 1, jdim
                rand(idim, jdim) = GET_SCALED(rand(idim, jdim))
            end do
        end do
#endif
        !block
        !use pm_io, only: disp
        !call disp%show("rand")
        !call disp%show( rand )
        !end block
        ! copy to the lower triangle.
        call setMatCopy(rand, rdpack, rand, rdpack, upp, transHerm)
#undef  TYPE_OF_RAND
#undef  GET_SCALED
#undef  GET_CONJG
#undef  GET_RE

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setCovRand_ENABLED && (DVINE_ENABLED || ONION_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GET_INDEX(i,j) i, j
!#if     UXD_ENABLED
!#define GET_INDEX(i,j) j, i
!#elif   XLD_ENABLED
!#define GET_INDEX(i,j) i, j
!#else
!#error  "Unrecognized interface."
!#endif
        real(TKG) :: beta
        integer(IK) :: ndim, idim, jdim, kdim
        ! Set the scale
#if     SD_ENABLED
        real(TKG), parameter :: scaleSq = 1._TKG
#define GET_SCALED(x, i, j)
#else
#if     S0_ENABLED
        real(TKG) :: scaleSq
        scaleSq = scale * scale
#define GET_SCALED(x, i, j) x = x * scaleSq
#elif   S1_ENABLED
        CHECK_ASSERTION(__LINE__, size(scale, 1, IK) == size(rand, 1, IK), SK_"@setCovRand(): The condition `size(scale) == size(rand, 1, IK)` must hold. size(scale), shape(rand) = "//getStr([size(scale, 1, IK), shape(rand, IK)]))
#define GET_SCALED(x, i, j) x = x * scale(i) * scale(j)
#else
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, all([0._TKG < scale]), SK_"@setCovRand(): The condition `all([0. < scale])` must hold. scale = "//getStr(scale))
#endif
        CHECK_ASSERTION(__LINE__, 0._TKG <= eta, SK_"@setCovRand(): The condition `0. <= eta` must hold. eta = "//getStr(eta))
       !CHECK_ASSERTION(__LINE__, 0_IK < size(rand, 1, IK), SK_"@setCovRand(): The condition `0 < size(rand, 1)` must hold. shape(rand) = "//getStr(shape(rand)))
        CHECK_ASSERTION(__LINE__, size(rand, 1, IK) == size(rand, 2, IK), SK_"@setCovRand(): The condition `size(rand, 1) <= size(rand, 2)` must hold. shape(rand) = "//getStr(shape(rand)))
        ndim = size(rand, 1, IK)
        if (1_IK < ndim) then
#if         ONION_ENABLED
            beta = eta + 0.5_TKG * ndim
            block
                integer(IK) :: kdimp1
                real(TKG) :: brand, ssnrandInv, chol(ndim - 1, ndim - 1), wrand(ndim - 1) ! sphere rv.
                do
                    call setBetaRand(rng, brand, beta, beta)
                    if (0._TKG < brand .and. brand < 1._TKG) exit
                end do
                rand(1, 1) = 1._TKG
                rand(2, 2) = 1._TKG
                rand(1, 2) = 2._TKG * brand - 1._TKG
                do kdim = 2_IK, ndim - 1_IK
                    call setNormRand(rng, wrand(1 : kdim))
                    ssnrandInv = 1._TKG / norm2(wrand(1 : kdim)) ! `wrand(1 : kdim) * ssnrandInv` would be a uniform random point on the surface of the kdim-sphere.
                    beta = beta - 0.5_TKG
                    do
                        call setBetaRand(rng, brand, 0.5_TKG * kdim, beta)
                        if (0._TKG < brand .and. brand < 1._TKG) exit
                    end do
                    ssnrandInv = ssnrandInv * sqrt(brand)
                    call setMatChol(rand, uppDia, method%info, chol, nothing, kdim, 0_IK, 0_IK, 0_IK, 0_IK)
                    if (method%info /= 0_IK) return
                    ! Compute z = cholow * wrand.
                    kdimp1 = kdim + 1
                    rand(1 : kdim, kdimp1) = 0._TKG
                    do jdim = 1, kdim
                        rand(jdim, kdimp1) = rand(jdim, kdimp1) + chol(jdim, jdim) * wrand(jdim) * ssnrandInv
                        do idim = jdim + 1, kdim
                            rand(idim, kdimp1) = rand(idim, kdimp1) + rand(idim, jdim) * rand(jdim, kdimp1)
                        end do
                    end do
                    rand(kdimp1, kdimp1) = 1._TKG
                    !block
                        !use pm_io, only: display_type
                        !type(display_type) :: disp
                        !call disp%skip
                        !call disp%show(rand(1 : kdimp1, 1 : kdimp1))
                        !call disp%skip
                    !end block
                end do
                do kdim = ndim, 1, -1
#if                 SD_ENABLED || S0_ENABLED
                    rand(kdim, kdim) = scaleSq
#elif               S1_ENABLED
                    rand(kdim, kdim) = rand(kdim, kdim) * scale(kdim) * scale(kdim)
#endif
                    do idim = kdim + 1, ndim
#if                     SD_ENABLED || S0_ENABLED
                        rand(kdim, idim) = rand(kdim, idim) * scaleSq
#elif                   S1_ENABLED
                        rand(kdim, idim) = rand(kdim, idim) * scale(idim) * scale(kdim)
#endif
                        rand(idim, kdim) = rand(kdim, idim)
                    end do
                end do
                !call setMatCopy(rand(2 : ndim, 1 : ndim), rdpack, rand(1 : ndim, 2 : ndim), rdpack, uppDia)
            end block
            !block
                !use pm_io, only: display_type
                !type(display_type) :: disp
                !call disp%skip
                !call disp%show(rand)
                !call disp%skip
            !end block
#elif       DVINE_ENABLED
            beta = eta + 0.5_TKG * (ndim + 1_IK)
            block
                real(TKG) :: parCorMat(size(rand, 1, IK), size(rand, 1, IK))
                parCorMat = 0._TKG
                do kdim = 1_IK, ndim - 1_IK
                    beta = beta - 0.5_TKG
#if                 SD_ENABLED || S0_ENABLED
                    rand(kdim, kdim) = scaleSq
#elif               S1_ENABLED
                    rand(kdim, kdim) = scale(kdim) * scale(kdim)
#endif
                    do idim = kdim + 1_IK, ndim
                        loopSensibleBeta: do
                            call setBetaRand(rng, parCorMat(kdim, idim), beta, beta)
                            if (0._TKG < parCorMat(kdim, idim) .and. parCorMat(kdim, idim) < 1._TKG) then
                                parCorMat(kdim, idim) = 2 * parCorMat(kdim, idim) - 1._TKG ! Linearly shift to the range (-1, 1).
                                ! Convert the partial correlations to full correlation.
                                rand(GET_INDEX(kdim, idim)) = parCorMat(kdim, idim)
                                do jdim = kdim - 1_IK, 1_IK, -1_IK
                                    rand(GET_INDEX(kdim, idim)) = parCorMat(jdim, idim) * parCorMat(jdim, kdim) + & ! LCOV_EXCL_LINE
                                    rand(GET_INDEX(kdim, idim)) * sqrt((1._TKG - parCorMat(jdim, idim)**2) * (1._TKG - parCorMat(jdim, kdim)**2))
                                end do
                                GET_SCALED(rand(GET_INDEX(kdim, idim)), idim, kdim)
                                rand(GET_INDEX(idim, kdim)) = rand(GET_INDEX(kdim, idim))
                                exit loopSensibleBeta
                            end if
                        end do loopSensibleBeta
                    end do
                end do
#if             SD_ENABLED || S0_ENABLED
                rand(kdim, kdim) = scaleSq
#elif           S1_ENABLED
                rand(kdim, kdim) = scale(kdim) * scale(kdim)
#endif
            end block
#else
#error  "Unrecognized interface."
#endif
        elseif (1_IK == ndim) then
#if         SD_ENABLED || S0_ENABLED
            rand(1,1) = scaleSq
#elif       S1_ENABLED
            rand(1,1) = scale(1) * scale(1)
#endif
        end if
#if     ONION_ENABLED
        method%info = 0_IK
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  GET_SCALED
#undef  GET_INDEX