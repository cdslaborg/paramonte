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
!>  This include file contains the implementation of procedures in [pm_distChol](@ref pm_distChol).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the output kind.
#if     CK_ENABLED
#define TYPE_OF_RAND complex(TKG)
        complex(TKG), parameter :: LB = -cmplx(1._TKG, 1._TKG, TKG), UB = cmplx(1._TKG, 1._TKG, TKG), ZERO = (0._TKG, 0._TKG), ONE = (1._TKG, 0._TKG)
#elif   RK_ENABLED
#define TYPE_OF_RAND real(TKG)
        real(TKG), parameter :: LB = -1._TKG, UB = 1._TKG, ZERO = 0._TKG, ONE = 1._TKG
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%
#if     getCholRand_ENABLED
        !%%%%%%%%%%%%%%%%%%

        type(rngf_type) :: rng
        call setCholRand(rng, rand, subset)
        select type (subset)
        type is (uppDia_type)
            call setMatInit(rand, low, vlow = ZERO)
        type is (lowDia_type)
            call setMatInit(rand, upp, vupp = ZERO)
        class default
            error stop "Unrecognized unsupported value specified for the input argument `subset`."
        end select

        !%%%%%%%%%%%%%%%%%%
#elif   setCholRand_ENABLED
        !%%%%%%%%%%%%%%%%%%

        real(TKG) :: normfac
        integer(IK) :: idim, ndim
        CHECK_ASSERTION(__LINE__, size(rand, 1, IK) == size(rand, 2, IK), SK_"@setCholRand(): The condition `size(rand, 1) == size(rand, 2)` must hold. shape(rand) = "//getStr(shape(rand, IK)))
        ndim = size(rand, 1, IK)
        if (ndim < 1_IK) return
        select type (subset)
        type is (uppDia_type)
            do idim = 1, ndim
                do
                    call setUnifRand(rng, rand(1 : idim, idim), LB, UB)
#if                 CK_ENABLED
                    rand(idim, idim)%im = 0._TKG
                    rand(idim, idim)%re = abs(rand(idim, idim)%re)
                    if (rand(idim, idim)%re < epsilon(0._TKG)) rand(idim, idim)%re = 1._TKG
#elif               RK_ENABLED
                    rand(idim, idim) = abs(rand(idim, idim))
                    if (rand(idim, idim) < epsilon(0._TKG)) rand(idim, idim) = 1._TKG
#endif
                    normfac = real(sqrt(dot_product(rand(1 : idim, idim), rand(1 : idim, idim))), TKG)
                    ! \todo: The performance of the above expression can be improved by replacing `dot_product` with `absq()`.
                    if (0._TKG == normfac) cycle
                    exit
                end do
                normfac = 1._TKG / normfac
                rand(1 : idim, idim) = rand(1 : idim, idim) * normfac
            end do
        type is (lowDia_type)
            do idim = 1, ndim
                do
                    call setUnifRand(rng, rand(idim, 1 : idim), LB, UB)
#if                 CK_ENABLED
                    rand(idim, idim)%im = 0._TKG
                    rand(idim, idim)%re = abs(rand(idim, idim)%re)
                    if (rand(idim, idim)%re < epsilon(0._TKG)) rand(idim, idim)%re = 1._TKG
#elif               RK_ENABLED
                    rand(idim, idim) = abs(rand(idim, idim))
                    if (rand(idim, idim) < epsilon(0._TKG)) rand(idim, idim) = 1._TKG
#endif
                    normfac = real(sqrt(dot_product(rand(idim, 1 : idim), rand(idim, 1 : idim))), TKG)
                    ! \todo: The performance of the above expression can be improved by replacing `dot_product` with `absq()`.
                    if (0._TKG == normfac) cycle
                    exit
                end do
                normfac = 1._TKG / normfac
                rand(idim, 1 : idim) = rand(idim, 1 : idim) * normfac
            end do
        class default
            error stop "Unrecognized unsupported value specified for the input argument `subset`."
        end select

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  TYPE_OF_RAND