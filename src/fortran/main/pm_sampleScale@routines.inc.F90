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
!>  This file contains the implementation details of the routines under the generic interfaces of [pm_sampleScale](@ref pm_sampleScale).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Saturday 2:33 AM, August 22, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the runtime checks.
#define CHECK_VAL_DIM \
CHECK_ASSERTION(__LINE__, 1 <= dim .and. dim <= rank(sample), \
SK_"@setScaled(): The condition `1 <= dim .and. dim <= rank(sample)` must hold. dim, rank(sample) = "//getStr([integer(IK) :: dim, rank(sample)]))
#define CHECK_LEN_AMOUNT(DIM) \
CHECK_ASSERTION(__LINE__, size(sample, DIM, IK) == size(amount, 1, IK), \
SK_"@setScaled(): The condition `size(sample, dim, 1) == size(amount, 1)` must hold. dim, shape(sample), size(amount) = "\
//getStr([DIM, shape(sample, IK), size(amount, 1, IK)]))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getScaled_ENABLED && D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        sampleScaled = sample * amount

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setScaled_ENABLED && D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        sample = sample * amount

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getScaled_ENABLED && D2_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        ! Set the indexing rule.
#if     ONO_ENABLED
#define GET_INDEX(I,J)I,J
#elif   OTH_ENABLED
#define GET_INDEX(I,J)J,I
#else
#error  "Unrecognized interface."
#endif
        ! Set the conjugation rule.
#if     OTH_ENABLED && CK_ENABLED
#define GET_CONJG(X)conjg(X)
#elif   ONO_ENABLED || (OTH_ENABLED && RK_ENABLED)
#define GET_CONJG(X)X
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: idim, isam, ndim, nsam
        CHECK_VAL_DIM
        CHECK_LEN_AMOUNT(3 - dim)
        nsam = size(sample, dim, IK)
        ndim = size(sample, 3 - dim, IK)
        if (dim == 2_IK) then
            do concurrent(idim = 1 : ndim, isam = 1 : nsam)
                sampleScaled(GET_INDEX(idim, isam)) = GET_CONJG(sample(idim, isam) * amount(idim))
            end do
        else
            do concurrent(idim = 1 : ndim, isam = 1 : nsam)
                sampleScaled(GET_INDEX(isam, idim)) = GET_CONJG(sample(isam, idim) * amount(idim))
            end do
        end if
#undef  GET_CONJG
#undef  GET_INDEX

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setScaled_ENABLED && D2_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ell, ndim, nsam
        CHECK_VAL_DIM
        CHECK_LEN_AMOUNT(3 - dim)
        nsam = size(sample, dim, IK)
        ndim = size(sample, 3 - dim, IK)
        if (dim == 2_IK) then
            do concurrent(ell = 1 : nsam)
                sample(1 : ndim, ell) = sample(1 : ndim, ell) * amount
            end do
        else
            do concurrent(ell = 1 : ndim)
                sample(1 : nsam, ell) = sample(1 : nsam, ell) * amount(ell)
            end do
        end if

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  CHECK_LEN_AMOUNT
#undef  CHECK_VAL_DIM