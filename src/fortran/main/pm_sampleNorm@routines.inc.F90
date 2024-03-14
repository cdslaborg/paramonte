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
!>  This file contains the implementation details of the routines under the generic interfaces of [pm_sampleNorm](@ref pm_sampleNorm).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Saturday 2:33 AM, August 22, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the runtime checks.
#define CHECK_VAL_DIM \
CHECK_ASSERTION(__LINE__, 1 <= dim .and. dim <= rank(sample), \
SK_"@setNormed(): The condition `1 <= dim .and. dim <= rank(sample)` must hold. dim, rank(sample) = "//getStr([integer(IK) :: dim, rank(sample)]))
#define CHECK_LEN_SCALE(DIM) \
CHECK_ASSERTION(__LINE__, size(sample, DIM, IK) == size(scale, 1, IK), \
SK_"@setNormed(): The condition `size(sample, dim, 1) == size(scale, 1)` must hold. dim, shape(sample), size(scale) = "\
//getStr([DIM, shape(sample, IK), size(scale, 1, IK)]))
#define CHECK_LEN_SHIFT(DIM) \
CHECK_ASSERTION(__LINE__, size(sample, DIM, IK) == size(shift, 1, IK), \
SK_"@setNormed(): The condition `size(sample, dim, 1) == size(shift, 1)` must hold. dim, shape(sample), size(shift) = "\
//getStr([DIM, shape(sample, IK), size(shift, 1, IK)]))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getNormed_ENABLED && D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        sampleNormed = (sample + shift) * scale

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setNormed_ENABLED && D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        sample = (sample + shift) * scale

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getNormed_ENABLED && D2_ENABLED
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
        CHECK_LEN_SCALE(3 - dim)
        CHECK_LEN_SHIFT(3 - dim)
        nsam = size(sample, dim, IK)
        ndim = size(sample, 3 - dim, IK)
        if (dim == 2_IK) then
            do concurrent(idim = 1 : ndim, isam = 1 : nsam)
                sampleNormed(GET_INDEX(idim, isam)) = GET_CONJG((sample(idim, isam) + shift(idim)) * scale(idim))
            end do
        else
            do concurrent(idim = 1 : ndim, isam = 1 : nsam)
                sampleNormed(GET_INDEX(isam, idim)) = GET_CONJG((sample(isam, idim) + shift(idim)) * scale(idim))
            end do
        end if
#undef  GET_CONJG
#undef  GET_INDEX

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setNormed_ENABLED && D2_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: idim, isam, ndim, nsam
        CHECK_VAL_DIM
        CHECK_LEN_SCALE(3 - dim)
        CHECK_LEN_SHIFT(3 - dim)
        nsam = size(sample, dim, IK)
        ndim = size(sample, 3 - dim, IK)
        if (dim == 2_IK) then
            do concurrent(isam = 1 : nsam, idim = 1 : ndim)
                sample(idim, isam) = (sample(idim, isam) + shift(idim)) * scale(idim)
            end do
        else
            do concurrent(idim = 1 : ndim, isam = 1 : nsam)
                sample(isam, idim) = (sample(isam, idim) + shift(idim)) * scale(idim)
            end do
        end if

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  CHECK_LEN_SCALE
#undef  CHECK_LEN_SHIFT
#undef  CHECK_VAL_DIM