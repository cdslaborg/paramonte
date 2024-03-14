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
!>  This file contains the implementation details of the routines under the generic interfaces of [pm_sampleShift](@ref pm_sampleShift).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Saturday 2:33 AM, August 22, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the runtime checks.
#define CHECK_VAL_DIM \
CHECK_ASSERTION(__LINE__, 1 <= dim .and. dim <= rank(sample), \
SK_"@setShifted(): The condition `1 <= dim .and. dim <= rank(sample)` must hold. dim, rank(sample) = "//getStr([integer(IK) :: dim, rank(sample)]))
#define CHECK_LEN_AMOUNT(DIM) \
CHECK_ASSERTION(__LINE__, size(sample, DIM, IK) == size(amount, 1, IK), \
SK_"@setShifted(): The condition `size(sample, dim, 1) == size(amount, 1)` must hold. dim, shape(sample), size(amount) = "\
//getStr([DIM, shape(sample, IK), size(amount, 1, IK)]))
        ! Set the indexing rule.
#if     ONO_ENABLED
#define GET_INDEX(I,J)I,J
#elif   OTH_ENABLED
#define GET_INDEX(I,J)J,I
#elif   getShifted_ENABLED && D2_ENABLED
#error  "Unrecognized interface."
#endif
        ! Set the conjugation rule.
#if     OTH_ENABLED && CK_ENABLED
#define GET_CONJG(X)conjg(X)
#elif   ONO_ENABLED || (OTH_ENABLED && RK_ENABLED)
#define GET_CONJG(X)X
#elif   getShifted_ENABLED && D2_ENABLED
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getShifted_ENABLED && D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        sampleShifted = sample + amount

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setShifted_ENABLED && (D1_ENABLED || (D2_ENABLED && ALL_ENABLED))
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        sample = sample + amount

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getShifted_ENABLED && D2_ENABLED && ALL_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: idim, jdim, ndim, mdim
        ndim = size(sample, 1, IK)
        mdim = size(sample, 2, IK)
        do concurrent(idim = 1 : ndim, jdim = 1 : mdim)
#if         ONO_ENABLED
            sampleShifted(idim, jdim) = GET_CONJG(sample(idim, jdim) + amount)
#elif       OTH_ENABLED
            sampleShifted(jdim, idim) = GET_CONJG(sample(idim, jdim) + amount)
#else
#error      "Unrecognized interface."
#endif
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getShifted_ENABLED && D2_ENABLED && DIM_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: idim, isam, ndim, nsam
        CHECK_VAL_DIM
        CHECK_LEN_AMOUNT(3 - dim)
        nsam = size(sample, dim, IK)
        ndim = size(sample, 3 - dim, IK)
        if (dim == 2_IK) then
            do concurrent(idim = 1 : ndim, isam = 1 : nsam)
                sampleShifted(GET_INDEX(idim, isam)) = GET_CONJG(sample(idim, isam) + amount(idim))
            end do
        else
            do concurrent(idim = 1 : ndim, isam = 1 : nsam)
                sampleShifted(GET_INDEX(isam, idim)) = GET_CONJG(sample(isam, idim) + amount(idim))
            end do
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setShifted_ENABLED && D2_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ell, ndim, nsam
        CHECK_VAL_DIM
        CHECK_LEN_AMOUNT(3 - dim)
        nsam = size(sample, dim, IK)
        ndim = size(sample, 3 - dim, IK)
        if (dim == 2_IK) then
            do concurrent(ell = 1 : nsam)
                sample(1 : ndim, ell) = sample(1 : ndim, ell) + amount
            end do
        else
            do concurrent(ell = 1 : ndim)
                sample(1 : nsam, ell) = sample(1 : nsam, ell) + amount(ell)
            end do
        end if

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  CHECK_LEN_AMOUNT
#undef  CHECK_VAL_DIM
#undef  GET_CONJG
#undef  GET_INDEX