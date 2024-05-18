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
!>  This file contains the implementation details of the routines under the generic interfaces of [pm_sampleAffinity](@ref pm_sampleAffinity).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Saturday 2:33 AM, August 22, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the default mean.
#if     ATL_ENABLED
#define CHECK_LEN_TLATE(DIM) \
CHECK_ASSERTION(__LINE__, size(sample, DIM, IK) == size(tlate, 1, IK), \
SK_"@setAffinity(): The condition `size(sample, dim, 1) == size(tlate, 1)` must hold. dim, shape(sample), size(tlate) = "\
//getStr([DIM, shape(sample, IK), size(tlate, 1, IK)]))
#define PLUS(X)+ X
#elif   DTL_ENABLED
#define CHECK_LEN_TLATE(DIM)
#define PLUS(X)
#elif   setAffinity_ENABLED
#error  "Unrecognized interface."
#endif
        ! Define the runtime checks.
#define CHECK_VAL_DIM \
CHECK_ASSERTION(__LINE__, 1 <= dim .and. dim <= rank(sample), \
SK_"@setAffinity(): The condition `1 <= dim .and. dim <= rank(sample)` must hold. dim, rank(sample) = "//getStr([integer(IK) :: dim, rank(sample)]))
#define CHECK_SHAPE_SAMAFF \
CHECK_ASSERTION(__LINE__, all(shape(affinity, IK) == shape(sample, IK)), \
SK_"@setAffinity(): The condition `all(shape(affinity) == shape(sample))` must hold. shape(affinity), shape(sample) = "\
//getStr([shape(affinity, IK), shape(sample, IK)]))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getAffinity_ENABLED && D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CALL_SET_AFFINE(CLASS) \
if (present(tlate)) then; \
call setAffinity(affinity, sample, tform, CLASS, tlate); \
else; \
call setAffinity(affinity, sample, tform, CLASS); \
end if;
        if (present(class)) then
            select type (class)
            type is (genrecmat_type)
                CALL_SET_AFFINE(class)
            type is (upperDiag_type)
                CALL_SET_AFFINE(class)
            type is (lowerDiag_type)
                CALL_SET_AFFINE(class)
            type is (upperUnit_type)
                CALL_SET_AFFINE(class)
            type is (lowerUnit_type)
                CALL_SET_AFFINE(class)
            class default
                error stop "Unrecognized `class`."
            end select
            return
        end if
        CALL_SET_AFFINE(genrecmat)
#undef  CALL_SET_AFFINE

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getAffinity_ENABLED && D2_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CALL_SET_AFFINE(CLASS) \
if (present(tlate)) then; \
call setAffinity(affinity, sample, dim, tform, CLASS, tlate); \
else; \
call setAffinity(affinity, sample, dim, tform, CLASS); \
end if;
        if (present(class)) then
            select type (class)
            type is (genrecmat_type)
                CALL_SET_AFFINE(class)
            type is (upperDiag_type)
                CALL_SET_AFFINE(class)
            type is (lowerDiag_type)
                CALL_SET_AFFINE(class)
            type is (upperUnit_type)
                CALL_SET_AFFINE(class)
            type is (lowerUnit_type)
                CALL_SET_AFFINE(class)
            class default
                error stop "Unrecognized `class`."
            end select
            return
        end if
        CALL_SET_AFFINE(genrecmat)
#undef  CALL_SET_AFFINE

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setAffinity_ENABLED && D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     !CGR_ENABLED
        integer(IK) :: idim
#endif
        integer(IK) :: ndim
        ndim = size(sample, 1, IK)
        CHECK_LEN_TLATE(1)
        CHECK_SHAPE_SAMAFF
        if (ndim < 1_IK) return
#if     CGR_ENABLED
        affinity = matmul(tform, sample) PLUS(tlate)
#elif   CUD_ENABLED
        !do idim = 1, ndim
        !    affinity(idim) = dot_product(sample(idim : ndim), tform(idim, idim : ndim)) PLUS(tlate(idim))
        !end do
        affinity(1 : ndim) = sample(ndim) * tform(1 : ndim, ndim) PLUS(tlate)
        do idim = ndim - 1, 1, -1
            affinity(1 : idim) = affinity(1 : idim) + sample(idim) * tform(1 : idim, idim)
        end do
#elif   CUU_ENABLED
        !do idim = 1, ndim
        !    affinity(idim) = dot_product(sample(idim : ndim), tform(idim, idim : ndim)) PLUS(tlate(idim))
        !end do
        affinity(1 : ndim - 1) = sample(ndim) * tform(1 : ndim - 1, ndim) PLUS(tlate(1 : ndim - 1))
        affinity(ndim) = sample(ndim) PLUS(tlate(ndim))
        do idim = ndim - 1, 1, -1
            affinity(1 : idim - 1) = affinity(1 : idim - 1) + sample(idim) * tform(1 : idim - 1, idim)
            affinity(idim) = affinity(idim) + sample(idim) * tform(idim, idim)
        end do
#elif   CLD_ENABLED
        affinity(1 : ndim) = sample(1) * tform(1 : ndim, 1) PLUS(tlate)
        do idim = 2, ndim
            affinity(idim : ndim) = affinity(idim : ndim) + sample(idim) * tform(idim : ndim, idim)
        end do
#elif   CLU_ENABLED
        affinity(1) = sample(1) PLUS(tlate(1))
        affinity(2 : ndim) = sample(1) * tform(2 : ndim, 1) PLUS(tlate(2 : ndim))
        do idim = 2, ndim
            affinity(idim) = affinity(idim) + sample(idim)
            affinity(idim + 1 : ndim) = affinity(idim + 1 : ndim) + sample(idim) * tform(idim + 1 : ndim, idim)
        end do
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setAffinity_ENABLED && D2_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     !CGR_ENABLED
        integer(IK) :: idim
#endif
        integer(IK) :: ndim, isam, nsam
        CHECK_VAL_DIM
        CHECK_SHAPE_SAMAFF
        CHECK_LEN_TLATE(3 - dim)
        nsam = size(sample, dim, IK)
        ndim = size(sample, 3 - dim, IK)
        if (ndim < 1_IK) return
        if (dim == 2_IK) then
            do concurrent(isam = 1 : nsam)
#if             CGR_ENABLED
                affinity(1 : ndim, isam) = matmul(tform, sample(1 : ndim, isam)) PLUS(tlate)
#elif           CUD_ENABLED
                !do idim = 1, ndim
                !    affinity(idim, isam) = dot_product(sample(idim : ndim, isam), tform(idim, idim : ndim)) PLUS(tlate(idim))
                !end do
                affinity(1 : ndim, isam) = sample(ndim, isam) * tform(1 : ndim, ndim) PLUS(tlate)
                do idim = ndim - 1, 1, -1
                    affinity(1 : idim, isam) = affinity(1 : idim, isam) + sample(idim, isam) * tform(1 : idim, idim)
                end do
#elif           CUU_ENABLED
                !do idim = 1, ndim
                !    affinity(idim, isam) = dot_product(sample(idim : ndim, isam), tform(idim, idim : ndim)) PLUS(tlate(idim))
                !end do
                affinity(1 : ndim - 1, isam) = sample(ndim, isam) * tform(1 : ndim - 1, ndim) PLUS(tlate(1 : ndim - 1))
                affinity(ndim, isam) = sample(ndim, isam) PLUS(tlate(ndim))
                do idim = ndim - 1, 1, -1
                    affinity(1 : idim - 1, isam) = affinity(1 : idim - 1, isam) + sample(idim, isam) * tform(1 : idim - 1, idim)
                    affinity(idim, isam) = affinity(idim, isam) + sample(idim, isam) * tform(idim, idim)
                end do
#elif           CLD_ENABLED
                affinity(1 : ndim, isam) = sample(1, isam) * tform(1 : ndim, 1) PLUS(tlate)
                do idim = 2, ndim
                    affinity(idim : ndim, isam) = affinity(idim : ndim, isam) + sample(idim, isam) * tform(idim : ndim, idim)
                end do
#elif           CLU_ENABLED
                affinity(1, isam) = sample(1, isam) PLUS(tlate(1))
                affinity(2 : ndim, isam) = sample(1, isam) * tform(2 : ndim, 1) PLUS(tlate(2 : ndim))
                do idim = 2, ndim
                    affinity(idim, isam) = affinity(idim, isam) + sample(idim, isam)
                    affinity(idim + 1 : ndim, isam) = affinity(idim + 1 : ndim, isam) + sample(idim, isam) * tform(idim + 1 : ndim, idim)
                end do
#else
#error          "Unrecognized interface."
#endif
            end do
        else
            ! \todo:
            ! The performance of this block can be improved by looping over `nsam` in the innermost layer.
            do concurrent(isam = 1 : nsam)
#if             CGR_ENABLED
                affinity(isam, 1 : ndim) = matmul(tform, sample(isam, 1 : ndim)) PLUS(tlate)
#elif           CUD_ENABLED
                affinity(isam, 1 : ndim) = sample(isam, ndim) * tform(1 : ndim, ndim) PLUS(tlate)
                do idim = ndim - 1, 1, -1
                    affinity(isam, 1 : idim) = affinity(isam, 1 : idim) + sample(isam, idim) * tform(1 : idim, idim)
                end do
#elif           CUU_ENABLED
                affinity(isam, 1 : ndim - 1) = sample(isam, ndim) * tform(1 : ndim - 1, ndim) PLUS(tlate(1 : ndim - 1))
                affinity(isam, ndim) = sample(isam, ndim) PLUS(tlate(ndim))
                do idim = ndim - 1, 1, -1
                    affinity(isam, 1 : idim - 1) = affinity(isam, 1 : idim - 1) + sample(isam, idim) * tform(1 : idim - 1, idim)
                    affinity(isam, idim) = affinity(isam, idim) + sample(isam, idim) * tform(idim, idim)
                end do
#elif           CLD_ENABLED
                affinity(isam, 1 : ndim) = sample(isam, 1) * tform(1 : ndim, 1) PLUS(tlate)
                do idim = 2, ndim
                    affinity(isam, idim : ndim) = affinity(isam, idim : ndim) + sample(isam, idim) * tform(idim : ndim, idim)
                end do
#elif           CLU_ENABLED
                affinity(isam, 1) = sample(isam, 1) PLUS(tlate(1))
                affinity(isam, 2 : ndim) = sample(isam, 1) * tform(2 : ndim, 1) PLUS(tlate(2 : ndim))
                do idim = 2, ndim
                    affinity(isam, idim) = affinity(isam, idim) + sample(isam, idim)
                    affinity(isam, idim + 1 : ndim) = affinity(isam, idim + 1 : ndim) + sample(isam, idim) * tform(idim + 1 : ndim, idim)
                end do
#else
#error          "Unrecognized interface."
#endif
            end do
        end if

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  CHECK_LEN_TLATE
#undef  CHECK_VAL_DIM
#undef  PLUS