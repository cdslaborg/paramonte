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
!>  This include file contains procedure implementations of [pm_arrayRefine](@ref pm_arrayRefine).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Saturday 1:48 AM, August 20, 2016, Institute for Computational Engineering and Sciences, UT Austin, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define array indexing rule.
#if     SK_ENABLED && D0_ENABLED
#define GET_INDEX(i) i:i
#define GET_SIZE(X)len(X, IK)
#elif   D1_ENABLED || D2_ENABLED
#define GET_INDEX(i) i
#define GET_SIZE(X)size(X, 1, IK)
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%
#if     getRefined_ENABLED
        !%%%%%%%%%%%%%%%%%

        integer(IK) :: rsize, weisum, weightRefined(size(weight, 1, IK))
        weightRefined = weight
        arrayRefined = array
#if     D0_ENABLED || D1_ENABLED
        call setRefined(arrayRefined, weightRefined, skip, rsize)
        if (0_IK < rsize) then
            weisum = sum(weightRefined(1 : rsize), mask = weightRefined(1 : rsize) > 0_IK)
            arrayRefined = getVerbose(arrayRefined(1 : rsize), weightRefined(1 : rsize), weisum)
        else
            call setResized(arrayRefined, 0_IK)
        end if
#elif   D2_ENABLED
        call setRefined(arrayRefined, dim, weightRefined, skip, rsize)
        if (0_IK < rsize) then
            weisum = sum(weightRefined(1 : rsize), mask = weightRefined(1 : rsize) > 0_IK)
            if (dim == 1_IK) then
                arrayRefined = getVerbose(arrayRefined(1 : rsize, :), weightRefined(1 : rsize), weisum, dim)
            else
                arrayRefined = getVerbose(arrayRefined(:, 1 : rsize), weightRefined(1 : rsize), weisum, dim)
            end if
        else
            if (dim == 1_IK) then
                call setResized(arrayRefined, [0_IK, size(array, 2, IK)])
            else
                call setResized(arrayRefined, [size(array, 1, IK), 0_IK])
            end if
        end if
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%
#elif   setRefined_ENABLED
        !%%%%%%%%%%%%%%%%%

        integer(IK) :: isam
        if (GET_SIZE(array) == 0_IK) then
            rsize = 0_IK
            return
        end if
        rsize = 0_IK
        call setReweight(weight, skip)
#if     D0_ENABLED || D1_ENABLED
        CHECK_ASSERTION(__LINE__, size(weight, 1, IK) == GET_SIZE(array), SK_": The condition `size(weight) == size/len(array)` must hold. size(weight), size/len(array) = "//getStr([size(weight, 1, IK), GET_SIZE(array)]))
        do isam = 1, GET_SIZE(array)
            if (0_IK < weight(isam)) then
                rsize = rsize + 1_IK
                if (rsize < isam) then ! The only other possibility is equality.
                    weight(rsize) = weight(isam)
                    array(GET_INDEX(rsize)) = array(GET_INDEX(isam))
                end if
            end if
        end do
#elif   D2_ENABLED
        CHECK_ASSERTION(__LINE__, dim == 1 .or. dim == 2, SK_": The condition `dim == 1 .or. dim == 2` must hold. dim = "//getStr(dim))
        CHECK_ASSERTION(__LINE__, size(weight, 1, IK) == size(array, dim, IK), SK_": The condition `size(weight) == size(array, dim)` must hold. size(weight), size(array, dim) = "//getStr([size(weight, 1, IK), size(array, dim, IK)]))
        if (dim == 1_IK) then
            do isam = 1, size(array, dim, IK)
                if (0_IK < weight(isam)) then
                    rsize = rsize + 1_IK
                    if (rsize < isam) then ! The only other possibility is equality.
                        weight(rsize) = weight(isam)
                        array(rsize, :) = array(isam, :)
                    end if
                end if
            end do
        else
            do isam = 1, size(array, dim, IK)
                if (0_IK < weight(isam)) then
                    rsize = rsize + 1_IK
                    if (rsize < isam) then ! The only other possibility is equality.
                        weight(rsize) = weight(isam)
                        array(:, rsize) = array(:, isam)
                    end if
                end if
            end do
        end if
#else
#error  "Unrecognized interface."
#endif
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  GET_INDEX
#undef  GET_SIZE
