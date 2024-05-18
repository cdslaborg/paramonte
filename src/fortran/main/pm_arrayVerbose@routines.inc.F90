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
!>  This include file contains procedure implementation of [pm_arrayVerbose](@ref pm_arrayVerbose).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Saturday 1:30 AM, August 20, 2016, Institute for Computational Engineering and Sciences, UT Austin, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CHECK_ENABLED
#define CHECK_SUM_WEIGHT(LINE) \
CHECK_ASSERTION(LINE, weisum == sum(weight, mask = weight > 0_IK), \
SK_"The condition `weisum == sum(weight, mask = weight > 0_IK)` must hold: "//\
getStr([weisum, sum(weight, mask = weight > 0_IK)]))
#else
#define CHECK_SUM_WEIGHT(LINE)
#endif
        ! Define indexing style.
#if     D0_ENABLED && SK_ENABLED
#define GET_INDEX(i) i:i
#define GET_SIZE len
#elif   D1_ENABLED || D2_ENABLED
#define GET_INDEX(i) i
#define GET_SIZE size
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%
#if     getVerbose_ENABLED
        !%%%%%%%%%%%%%%%%%

        integer(IK) :: ipnt, iweight, counter
#if     D0_ENABLED || D1_ENABLED
        CHECK_SUM_WEIGHT(__LINE__) ! fpp
        CHECK_ASSERTION(__LINE__, size(weight) == GET_SIZE(array), \
        SK_"@getVerbose(): The size of `weight` must equal the size of `array`. size(array), size(weight) = "//\
        getStr([GET_SIZE(array), size(weight)])) ! fpp
        counter = 0_IK
        do ipnt = 1_IK, GET_SIZE(array, kind = IK)
            do iweight = 1_IK, weight(ipnt)
                counter = counter + 1_IK
                verbose(GET_INDEX(counter)) = array(GET_INDEX(ipnt))
            end do
        end do
#elif   D2_ENABLED
        integer(IK) :: ndim, npnt
        CHECK_SUM_WEIGHT(__LINE__) ! fpp
        CHECK_ASSERTION(__LINE__, dim == 1_IK .or. dim == 2_IK, SK_"@getVerbose(): The input `dim` must be either 1 or 2. dim = "//getStr(dim)) ! fpp
        CHECK_ASSERTION(__LINE__, size(weight) == size(array, dim), SK_"@getVerbose(): The size of `weight` must equal the size of `array` along dimension `dim`. dim, size(array, dim), size(weight) = "//\
        getStr([dim, size(array, dim, IK), size(weight, 1, IK)])) ! fpp
        ndim = size(array, 3 - dim, IK)
        npnt = size(array, dim, IK)
        if (dim == 2_IK) then
            counter = 0_IK
            do ipnt = 1_IK, npnt
                do iweight = 1_IK, weight(ipnt)
                    counter = counter + 1_IK
                    verbose(1:ndim, counter) = array(1:ndim,ipnt)
                end do
            end do
        elseif (dim == 1_IK) then
            ! \todo The memory access pattern can be improved by iterating over 1:ndim.
            counter = 0_IK
            do ipnt = 1_IK, npnt
                do iweight = 1_IK, weight(ipnt)
                    counter = counter + 1_IK
                    verbose(counter, 1:ndim) = array(ipnt, 1:ndim)
                end do
            end do
        end if
#endif
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  CHECK_SUM_WEIGHT
#undef  GET_INDEX
#undef  GET_SIZE