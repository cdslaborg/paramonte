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
!>  This file contains implementations of procedures [pm_knn](@ref pm_knn).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Thursday 8:40 PM, July 20, 2023, Dallas, TX
!>  \AmirShahmoradi, Saturday 1:00 AM, September, 1, 2018, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%
#if     setKnnSorted_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ipnt, nref, npnt
#if     Kth_ENABLED
        real(TKG) :: selection
        CHECK_ASSERTION(__LINE__, 0_IK < k .and. k <= size(distance, 1, IK), SK_"@setKnnSorted(): The condition `0 < k .and. k <= size(distance, 1)` must hold. k, shape(distance) = "//getStr([k, shape(distance, IK)]))
#elif   Ind_ENABLED
        CHECK_ASSERTION(__LINE__, all(shape(distance) == shape(rank)), SK_"@setKnnSorted(): The condition `all(shape(distance, IK) == shape(rank, IK))` must hold. shape(distance), shape(rank) = "//getStr([shape(distance), shape(rank)]))
#endif
        nref = size(distance, 1, IK)
        npnt = size(distance, 2, IK)
        DO_CONCURRENT(ipnt,1,npnt)
#if         Val_ENABLED
            call setSorted(distance(1:nref, ipnt))
#elif       Ind_ENABLED
            call setSorted(distance(1:nref, ipnt), rank(1:nref, ipnt))
#elif       Kth_ENABLED
            call setSelected(selection, distance(1:nref, ipnt), k)
#else
#error      "Unrecognized interface."
#endif
        end do

        !%%%%%%%%%%%%%%%%
#elif   setKnnVal_ENABLED
        !%%%%%%%%%%%%%%%%

        integer(IK) :: ipnt, nref, npnt
        nref = size(distance, 1, IK)
        npnt = size(distance, 2, IK)
        CHECK_ASSERTION(__LINE__, size(val, 1, IK) == size(distance, 2, IK), SK_"@setKnnVal(): The condition `size(val) == size(distance, 2)` must hold. size(val), shape(distance) = "//getStr([shape(val, IK), shape(distance, IK)]))
        DO_CONCURRENT(ipnt,1,npnt)
            call setSelected(selection, distance(1:nref, ipnt), k)
        end do

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
