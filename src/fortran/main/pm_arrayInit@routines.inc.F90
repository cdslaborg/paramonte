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
!>  This file contains the implementation details of the routines under the generic interface of [pm_arrayInit](@ref pm_arrayInit).
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%
#if     getCoreHalo_ENABLED
        !%%%%%%%%%%%%%%%%%%

        call setCoreHalo(array, core, halo, coffset)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setCoreHalo_ENABLED && (D0_ENABLED || D1_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i
        ! Set the size inference macros for scalar string vs. vectors of arbitrary types.
#if     D0_ENABLED
#define GET_SIZE(array) len(array, IK)
#elif   D1_ENABLED
#define GET_SIZE(array) size(array, 1, IK)
#else
#error  "Unrecognized interface."
#endif
        ! Set the core bounds `csize` for array vs. scalar `core` argument.
#if     Arr_ENABLED
#define GET_CORE(i) core(i : i)
#define GET_CSIZE GET_SIZE(core)
        CHECK_ASSERTION(__LINE__, coffset + GET_SIZE(core) <= GET_SIZE(array), \
        SK_"@setCoreHalo(): The condition `coffset + size(core) <= size(array)` must hold. rank(array), coffset, size(core), size(array) = "\
        //getStr([int(rank(array), IK), coffset, GET_SIZE(core), GET_SIZE(array)])) ! fpp
#elif   Sca_ENABLED
#define GET_CORE(i) core
#define GET_CSIZE csize
        CHECK_ASSERTION(__LINE__, 0_IK <= csize, \
        SK_"@setCoreHalo(): The condition `0 <= csize` must hold. csize = "//getStr(csize)) ! fpp
        CHECK_ASSERTION(__LINE__, coffset + csize <= GET_SIZE(array), \
        SK_"@setCoreHalo(): The condition `coffset + csize <= size(array)` must hold. rank(array), coffset, csize, size(array) = "\
        //getStr([int(rank(array), IK), coffset, csize, GET_SIZE(array)])) ! fpp
#else
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, 0_IK <= coffset, \
        SK_"@setCoreHalo(): The condition `0 <= coffset` must hold. coffset = "//getStr(coffset)) ! fpp
        do concurrent(i = 1_IK : coffset)
            array(i:i) = halo
        end do
        do concurrent(i = 1_IK : GET_CSIZE)
            array(coffset + i : coffset + i) = GET_CORE(i)
        end do
        do concurrent(i = coffset + GET_CSIZE + 1_IK : GET_SIZE(array))
            array(i:i) = halo
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setCoreHalo_ENABLED && (D2_ENABLED || D3_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Set the core bounds `csize` for array vs. scalar `core` argument.
#if     Arr_ENABLED
#define GET_CSIZE(idim) size(core, idim, IK)
#elif   Sca_ENABLED
#define GET_CSIZE(idim) csize(idim)
#else
#error  "Unrecognized interface."
#endif
        ! Set the rank subset notation.
#if     D2_ENABLED
#define SUBSET :
#define GET_SUBSET(X) X(1)
#elif   D3_ENABLED
#define SUBSET :, :
#define GET_SUBSET(X) X(1 : rankArrM1)
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: i, rankArray, rankArrM1
        rankArray = int(rank(array), IK)
        rankArrM1 = rankArray - 1_IK
        array(SUBSET, 1 : coffset(rankArray)) = halo
        do i = 1_IK, GET_CSIZE(rankArray)
            call setCoreHalo( halo = halo & ! LCOV_EXCL_LINE
                            , array = array(SUBSET, coffset(rankArray) + i) & ! LCOV_EXCL_LINE
                            , coffset = GET_SUBSET(coffset) & ! LCOV_EXCL_LINE
#if                         Arr_ENABLED
                            , core = core(SUBSET, i) & ! LCOV_EXCL_LINE
#elif                       Sca_ENABLED
                            , csize = GET_SUBSET(csize) & ! LCOV_EXCL_LINE
                            , core = core & ! LCOV_EXCL_LINE
#else
#error                      "Unrecognized interface."
#endif
                            )
        end do
        array(SUBSET, coffset(rankArray) + GET_CSIZE(rankArray) + 1_IK : size(array, rankArray, IK)) = halo
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif

#undef  GET_SUBSET
#undef  GET_CSIZE
#undef  GET_SIZE
#undef  GET_CORE
#undef  SUBSET