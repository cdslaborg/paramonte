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
!>  This file contains the implementation details of the routines under the generic interface
!>  [setReversed](@ref pm_arrayReverse::setReversed).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Set the indexing rules.
#if     SK_ENABLED && D0_ENABLED
#define GET_INDEX(i) i:i
#define GET_SIZE len
#elif   D1_ENABLED
#define GET_SIZE size
#define GET_INDEX(i) i
#else
#error  "Unrecognized interface."
#endif
        ! Set the array bounds.
#if     SK_ENABLED && D0_ENABLED && Old_ENABLED
#define GEN_UBOUND(array) len(array, kind = IK)
#define GET_LBOUND(array) 1_IK
#else
#define GEN_UBOUND(array) ubound(array, 1, kind = IK)
#define GET_LBOUND(array) lbound(array, 1, kind = IK)
#endif
        ! Declare the temporary storage.
#if     Old_ENABLED
#if     SK_ENABLED && D0_ENABLED
        character(1,SKG) :: temp
#elif   SK_ENABLED && D1_ENABLED
        character(len(array),SKG) :: temp
#elif   IK_ENABLED && D1_ENABLED
        integer(IKG) :: temp
#elif   LK_ENABLED && D1_ENABLED
        logical(LKG) :: temp
#elif   CK_ENABLED && D1_ENABLED
        complex(CKG) :: temp
#elif   RK_ENABLED && D1_ENABLED
        real(RKG) :: temp
#elif   PSSK_ENABLED && D1_ENABLED
        use pm_container, only: css_pdt
        type(css_pdt(SKG)) :: temp
#else
#error  "Unrecognized interface."
#endif
#elif   !New_ENABLED
#error  "Unrecognized interface."
#endif
        ! Begin the reversal.
#if     Old_ENABLED
        integer(IK) :: lb, ub
        lb = GET_LBOUND(array)
        ub = GEN_UBOUND(array)
        do
            if (lb >= ub) exit
            temp = array(GET_INDEX(lb))
            array(GET_INDEX(lb)) = array(GET_INDEX(ub))
            array(GET_INDEX(ub)) = temp
            ub = ub - 1_IK
            lb = lb + 1_IK
        end do
#elif   New_ENABLED
        integer(IK) :: i, lenArray
        lenArray = GET_SIZE(array, kind = IK)
#if     setReversed_ENABLED
        CHECK_ASSERTION(__LINE__, GET_SIZE(array, kind = IK) == GET_SIZE(ArrayReversed, kind = IK), \
        SK_"@setReversed(): The sizes of the input `array` and `ArrayReversed` must equal. lenArray, lenArrayNew = "\
        //getStr([GET_SIZE(array, kind = IK), GET_SIZE(ArrayReversed, kind = IK)]))
#endif
        do i = 1_IK, lenArray
            ArrayReversed(GET_INDEX(i)) = array(GET_INDEX(lenArray - i + 1_IK))
        end do
#else
#error  "Unrecognized interface."
#endif
#undef  GEN_UBOUND
#undef  GET_LBOUND
#undef  GET_INDEX
#undef  GET_SIZE