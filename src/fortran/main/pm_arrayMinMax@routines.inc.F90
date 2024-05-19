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
!>  This include file contains procedure implementations of [pm_arrayMinMax](@ref pm_arrayMinMax).
!>
!>  \author
!>  \AmirShahmoradi, Sunday 3:33 AM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! DEfine the comparison operation.
#if     LK_ENABLED
#define IS_MORE(a,b) a .and. .not. b
#define IS_LESS(a,b) b .and. .not. a
#elif   CK_ENABLED
#define IS_MORE(a,b) a%re > b%re .or. (a%re == b%re .and. a%im > b%im)
#define IS_LESS(a,b) a%re < b%re .or. (a%re == b%re .and. a%im < b%im)
#elif   BSSK_ENABLED || PSSK_ENABLED
#define IS_MORE(a,b) a%val > b%val
#define IS_LESS(a,b) a%val < b%val
#else
#define IS_MORE(a,b) a > b
#define IS_LESS(a,b) a < b
#endif
        !%%%%%%%%%%%%%%%%%%%
#if     getMinMaxVal_ENABLED
        !%%%%%%%%%%%%%%%%%%%

#if     SK_ENABLED && D0_ENABLED
        call setMinMaxVal(array, minMaxVal(1:1), minMaxVal(2:2))
#else
        call setMinMaxVal(array, minMaxVal(1), minMaxVal(2))
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMinMaxVal_ENABLED && D0_ENABLED && SK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: iell
        character(1,SKG), allocatable :: empty(:)
        allocate(empty(0))
        vmin = minval(empty)
        vmax = maxval(empty)
        do iell = 1, len(array, IK)
            if (IS_MORE(vmin, array(iell : iell))) vmin = array(iell : iell)
            if (IS_LESS(vmax, array(iell : iell))) vmax = array(iell : iell)
        end do
        deallocate(empty)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMinMaxVal_ENABLED && D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: iell
#if     BSSK_ENABLED || PSSK_ENABLED
        if (size(array, 1, IK) == 0_IK) return
        vmin = array(1)
        vmax = array(1)
#elif   SK_ENABLED
        character(len(vmin, IK),SKG), allocatable :: empty(:)
        if (allocated(empty)) deallocate(empty)
        allocate(empty(0))
        vmin = minval(empty)
        vmax = maxval(empty)
#elif   IK_ENABLED
        vmin = +huge(0_IKG)
        vmax = -huge(0_IKG)
#elif   LK_ENABLED
        vmin = .true._LKG
        vmax = .false._LKG
#elif   CK_ENABLED
        complex(CKG), parameter :: POSINF = cmplx(+huge(0._CKG), +huge(0._CKG), CKG)
        complex(CKG), parameter :: NEGINF = cmplx(-huge(0._CKG), -huge(0._CKG), CKG)
        vmin = POSINF
        vmax = NEGINF
#elif   RK_ENABLED
        vmin = +huge(0._RKG)
        vmax = -huge(0._RKG)
#else
#error  "Unrecognized interface."
#endif
        do iell = 1, size(array, 1, IK)
            if (IS_MORE(vmin, array(iell))) vmin = array(iell)
            if (IS_LESS(vmax, array(iell))) vmax = array(iell)
        end do

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  IS_MORE
#undef  IS_LESS