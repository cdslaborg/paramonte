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
!>  This file contains the implementation details of the routines of [pm_arrayRemap](@ref pm_arrayRemap).
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK_ENABLED && D0_ENABLED
!#define ALLOCATE_ARRAYNEW allocate(character(lenArray)::arrayNew)
#define GET_LBOUND(array) 1
#define GET_INDEX(i) i:i
#define GET_SIZE(X)len(X, IK)
#else
!#define ALLOCATE_ARRAYNEW allocate(arrayNew, mold = array)
#define GET_INDEX(i) i
#define GET_SIZE(X)size(X, 1, IK)
#endif
        ! Declare the temporary array.
#if     Old_ENABLED
#if     SK_ENABLED && D0_ENABLED
        character(:,SKG), allocatable :: arrayNew
#elif   SK_ENABLED && D1_ENABLED
        character(len(array, IK),SKG), allocatable :: arrayNew(:)
#elif   IK_ENABLED && D1_ENABLED
        integer(IKG), allocatable :: arrayNew(:)
#elif   LK_ENABLED && D1_ENABLED
        logical(LKG), allocatable :: arrayNew(:)
#elif   CK_ENABLED && D1_ENABLED
        complex(CKG), allocatable :: arrayNew(:)
#elif   RK_ENABLED && D1_ENABLED
        real(RKG), allocatable :: arrayNew(:)
#else
#error  "Unrecognized interface."
#endif
#elif   !New_ENABLED
#error  "Unrecognized interface."
#endif
        integer(IK) :: i, lenArray
#if     New_ENABLED || (SK_ENABLED && D0_ENABLED)
        integer(IK), parameter :: offset = 0_IK
#elif   Old_ENABLED
        integer(IK) :: offset
        offset = lbound(array, 1, IK) - 1_IK
#else
#error  "Unrecognized interface."
#endif
        lenArray = GET_SIZE(array)
#if     setRemapped_ENABLED && New_ENABLED
        CHECK_ASSERTION(__LINE__, lenArray == GET_SIZE(arrayNew), SK_"@setRemapped(): The lengths of the arguments `array` and `arrayNew` must equal. lenArray, lenArrayNew = "//getStr([lenArray, GET_SIZE(arrayNew)]))
#endif
        CHECK_ASSERTION(__LINE__, all(1_IK + offset <= index) .and. all(index <= lenArray + offset), SK_"@setRemapped(): All `index` values must be within the lower and upper bounds of the input `array`. index, lb, ub = "//getStr([index, 1_IK + offset, lenArray + offset]))
        CHECK_ASSERTION(__LINE__, lenArray == size(index, 1, IK), SK_"@setRemapped(): The size of the arguments `array` and `index` must equal. lenArray, lenIndex = "//getStr([lenArray, size(index, 1, IK)]))
#if     Old_ENABLED
        allocate(arrayNew, mold = array)
#endif
        do i = 1_IK, lenArray
#if             Rev_ENABLED
                arrayNew(GET_INDEX(i + offset)) = array(GET_INDEX(index(lenArray)))
                lenArray = lenArray - 1_IK
#elif           For_ENABLED
                arrayNew(GET_INDEX(i + offset)) = array(GET_INDEX(index(i)))
#else
#error          "Unrecognized interface."
#endif
        end do
#if     Old_ENABLED
        call move_alloc(from = arrayNew, to = array)
#if     __GFORTRAN__
        !>  \todo
        !>  The following bug bypass must be resolved once the Gfortran bug is fixed.
        !>  \bug gfortran 10.3 does not deallocate `arrayNew` upon return.
        if (allocated(arrayNew)) deallocate(arrayNew)
#endif
#endif
#undef  GET_INDEX
#undef  GET_SIZE