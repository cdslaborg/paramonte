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
!>  This file contains the implementation details of the routines under the generic interfaces in [pm_arrayRemove](@ref pm_arrayRemove).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the temporary new array for cases where the result is to be returned in the input array.
#if     setRemoved_ENABLED
#if     SK_ENABLED && D0_D0_ENABLED
        character(:,SKG)            , allocatable :: ArrayRemoved
#elif   SK_ENABLED && (D1_D0_ENABLED || D1_D1_ENABLED)
        character(len(array,IK),SKG), allocatable :: ArrayRemoved(:)
#elif   IK_ENABLED && (D1_D0_ENABLED || D1_D1_ENABLED)
        integer(IKG)                , allocatable :: ArrayRemoved(:)
#elif   LK_ENABLED && (D1_D0_ENABLED || D1_D1_ENABLED)
        logical(LKG)                , allocatable :: ArrayRemoved(:)
#elif   CK_ENABLED && (D1_D0_ENABLED || D1_D1_ENABLED)
        complex(CKG)                , allocatable :: ArrayRemoved(:)
#elif   RK_ENABLED && (D1_D0_ENABLED || D1_D1_ENABLED)
        real(RKG)                   , allocatable :: ArrayRemoved(:)
#else
#error  "Unrecognized interface."
#endif
#elif   !getRemoved_ENABLED
#error  "Unrecognized interface."
#endif
        ! Define logical vs. normal equivalence operators. This becomes relevant only when user-specified comparison function iseq() is missing.
#if     LK_ENABLED
#define IS_EQUAL .eqv.
#elif   SK_ENABLED || IK_ENABLED || CK_ENABLED || RK_ENABLED
#define IS_EQUAL ==
#else
#error  "Unrecognized interface."
#endif
        ! Determine assumed-length scalar character vs. array input arguments.
#if     D0_D0_ENABLED
        integer(IK) :: lenPattern
#define GET_INDEX(i) i : i + lenPattern - 1_IK
#define GET_SIZE len
#if     CusCom_ENABLED
#define ISEQ(segment,pattern) iseq(segment,pattern)
#else
#define ISEQ(segment,pattern) segment == pattern
#endif
#elif   D1_D1_ENABLED
        integer(IK) :: lenPattern
#define GET_INDEX(i) i : i + lenPattern - 1_IK
#define GET_SIZE size
#if     CusCom_ENABLED
#define ISEQ(Segment,pattern) iseq(Segment, pattern, lenPattern)
#else
#define ISEQ(Segment,pattern) all(Segment IS_EQUAL pattern)
#endif
#elif   D1_D0_ENABLED
        integer(IK), parameter :: lenPattern = 1_IK
#define GET_INDEX(i) i
#define GET_SIZE size
#if     CusCom_ENABLED
#define ISEQ(segment,pattern) iseq(segment, pattern)
#elif   DefCom_ENABLED
#define ISEQ(segment,pattern) segment IS_EQUAL pattern
#else
#error  "Unrecognized interface."
#endif
#else
#error  "Unrecognized interface."
#endif
        ! Set the array offset.
#if     D0_D0_ENABLED || getRemoved_ENABLED
        integer(IK) , parameter     :: offset = 0_IK
#else
        integer(IK)                 :: offset
#endif
        ! This `lenArrayOld` serves as the array index offset, to be also used later.
        integer(IK) , allocatable   :: DOP(:) ! pattern Occurrence Position in the array.
        integer(IK)                 :: lenArray, i, iLast
        integer(IK)                 :: lenDOP, lenDOPMax, tokenStart
        integer(IK)                 :: lenArrayOld, lenArrayRemoved, lenArrayCurrent
#if     CusIns_ENABLED
        integer(IK)                 :: lenInstance, lenInstanceNew, maxInstance
        integer(IK) , allocatable   :: InstanceNew(:)
        logical(LK)                 :: sorted_def
        logical(LK)                 :: unique_def
        lenInstance = size(instance, kind = IK)
        if (lenInstance == 0_IK) then
#if         getRemoved_ENABLED
            ArrayRemoved = array
#endif
            return
        end if
#endif
        ! Set the non-default array offset.
#if     !(D0_D0_ENABLED || getRemoved_ENABLED)
        offset = lbound(array,1,IK) - 1_IK
#endif
        ! Set the pattern length.
#if     D0_D0_ENABLED || D1_D1_ENABLED
        lenPattern = GET_SIZE(pattern, kind = IK) ! \warning GET_SIZE is a preprocessor macro.
        if (lenPattern == 0_IK) then
#if         getRemoved_ENABLED
            ArrayRemoved = array
#endif
            return
        end if
#endif
        lenArray = GET_SIZE(array, kind = IK) ! \warning GET_SIZE is a preprocessor macro.
        if (lenArray > lenPattern) then
            lenDOPMax = lenArray / lenPattern + 1_IK
#if         CusIns_ENABLED
            maxInstance = maxval(instance)
            if (minval(instance) > 0_IK .and. maxInstance < lenDOPMax) lenDOPMax = maxInstance
#elif       !DefIns_ENABLED
#error      "Unrecognized interface."
#endif
            ! Find all requested instances of pattern.
            allocate(DOP(lenDOPMax))
            lenDOP = 0_IK
            i = offset + 1_IK
            iLast = offset + lenArray - lenPattern + 1_IK
            loopFindDOP: do
                if (ISEQ(array(GET_INDEX(i)), pattern)) then ! fpp
                    !                if (& ! LCOV_EXCL_LINE
                    !#if             setRemovedDefComDefIns_D1_D0_ENABLED || getRemovedDefComDefIns_D1_D0_ENABLED
                    !#if             CusCom_ENABLED
                    !                iseq(array(i), pattern) & ! \warning ALL is a preprocessor macro. ! LCOV_EXCL_LINE
                    !#else
                    !                array(i) IS_EQUAL pattern & ! \warning ALL is a preprocessor macro. ! LCOV_EXCL_LINE
                    !#endif
                    !#elif           setRemovedDefComDefIns_D1_D1_ENABLED || getRemovedDefComDefIns_D1_D1_ENABLED
                    !#if             CusCom_ENABLED
                    !                iseq(array(i : i + lenPattern - 1), pattern, lenPattern) & ! \warning ALL is a preprocessor macro. ! LCOV_EXCL_LINE
                    !#else
                    !                ALL (array(i : i + lenPattern - 1) IS_EQUAL pattern) & ! \warning ALL is a preprocessor macro. ! LCOV_EXCL_LINE
                    !#endif
                    !#endif
                    !               ) then
                    lenDOP = lenDOP + 1_IK
                    DOP(lenDOP) = i
                    i = i + lenPattern
                    !if (lenDOP == lenDOPMax) exit loopFindDOP
                else
                    i = i + 1_IK
                end if
                if (i > iLast) exit loopFindDOP
            end do loopFindDOP
            ! Remove array at all requested instances of pattern.
            blockInstanceExists: if (lenDOP > 0_IK) then
#if             CusIns_ENABLED
                ! Convert all negative and positive instances to counts from the beginning within the possible range [1, lenDOP].
                !lenInstance = size(instance, kind = IK) ! this is now moved up to quit if zero-length instance is encountered.
                allocate(InstanceNew(lenInstance))
                lenInstanceNew = 0_IK
                i = 0_IK
                ! This loop requires lenInstance to be at least 1, which is guaranteed by the condition after `lenInstance` definition in the above.
                do
                    i = i + 1_IK
                    if (instance(i) > 0_IK .and. instance(i) <= lenDOP) then
                        lenInstanceNew = lenInstanceNew  + 1_IK
                        InstanceNew(lenInstanceNew) = instance(i)
                    elseif (instance(i) < 0_IK .and. instance(i) + lenDOP + 1_IK > 0_IK) then
                        lenInstanceNew = lenInstanceNew  + 1_IK
                        InstanceNew(lenInstanceNew) = instance(i) + lenDOP + 1_IK
                    end if
                    if (i == lenInstance) exit
                end do
                sorted_def = .false._LK
                if (present(sorted)) sorted_def = sorted
                if (.not. sorted_def) call setSorted(InstanceNew(1:lenInstanceNew))
                unique_def = .false._LK
                if (present(unique)) unique_def = unique
                if (unique_def) then
                    lenDOP = lenInstanceNew
                else
                    InstanceNew = getUnique(InstanceNew(1:lenInstanceNew))
                    lenDOP = size(InstanceNew, kind = IK)
                end if
                if (lenDOP == 0_IK) then ! instance is empty, return the input array, untouched.
#if                 getRemoved_ENABLED
                    ArrayRemoved = array
#endif
                    ! The following deallocations are essential since gfortran,
                    ! as of version 10.3, cannot automatically deallocate array upon return.
#if                 CusIns_ENABLED
                    deallocate(InstanceNew)
#endif
                    deallocate(DOP)
                    return
                end if
#define         INSTANCENEW(i) InstanceNew(i)
#else
                !CusIns_ENABLED
#define         INSTANCENEW(i) i
#endif
                !CusIns_ENABLED
                lenArrayRemoved = lenArray - lenDOP * lenPattern
#if             SK_ENABLED && D0_D0_ENABLED
                allocate(character(lenArrayRemoved,SKG) :: ArrayRemoved)
#elif           SK_ENABLED && (D1_D0_ENABLED || D1_D1_ENABLED) && getRemoved_ENABLED
                ! \bug
                ! An Intel ifort compiler bug as of version 2021.4 prevents
                ! the merging of the following allocation with the one after.
                allocate(character(len(array),SKG) :: ArrayRemoved(lenArrayRemoved))
#else
                allocate(ArrayRemoved(offset + 1_IK : offset + lenArrayRemoved))
#endif
                tokenStart = offset + 1_IK
                lenArrayOld = offset
                do i = 1, lenDOP
                    lenArrayCurrent = lenArrayOld + DOP(INSTANCENEW(i)) - tokenStart
                    ArrayRemoved(lenArrayOld+1:lenArrayCurrent) = array(tokenStart : DOP(INSTANCENEW(i)) - 1)
                    tokenStart = DOP(INSTANCENEW(i)) + lenPattern
                    lenArrayOld = lenArrayCurrent
                end do
                ArrayRemoved(lenArrayOld + 1_IK : offset + lenArrayRemoved) = array(tokenStart : offset + lenArray)
#if             CusIns_ENABLED
                ! This is essential since gfortran, as of version 10.3,
                ! cannot automatically deallocate array upon return.
                deallocate(InstanceNew)
#endif
            else blockInstanceExists
#if             getRemoved_ENABLED
                ArrayRemoved = array
#endif
                deallocate(DOP)
                return
            end if blockInstanceExists
            deallocate(DOP)
        elseif (lenArray == lenPattern) then
            if (ISEQ(array(GET_INDEX(offset + 1_IK)), pattern) & ! LCOV_EXCL_LINE
#if         CusIns_ENABLED
            .and. any(abs(instance) == 1_IK) & ! LCOV_EXCL_LINE
#endif
           ) then
#if             SK_ENABLED && D0_D0_ENABLED
                allocate(character(0,SKG) :: ArrayRemoved)
#elif           SK_ENABLED && (D1_D0_ENABLED || D1_D1_ENABLED)
                allocate(character(len(array),SKG) :: ArrayRemoved(0))
#else
                allocate(ArrayRemoved(0))
#endif
            else
                ArrayRemoved = array
            end if
        else ! array is smaller than pattern
#if         getRemoved_ENABLED
            ArrayRemoved = array
#endif
            return
        end if
#if     setRemoved_ENABLED
        call move_alloc(from = ArrayRemoved, to = array)
#endif
#undef  INSTANCENEW
#undef  GET_INDEX
#undef  GET_SIZE
#undef  IS_EQUAL
#undef  ISEQ