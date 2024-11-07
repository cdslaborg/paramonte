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
!>  This file contains the implementation details of the routines under the generic interfaces in [pm_arraySplit](@ref pm_arraySplit).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define vector vs. scalar array + sep and operations.
#if     D0_D0_ENABLED
#define GET_SIZE len
#define ALL(item) item
#define GET_INDEX(i) i : i + sepLen - 1_IK
#elif   D1_D0_ENABLED
#define GET_SIZE size
#define ALL(item) all(item)
#define GET_INDEX(i) i : i + sepLen - 1_IK
#elif   D1_D1_ENABLED
#define GET_SIZE size
#define ALL(item) all(item)
#define GET_INDEX(i) i : i + sepLen - 1_IK
#else
#error  "Unrecognized interface."
#endif
        ! Define logical vs. normal equivalence operators.
#if     LK_ENABLED
#define IS_EQUAL .eqv.
#elif   SK_ENABLED || IK_ENABLED || CK_ENABLED || RK_ENABLED
#define IS_EQUAL ==
#else
#error  "Unrecognized interface."
#endif

        ! Define custom-comparison macro.
#if     CusCom_ENABLED && D1_D1_ENABLED
#define ISEQ(Segment,sep) iseq(Segment, sep, sepLen)
#elif   CusCom_ENABLED && (D0_D0_ENABLED || D1_D0_ENABLED)
#define ISEQ(segment,sep) iseq(segment, sep)
#elif   DefCom_ENABLED
#define ISEQ(segment,sep) ALL(segment IS_EQUAL sep)
#else
#error  "Unrecognized interface."
#endif

        ! Define the scalar vs. array interface for container component assignment.
#if     Jagged_ENABLED && D1_D0_ENABLED
#define GET_ARRAY(object) [object]
#elif   Jagged_ENABLED && (D0_D0_ENABLED || D1_D1_ENABLED)
#define GET_ARRAY(object) object
#elif   !(Index_ENABLED || Fixed_ENABLED)
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     setSplit_ENABLED && Fixed_ENABLED && DefIns_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the return value when no splitting occurs.

        logical(LK) :: keep_def
        integer(IK) :: lenArray, i, iLast, offset
#if     D1_D0_ENABLED
        integer(IK), parameter :: sepLen = 1_IK
#elif   D0_D0_ENABLED || D1_D1_ENABLED
        integer(IK) :: sepLen
        sepLen = GET_SIZE(sep, kind = IK)
#else
#error  "Unrecognized interface."
#endif
        offset = 1_IK
        if (present(keep)) then
            if (keep) offset = sepLen
        end if
        lenArray = GET_SIZE(array, kind = IK)
        CHECK_ASSERTION(__LINE__, 1_IK < size(field, 1, IK), SK_"@setSplit(): The condition `1 < size(field)` must hold. shape(field) = "//getStr(shape(field, IK)))
#if     D0_D0_ENABLED || D1_D1_ENABLED
        if (sepLen == 0_IK) then
            nsplit = 1_IK
            field(1) = 1_IK
            field(2) = lenArray + 1_IK
            return
        end if
#endif
        if (sepLen < lenArray) then
            i = 1_IK
            iLast = lenArray - sepLen + 1_IK
            if (present(keep)) then
                if (keep) then
                    loopFindSplitStart1: do
                        if (ISEQ(array(GET_INDEX(i)), sep)) then
                            CHECK_ASSERTION(__LINE__, nsplit < size(field, 1, IK), SK_"@setSplit(): The condition `nsplit < size(field, 1)` must hold. shape(field) = "//getStr(shape(field, IK)))
                            nsplit = nsplit + 1_IK
                            field(nsplit) = i
                            i = i + sepLen
                            !if (nsplit == nsplitMax) exit loopFindSplitStart1
                        else
                            i = i + 1_IK
                        end if
                        if (i > iLast) exit loopFindSplitStart1
                    end do loopFindSplitStart1

                    return
                end if
            end if
            loopFindSepPos: do
                if (ISEQ(array(GET_INDEX(i)), sep)) then
                    CHECK_ASSERTION(__LINE__, nsplit < size(field, 1, IK), SK_"@setSplit(): The condition `nsplit < size(field, 1)` must hold. shape(field) = "//getStr(shape(field, IK)))
                    nsplit = nsplit + 1_IK
                    field(nsplit) = i
                    i = i + sepLen
                    !if (nsplit == nsplitMax) exit loopFindSepPos
                else
                    i = i + 1_IK
                end if
                if (i > iLast) exit loopFindSepPos
            end do loopFindSepPos

            ! Split array at all requested instances of sep.

            blockInstanceExists: if (nsplit > 0_IK) then
                if (present(keep)) then
                    keep_def = keep
                else
                    keep_def = .false._LK
                end if
                if (keep_def) then
                    !lenArraySplit = nsplit + 1_IK
                    allocate(field(2, 2 * nsplit + 1))
                    splitCounter = 1_IK
                    field(1, splitCounter) = 1_IK
                    do i = 1_IK, nsplit
                        field(2, splitCounter) = sepLoc(i) - 1_IK
                        field(1, splitCounter + 1) = sepLoc(i)
                        field(2, splitCounter + 1) = sepLoc(i) + sepLen - 1_IK
                        splitCounter = splitCounter + 2_IK
                        field(1, splitCounter) = sepLoc(i) + sepLen
                    end do
                    field(2, splitCounter) = lenArray
                else
                    allocate(field(2, nsplit + 1))
                    field(1, 1) = 1_IK
                    do i = 1_IK, nsplit
                        field(2, i) = sepLoc(i) - 1_IK
                        field(1, i + 1) = sepLoc(i) + sepLen
                    end do
                    field(2, i) = lenArray
                end if
            else blockInstanceExists
                nsplit = 1_IK
                field(1) = 1_IK
                field(2) = lenArray + sepLen + 1_IK
            end if blockInstanceExists
            deallocate(sepLoc)
        elseif (lenArray == sepLen) then
            if (ALL(array IS_EQUAL sep)) then
                if (present(keep)) then
                    if (keep) then
                        nsplit = 3_IK
                        field(2) = 1_IK
                        field(3) = lenArray + 1_IK
                        field(4) = lenArray + 1_IK
                        return
                    end if
                end if
                nsplit = 2_IK
                field(1) = 1_IK
                field(2) = 0_IK
                field(1) = lenArray
                field(2) = lenArray - 1_IK
            else
                RETURN_WHOLE_ARRAY
            end if
        else ! `array` is smaller than `sep`. Return whole `array`.
            nsplit = 1_IK
            field(1) = 1_IK
            field(2) = lenArray + 1_IK
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setSplit_ENABLED && Fixed_ENABLED && CusIns_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setSplit_ENABLED && (Index_ENABLED || Jagged_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the return value when no splitting occurs.
#if     Jagged_ENABLED
        integer(IK)                 :: tokenStart
#define RETURN_WHOLE_ARRAY \
allocate(field(1)); field(1)%val = array
#elif   Index_ENABLED
#define RETURN_WHOLE_ARRAY \
allocate(field(2, 1)); field(1, 1) = 1_IK; field(2, 1) = GET_SIZE(array, kind = IK)
#else
#error  "Unrecognized interface."
#endif

#if     D0_D0_ENABLED || D1_D1_ENABLED
        integer(IK)                 :: sepLen
#elif   D1_D0_ENABLED
        integer(IK) , parameter     :: sepLen = 1_IK
#else
#error  "Unrecognized interface."
#endif

        integer(IK) , allocatable   :: sepLoc(:) ! sep Occurrence Position in the array.
        integer(IK)                 :: lenArray, i, iLast
        integer(IK)                 :: sepLocLen, sepLocLenMax
       !integer(IK)                 :: lenArraySplit
        logical(LK)                 :: keep_def
        integer(IK)                 :: splitCounter

#if     CusIns_ENABLED
        integer(IK)                 :: lenInstance, lenInstanceNew, instanceMax
        integer(IK) , allocatable   :: instanceNew(:)
        logical(LK)                 :: sorted_def
        logical(LK)                 :: unique_def

        lenInstance = size(instance, kind = IK)
        if (lenInstance == 0_IK) then
            RETURN_WHOLE_ARRAY
            return
        end if
#elif   !DefIns_ENABLED
#error  "Unrecognized interface."
#endif

#if     D0_D0_ENABLED || D1_D1_ENABLED
        sepLen = GET_SIZE(sep, kind = IK)
        if (sepLen == 0_IK) then
            RETURN_WHOLE_ARRAY
            return
        end if
#endif
        lenArray = GET_SIZE(array, kind = IK)
        if (lenArray > sepLen) then

            sepLocLenMax = lenArray / sepLen + 1_IK
#if         CusIns_ENABLED
            instanceMax = maxval(instance)
            if (minval(instance) > 0_IK .and. instanceMax < sepLocLenMax) sepLocLenMax = instanceMax
#endif
            ! Find all requested instances of sep.

            allocate(sepLoc(sepLocLenMax))
            i = 1_IK
            sepLocLen = 0_IK
            iLast = lenArray - sepLen + 1_IK
            loopFindSepPos: do
                if (ISEQ(array(GET_INDEX(i)), sep)) then
                    sepLocLen = sepLocLen + 1_IK
                    sepLoc(sepLocLen) = i
                    i = i + sepLen
                    !if (sepLocLen == sepLocLenMax) exit loopFindSepPos
                else
                    i = i + 1_IK
                end if
                if (i > iLast) exit loopFindSepPos
            end do loopFindSepPos

            ! Split array at all requested instances of sep.

            blockInstanceExists: if (sepLocLen > 0_IK) then
#if             CusIns_ENABLED
                ! Convert all negative and positive instances to counts from the beginning within the possible range [1, sepLocLen].
                !lenInstance = size(instance, kind = IK) ! this is now moved up to quit if zero-length instance is encountered.
                allocate(instanceNew(lenInstance))
                lenInstanceNew = 0_IK
                i = 0_IK
                ! This loop requires lenInstance to be at least 1,
                ! which is guaranteed by the condition after `lenInstance` definition in the above.
                do
                    i = i + 1_IK
                    if (instance(i) > 0_IK .and. instance(i) <= sepLocLen) then
                        lenInstanceNew = lenInstanceNew  + 1_IK
                        instanceNew(lenInstanceNew) = instance(i)
                    elseif (instance(i) < 0_IK .and. instance(i) + sepLocLen + 1_IK > 0_IK) then
                        lenInstanceNew = lenInstanceNew  + 1_IK
                        instanceNew(lenInstanceNew) = instance(i) + sepLocLen + 1_IK
                    end if
                    if (i == lenInstance) exit
                end do

                sorted_def = .false._LK
                if (present(sorted)) sorted_def = sorted
                if (.not. sorted_def) call setSorted(instanceNew(1:lenInstanceNew))

                unique_def = .false._LK
                if (present(unique)) unique_def = unique
                if (unique_def) then
                    sepLocLen = lenInstanceNew
                else
                    instanceNew = getUnique(instanceNew(1:lenInstanceNew))
                    sepLocLen = size(instanceNew, kind = IK)
                end if

                if (sepLocLen == 0_IK) then ! instance is empty, return the input array, untouched.
                    RETURN_WHOLE_ARRAY
                    ! The following deallocations are essential since gfortran, as of version 10.3, cannot automatically deallocate array upon return.
                    deallocate(instanceNew)
                    deallocate(sepLoc)
                    return
                end if
#define         INSTANCENEW(i) instanceNew(i)
#else
#define         INSTANCENEW(i) i
#endif
                if (present(keep)) then
                    keep_def = keep
                else
                    keep_def = .false._LK
                end if
                if (keep_def) then
                    !lenArraySplit = sepLocLen + 1_IK
#if                 Jagged_ENABLED
                    allocate(field(2_IK * sepLocLen + 1_IK))
                    tokenStart = 1_IK
                    splitCounter = 1_IK
                    do i = 1_IK, sepLocLen
                        field(splitCounter)%val = array(tokenStart : sepLoc(INSTANCENEW(i)) - 1_IK)
                        field(splitCounter + 1)%val = GET_ARRAY(sep)
                        tokenStart = sepLoc(INSTANCENEW(i)) + sepLen
                        splitCounter = splitCounter + 2_IK
                    end do
                    field(splitCounter)%val = array(tokenStart : lenArray)
#elif               Index_ENABLED
                    allocate(field(2, 2 * sepLocLen + 1))
                    splitCounter = 1_IK
                    field(1, splitCounter) = 1_IK
                    do i = 1_IK, sepLocLen
                        field(2, splitCounter) = sepLoc(INSTANCENEW(i)) - 1_IK
                        field(1, splitCounter + 1) = sepLoc(INSTANCENEW(i))
                        field(2, splitCounter + 1) = sepLoc(INSTANCENEW(i)) + sepLen - 1_IK
                        splitCounter = splitCounter + 2_IK
                        field(1, splitCounter) = sepLoc(INSTANCENEW(i)) + sepLen
                    end do
                    field(2, splitCounter) = lenArray
#else
#error              "Unrecognized interface."
#endif
                else
                    !lenArraySplit = sepLocLen + 1_IK
#if                 Jagged_ENABLED
                    allocate(field(sepLocLen + 1_IK))
                    tokenStart = 1_IK
                    do i = 1, sepLocLen
                        field(i)%val = array(tokenStart : sepLoc(INSTANCENEW(i)) - 1)
                        tokenStart = sepLoc(INSTANCENEW(i)) + sepLen
                    end do
                    field(i)%val = array(tokenStart : lenArray)
#elif               Index_ENABLED
                    allocate(field(2, sepLocLen + 1))
                    field(1, 1) = 1_IK
                    do i = 1_IK, sepLocLen
                        field(2, i) = sepLoc(INSTANCENEW(i)) - 1_IK
                        field(1, i + 1) = sepLoc(INSTANCENEW(i)) + sepLen
                    end do
                    field(2, i) = lenArray
#else
#error              "Unrecognized interface."
#endif
                end if
#if             CusIns_ENABLED
                deallocate(instanceNew) ! This is essential since gfortran, as of version 10.3, cannot automatically deallocate array upon return.
#endif
            else blockInstanceExists
                !lenArraySplit = 1_IK
                RETURN_WHOLE_ARRAY
            end if blockInstanceExists

            deallocate(sepLoc)

        elseif (lenArray == sepLen) then

            if (ALL(array IS_EQUAL sep) & ! LCOV_EXCL_LINE ! \warning ALL is a preprocessor macro.
#if         CusIns_ENABLED
            .and. any(abs(instance) == 1_IK) & ! LCOV_EXCL_LINE
#endif
           ) then
                if (present(keep)) then
                    if (keep) then
                        !lenArraySplit = 2_IK
#if                     Jagged_ENABLED
                        allocate(field(3))
                        field(1)%val = array(1:0)
                        field(2)%val = array ! == [sep]
                        field(3)%val = array(1:0)
#elif                   Index_ENABLED
                        allocate(field(2,3))
                        field(1,1) = 1_IK
                        field(2,1) = 0_IK
                        field(1,2) = 1_IK
                        field(2,2) = lenArray
                        field(1,3) = lenArray
                        field(2,3) = lenArray - 1_IK
#else
#error                  "Unrecognized interface."
#endif
                        return
                    end if
                end if
                !lenArraySplit = 2_IK
#if             Jagged_ENABLED
                allocate(field(2))
                field(1)%val = array(1:0)
                field(2)%val = array(1:0)
#elif           Index_ENABLED
                allocate(field(2,2))
                field(1, 1) = 1_IK
                field(2, 1) = 0_IK
                field(1, 2) = lenArray
                field(2, 2) = lenArray - 1_IK
#else
#error          "Unrecognized interface."
#endif
            else
                !lenArraySplit = 1_IK
                RETURN_WHOLE_ARRAY
            end if

        else ! array is smaller than sep

            !lenArraySplit = 1_IK
            RETURN_WHOLE_ARRAY

        end if

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  RETURN_WHOLE_ARRAY
#undef  INSTANCENEW
#undef  GET_ARRAY
#undef  GET_INDEX
#undef  GET_SIZE
#undef  IS_EQUAL
#undef  ISEQ
#undef  ALL