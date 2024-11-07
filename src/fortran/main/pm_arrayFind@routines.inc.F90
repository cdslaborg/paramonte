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
!>  This file contains the implementation details of the routines under the generic interfaces in [pm_arrayFind](@ref pm_arrayFind).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define logical vs. normal equivalence operators.
#if     LK_ENABLED
#define IS_EQUAL .eqv.
#elif   SK_ENABLED || IK_ENABLED || CK_ENABLED || RK_ENABLED
#define IS_EQUAL ==
#else
#error  "Unrecognized interface."
#endif
        ! Determine assumed-length scalar character vs. array input arguments.
#if     D0_D0_ENABLED && SK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%
#define GET_INDEX(i) i : i + lenPattern - 1_IK
#define GET_SIZE len
#if     DefCom_ENABLED
#define ISEQ(segment,pattern) segment == pattern
#elif   CusCom_ENABLED
#define ISEQ(segment,pattern) iseq(segment, pattern)
#else
#error  "Unrecognized interface."
#endif
#elif   D1_D1_ENABLED
        !%%%%%%%%%%%%
#define GET_INDEX(i) i : i + lenPattern - 1_IK
#define GET_SIZE size
#if     DefCom_ENABLED
#define ISEQ(segment,pattern) all(segment IS_EQUAL pattern)
#elif   CusCom_ENABLED
#define ISEQ(segment,pattern) iseq(segment, pattern, lenPattern)
#else
#error  "Unrecognized interface."
#endif
#elif   D1_D0_ENABLED
        !%%%%%%%%%%%%
#define GET_INDEX(i) i
#define GET_SIZE size
#if     DefCom_ENABLED
#define ISEQ(segment,pattern) segment IS_EQUAL pattern
#elif   CusCom_ENABLED
#define ISEQ(segment,pattern) iseq(segment, pattern)
#else
#error  "Unrecognized interface."
#endif
#else
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
        !%%%%%%%%%%%%%%%%%%
#if     getCountLoc_ENABLED
        !%%%%%%%%%%%%%%%%%%

        integer(IK) :: lenArray, i, targetLen, blindness_def
#if     D1_D0_ENABLED
        integer(IK), parameter :: lenPattern = 1_IK
#elif   D0_D0_ENABLED || D1_D1_ENABLED
        integer(IK) :: lenPattern
#else
#error  "Unrecognized interface."
#endif
#if     DisBor_ENABLED
        logical(LK) :: isDiscrete
        isDiscrete = .true._LK
#endif
        count = 0_IK
#if     D0_D0_ENABLED || D1_D1_ENABLED
        lenPattern = GET_SIZE(pattern, kind = IK)
        if (lenPattern == 0_IK) return
#endif
        lenArray = GET_SIZE(array, kind = IK)
        if (lenPattern < lenArray) then
            ! lenArray shall not be zero in this block.
            if (present(blindness)) then
                CHECK_ASSERTION(__LINE__, blindness /= 0_IK, SK_"@setLoc()/@getLoc(): The input `blindness` must be non-zero. blindness = "//getStr([blindness]))
                blindness_def = blindness
            else
                blindness_def = 1_IK
            end if
            ! Find all requested instances of pattern.
            i = 1_IK
            targetLen = lenArray - lenPattern + 1_IK
            loopFindPattern: do
#if             DefBor_ENABLED
                if (ISEQ(array(GET_INDEX(i)), pattern)) then
                    count = count + 1_IK
                    i = i + blindness_def
                else
                    i = i + 1_IK
                end if
#elif           DisBor_ENABLED
                if (ISEQ(array(GET_INDEX(i)), pattern)) then
                    i = i + blindness_def
                    if (isDiscrete) then
                        isDiscrete = .false._LK
                        count = count + 1_IK
                    end if
                else
                    i = i + 1_IK
                    isDiscrete = .true._LK
                end if
#else
#error          "Unrecognized interface."
#endif
                if (targetLen < i) exit loopFindPattern
            end do loopFindPattern
        elseif (lenArray == lenPattern) then
            if (ISEQ(array(GET_INDEX(1)), pattern)) count = 1_IK
        end if

        !%%%%%%%%%%%%%
#elif   getLoc_ENABLED
        !%%%%%%%%%%%%%

        integer(IK)                 :: lenArray, i, targetLen, blindness_def
        integer(IK)                 :: count, countMax
#if     D1_D0_ENABLED
        integer(IK) , parameter     :: lenPattern = 1_IK
#elif   D0_D0_ENABLED || D1_D1_ENABLED
        integer(IK)                 :: lenPattern
#else
#error  "Unrecognized interface."
#endif
#if     CusIns_ENABLED
        integer(IK) , allocatable   :: locTemp(:) ! pattern Occurrence Position in the array.
        integer(IK)                 :: instanceMax, instanceLen, instanceNewLen
        integer(IK)                 :: instanceNew(size(instance))
        logical(LK)                 :: positive_def
        logical(LK)                 :: sorted_def
        instanceLen = size(instance, kind = IK)
        if (instanceLen == 0_IK) then
            allocate(loc(0))
            return
        end if
#elif   !DefIns_ENABLED
#error  "Unrecognized interface."
#endif
#if     D0_D0_ENABLED || D1_D1_ENABLED
        lenPattern = GET_SIZE(pattern, kind = IK)
        if (lenPattern == 0_IK) then
            allocate(loc(0))
            return
        end if
#endif
        lenArray = GET_SIZE(array, kind = IK)
        if (lenArray > lenPattern) then
            if (present(blindness)) then
                CHECK_ASSERTION(__LINE__, blindness /= 0_IK, SK_"@setLoc()/@getLoc(): The input `blindness` must be non-zero. blindness = "//getStr([blindness]))
                blindness_def = blindness
                countMax = lenArray / blindness_def + 1_IK
            else
                blindness_def = 1_IK
                countMax = lenArray - lenPattern + 1_IK
            end if
#if         CusIns_ENABLED
            positive_def = .false._LK
            if (present(positive)) positive_def = positive
            sorted_def = .false._LK
            if (present(sorted)) sorted_def = sorted
            if (sorted_def) then
                if (positive_def) then
                    if (instance(instanceLen) < countMax) countMax = instance(instanceLen)
                else
                    if (instance(instanceLen) < countMax .and. instance(instanceLen) > 0_IK) countMax = instance(instanceLen)
                end if
            else
                if (positive_def) then
                    instanceMax = maxval(instance)
                    if (instanceMax < countMax) countMax = instanceMax
                else
                    instanceMax = -huge(instanceMax)
                    do i = 1, instanceLen
                        if (instance(i) < 1_IK) then
                            instanceMax = -huge(instanceMax)
                            exit
                        elseif (instanceMax < instance(i)) then
                            instanceMax = instance(i)
                        end if
                    end do
                    if (instanceMax > 0_IK .and. instanceMax < countMax) countMax = instanceMax
                end if
            end if
#endif
            ! Find all requested instances of pattern.
            ! if (blindness_def > 0_IK) then
                allocate(loc(countMax))
                i = 1_IK
                count = 0_IK
                targetLen = lenArray - lenPattern + 1_IK
                loopFindPattern: do
                    if (ISEQ(array(GET_INDEX(i)),pattern)) then
                        count = count + 1_IK
                        loc(count) = i
                        i = i + blindness_def
                        if (count == countMax) exit loopFindPattern
                    else
                        i = i + 1_IK
                    end if
                    if (i > targetLen) exit loopFindPattern
                end do loopFindPattern
            ! else
            !     allocate(locTemp(countMax))
            !     i = lenArray - lenPattern + 1_IK
            !     count = 0_IK
            !     targetLen = 1_IK
            !     loopFindPatternReverse: do
            !         if (ISEQ(array(GET_INDEX(i)),pattern)) then
            !             count = count + 1_IK
            !             locTemp(count) = i
            !             i = i + blindness_def
            !             if (count == countMax) exit loopFindPatternReverse
            !         else
            !             i = i - 1_IK
            !         end if
            !         if (i < targetLen) exit loopFindPatternReverse
            !     end do loopFindPatternReverse
            !     loc = locTemp(count:1:-1)
            !#if CusIns_ENABLED
            !   return
            !#endif
            ! end if
#if         CusIns_ENABLED
            ! Find array at all requested instances of pattern.
            blockPatternFound: if (count > 0_IK) then
                ! Convert all negative and positive instances to counts from the beginning within the possible range [1, count].
                instanceNewLen = 0_IK
                i = 0_IK
                do ! This loop requires instanceLen to be at least 1, which is guaranteed by the condition after `instanceLen` definition in the above.
                    i = i + 1_IK
                    if (instance(i) > 0_IK .and. instance(i) <= count) then
                        instanceNewLen = instanceNewLen  + 1_IK
                        instanceNew(instanceNewLen) = instance(i)
                    elseif (instance(i) < 0_IK .and. instance(i) + count + 1_IK > 0_IK) then
                        instanceNewLen = instanceNewLen  + 1_IK
                        instanceNew(instanceNewLen) = instance(i) + count + 1_IK
                    end if
                    if (i == instanceLen) exit
                end do
                !sorted_def = .false._LK
                !if (present(sorted)) sorted_def = sorted
                !if (.not. sorted_def) call setSorted(instanceNew(1:instanceNewLen))
                !positive_def = .false._LK
                !if (present(positive)) positive_def = positive
                !if (positive_def) then
                    locTemp = loc(instanceNew(1:instanceNewLen))
                    call move_alloc(locTemp, loc)
                !else
                !    block
                !        use pm_arrayUnique, only: getUnique
                !        instanceNew = getUnique(instanceNew(1:instanceNewLen))
                !    end block
                !    count = size(instanceNew, kind = IK)
                !end if
                !deallocate(locTemp)
                return
            end if blockPatternFound
#endif
            loc = loc(1:count)
            return
        elseif (lenArray == lenPattern & ! LCOV_EXCL_LINE
#if     CusIns_ENABLED
        .and. any(abs(instance) == 1_IK) & ! LCOV_EXCL_LINE
#endif
       ) then
            if (ISEQ(array(GET_INDEX(1)),pattern)) then
                loc = [1_IK]
                return
            end if
        end if
        allocate(loc(0))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setLoc_ENABLED && DefIns_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: lenLoc, lenArray, i, targetLen
#if     D0_D0_ENABLED
        integer(IK), parameter :: offset = 0_IK ! for lower bound of allocatable `loc`.
#elif   D1_D0_ENABLED || D1_D1_ENABLED
        integer(IK) :: offset
#endif
#if     D1_D0_ENABLED
        integer(IK), parameter :: lenPattern = 1_IK
#elif   D0_D0_ENABLED || D1_D1_ENABLED
        integer(IK) :: lenPattern
        lenPattern = GET_SIZE(pattern, kind = IK)
        if (lenPattern == 0_IK) then
            nloc = 0_IK
            return
        end if
#else
#error  "Unrecognized interface."
#endif
        lenArray = GET_SIZE(array, kind = IK)
        CHECK_ASSERTION(__LINE__, 0_IK < blindness, SK_"@setLoc(): The condition `0 < blindness` must hold. blindness = "//getStr([blindness]))
        CHECK_ASSERTION(__LINE__, allocated(loc), SK_"@setLoc(): The condition `allocated(loc)` must hold. allocated(loc) = "//getStr(allocated(loc)))
        lenLoc = size(loc, 1, IK)
#if     D1_D0_ENABLED || D1_D1_ENABLED
        offset = lbound(loc, 1, IK) - 1_IK
#endif
        if (lenPattern < lenArray) then
            i = 1_IK
            nloc = 0_IK
            targetLen = lenArray - lenPattern + 1_IK
            loopFindPattern: do
                if (ISEQ(array(GET_INDEX(i)), pattern)) then
                    nloc = nloc + 1_IK
                    if (lenLoc < nloc) then
                        lenLoc = 2 * nloc
                        call setResized(loc, lenLoc)
                    end if
                    loc(nloc + offset) = i
                    i = i + blindness
                else
                    i = i + 1_IK
                end if
                if (targetLen < i) exit loopFindPattern
            end do loopFindPattern
            return
        elseif (lenArray == lenPattern) then
            if (ISEQ(array(GET_INDEX(1)), pattern)) then
                if (lenLoc < 1_IK) call setResized(loc, 1_IK)
                loc(1 + offset) = 1_IK
                nloc = 1_IK
                return
            end if
        end if
        nloc = 0_IK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setLoc_ENABLED && CusIns_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: lenLoc, lenArray, i, targetLen, instanceMax, instanceLen, instanceNew(size(instance))
#if     D1_D0_ENABLED
        integer(IK), parameter :: lenPattern = 1_IK
#elif   D0_D0_ENABLED || D1_D1_ENABLED
        integer(IK) :: lenPattern
#else
#error  "Unrecognized interface."
#endif
#if     D0_D0_ENABLED
        integer(IK), parameter :: offset = 0_IK ! for lower bound of allocatable `loc`.
#elif   D1_D0_ENABLED || D1_D1_ENABLED
        integer(IK) :: offset
#endif
        nloc = 0_IK
        CHECK_ASSERTION(__LINE__, allocated(loc), SK_"@setLoc(): The condition `allocated(loc)` must hold. allocated(loc) = "//getStr(allocated(loc)))
        CHECK_ASSERTION(__LINE__, 0_IK < blindness, SK_"@setLoc(): The condition `0 < blindness` must hold. blindness = "//getStr([blindness]))
#if     D1_D0_ENABLED || D1_D1_ENABLED
        offset = lbound(loc, 1, IK) - 1_IK
#endif
        instanceLen = size(instance, kind = IK)
        if (instanceLen == 0_IK) return
#if     D0_D0_ENABLED || D1_D1_ENABLED
        lenPattern = GET_SIZE(pattern, kind = IK)
        if (lenPattern == 0_IK) return
#endif
        lenLoc = size(loc, 1, IK)
        lenArray = GET_SIZE(array, kind = IK)
        if (lenPattern < lenArray) then
            if (positive) then
                CHECK_ASSERTION(__LINE__, all(0_IK < instance), SK_"@setLoc(): The condition `all(0 < instance)` must hold. instance = "//getStr([instance]))
                if (positive .and. sorted) then
                    CHECK_ASSERTION(__LINE__, isAscending(instance), SK_"@setLoc(): The condition `isSorted(instance)` must hold. instance = "//getStr([instance]))
                    instanceMax = instance(instanceLen)
                else
                    instanceMax = maxval(instance, 1)
                end if
                ! find all patterns up to `instanceMax`.
                i = 1_IK
                targetLen = lenArray - lenPattern + 1_IK
                loopFindPattern: do
                    if (instanceMax == nloc) exit loopFindPattern ! This is the only extra line compared to the non-positive `instance`.
                    if (ISEQ(array(GET_INDEX(i)), pattern)) then
                        nloc = nloc + 1_IK
                        if (lenLoc < nloc) then
                            lenLoc = 2 * nloc
                            call setResized(loc, lenLoc)
                        end if
                        loc(nloc + offset) = i
                        i = i + blindness
                    else
                        i = i + 1_IK
                    end if
                    if (targetLen < i) exit loopFindPattern
                end do loopFindPattern
                if (0_IK < nloc) then
                    ! Clean up `instance` from nonsense values.
                    targetLen = 0_IK
                    do i = 1, instanceLen
                        if (instance(i) <= nloc) then
                            targetLen = targetLen + 1_IK
                            instanceNew(targetLen) = instance(i) + offset
                        end if
                    end do
                    if (nloc < targetLen) call setResized(loc, targetLen)
                    loc(1 + offset : targetLen + offset) = loc(instanceNew(1 : targetLen))
                    nloc = targetLen
                end if
            else
                ! Assume the worst case scenario instanceMax = lenArray - lenPattern + 1_IK.
                i = 1_IK
                targetLen = lenArray - lenPattern + 1_IK
                loopFindPatternAll: do
                    if (ISEQ(array(GET_INDEX(i)), pattern)) then
                        nloc = nloc + 1_IK
                        if (lenLoc < nloc) then
                            lenLoc = 2 * nloc
                            call setResized(loc, lenLoc)
                        end if
                        loc(nloc + offset) = i
                        i = i + blindness
                    else
                        i = i + 1_IK
                    end if
                    if (targetLen < i) exit loopFindPatternAll
                end do loopFindPatternAll
                !print *, nloc
                if (0_IK < nloc) then
                    ! Clean up `instance` from nonsense values.
                    targetLen = 0_IK
                    do i = 1, instanceLen
                        if (0_IK < instance(i) .and. instance(i) <= nloc) then
                            targetLen = targetLen + 1_IK
                            instanceNew(targetLen) = instance(i) + offset
                        else ! ensure the specified instance is a valid negative value.
                            instanceNew(targetLen + 1) = instance(i) + nloc + 1_IK + offset
                            if (0_IK < instanceNew(targetLen + 1) .and. instanceNew(targetLen + 1) <= nloc) targetLen = targetLen + 1_IK
                        end if
                    end do
                    if (nloc < targetLen) call setResized(loc, targetLen)
                    loc(1 + offset : targetLen + offset) = loc(instanceNew(1 : targetLen))
                    nloc = targetLen
                    !print *, targetLen
                    !print *, instanceNew(1 : nloc)
                end if
            end if
        elseif (lenArray == lenPattern .and. any(abs(instance) == 1_IK)) then
            if (ISEQ(array(GET_INDEX(1)), pattern)) then
                if (size(loc, 1, IK) < 1_IK) call setResized(loc, 1_IK)
                loc(1 + offset) = 1_IK
                nloc = 1_IK
                return
            end if
        end if

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  INSTANCENEW
#undef  GET_INDEX
#undef  GET_SIZE
#undef  IS_EQUAL
#undef  ISEQ