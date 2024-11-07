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
!>  This file contains the implementation details of the routines under the generic interfaces of [pm_arrayReplace](@ref pm_arrayReplace).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define logical vs. normal equivalence operators
#if     LK_ENABLED
#define IS_EQUAL .eqv.
#elif   SK_ENABLED || IK_ENABLED || CK_ENABLED || RK_ENABLED
#define IS_EQUAL ==
#else
#error  "Unrecognized interface."
#endif
        ! Define scalar vs. vector operations.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     D0_D0_D0_ENABLED && SK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#define GET_SIZE len
#define GET_INDEX(i) i : i + lenPattern - 1_IK
#define D0_D0_D0_ENABLED 1
#if     CusCom_ENABLED
#define ISEQ(segment, pattern) iseq(segment, pattern)
#elif   DefCom_ENABLED
#define ISEQ(segment, pattern) segment == pattern
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   D1_D0_D0_ENABLED || D1_D0_D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#define GET_SIZE size
#define GET_INDEX(i) i
#if     CusCom_ENABLED
#define ISEQ(segment, pattern) iseq(segment, pattern)
#elif   DefCom_ENABLED
#define ISEQ(segment, pattern) segment IS_EQUAL pattern
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   D1_D1_D0_ENABLED || D1_D1_D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#define GET_SIZE size
#define GET_INDEX(i) i : i + lenPattern - 1_IK
#if     CusCom_ENABLED
#define ISEQ(segment,pattern) iseq(segment, pattern, lenPattern)
#elif   DefCom_ENABLED
#define ISEQ(segment,pattern) all(segment IS_EQUAL pattern)
#else
#error  "Unrecognized interface."
#endif
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
        ! Define the temporary new array for cases where the result is to be returned in the input array.
#if     setReplaced_ENABLED && D0_D0_D0_ENABLED && SK_ENABLED
        character(:,SKG)            , allocatable :: arrayNew
#elif   setReplaced_ENABLED && (D1_D0_D0_ENABLED || D1_D1_D0_ENABLED || D1_D0_D1_ENABLED || D1_D1_D1_ENABLED)
#if     SK_ENABLED
        character(len(array,IK),SKG), allocatable :: arrayNew(:)
#elif   IK_ENABLED
        integer(IKG)                , allocatable :: arrayNew(:)
#elif   LK_ENABLED
        logical(LKG)                , allocatable :: arrayNew(:)
#elif   CK_ENABLED
        complex(CKG)                , allocatable :: arrayNew(:)
#elif   RK_ENABLED
        real(RKG)                   , allocatable :: arrayNew(:)
#else
#error  "Unrecognized interface."
#endif
#elif   !getReplaced_ENABLED
#error  "Unrecognized interface."
#endif
        ! Declare local variables.
#if     CusIns_ENABLED
        integer(IK)                 :: lenInstance, lenInstanceNew, maxInstance!, minInstance
        integer(IK) , allocatable   :: instanceNew(:)
        logical(LK)                 :: sorted_def
        logical(LK)                 :: unique_def
#endif
        integer(IK) , allocatable   :: POP(:) ! pattern Occurrence Position in the array.
        integer(IK)                 :: lenArray, lenDiff, i, iLast
        integer(IK)                 :: lenArrayNew, newPOP, newPOPNext, lenPOP, lenPOPMax
        ! Declare the replacement length.
#if     D1_D0_D0_ENABLED || D1_D1_D0_ENABLED
        integer(IK) , parameter     :: lenReplacement = 1_IK
#elif   D0_D0_D0_ENABLED || D1_D1_D1_ENABLED || D1_D0_D1_ENABLED
#define lenReplacement_ENABLED 1
        integer(IK)                 :: lenReplacement
#else
#error  "Unrecognized interface."
#endif
        ! Declare the pattern length.
#if     D1_D0_D0_ENABLED || D1_D0_D1_ENABLED
        integer(IK) , parameter     :: lenPattern = 1_IK
#elif   D0_D0_D0_ENABLED || D1_D1_D1_ENABLED || D1_D1_D0_ENABLED
#define lenPattern_ENABLED 1
        integer(IK)                 :: lenPattern
#else
#error  "Unrecognized interface."
#endif
        ! Set the array offset.
#if     D0_D0_D0_ENABLED || getReplaced_ENABLED
        integer(IK) , parameter     :: offset = 0_IK
#elif   setReplaced_ENABLED
        integer(IK)                 :: offset
        offset = lbound(array,1,IK) - 1_IK
#else
#error  "Unrecognized interface."
#endif
        ! Set the replacement length.
#if     lenReplacement_ENABLED
        lenReplacement = GET_SIZE(replacement, kind = IK) ! \warning GET_SIZE is a preprocessor macro.
#endif
        ! Set the pattern length.
#if     lenPattern_ENABLED
        lenPattern = GET_SIZE(pattern, kind = IK) ! \warning GET_SIZE is a preprocessor macro.
#endif
        lenArray = GET_SIZE(array, kind = IK) ! \warning GET_SIZE is a preprocessor macro.
#if     CusIns_ENABLED
        lenInstance = size(instance, kind = IK)
        if (lenInstance == 0_IK) then
#if         getReplaced_ENABLED
            arrayNew = array
#elif       !setReplaced_ENABLED
#error      "Unrecognized Interface."
#endif
            return
        end if
#endif
        if (lenArray > lenPattern) then
            blockFullEmptyPattern: if (lenPattern > 0_IK) then
                lenPOPMax = lenArray / lenPattern + 1_IK
#if             CusIns_ENABLED
                !print *, "instance", instance
                maxInstance = maxval(instance)
                if (minval(instance) >= 0_IK .and. maxInstance < lenPOPMax) lenPOPMax = maxInstance
#endif
                !print *, "array", array
                !print *, "pattern", pattern
                !print *, "lenArray, offset, lenPattern", lenArray, offset, lenPattern
                ! Find all requested instances of pattern.
                allocate(POP(lenPOPMax))!, source = -huge(1_IK))
                i = 1_IK + offset
                lenPOP = 0_IK
                iLast = lenArray + offset - lenPattern + 1_IK
                loopFindPOP: do
#if                 getReplaced_ENABLED && CusCom_ENABLED && CusIns_ENABLED && D1_D1_D1_ENABLED && RK_ENABLED
                    !!  \bug
                    !!  gfortran 11 cannot correctly pass the length of the input `array` argument to `iseq()`
                    !!  via an explicit interface (and so why `iseq()` interface remains implicit.
                    !print *, array(GET_INDEX(i)), pattern
                    !print *, size(array(GET_INDEX(i))), size(pattern)
                    !print *, ISEQ(array(GET_INDEX(i)), pattern)
#endif
                    if (ISEQ(array(GET_INDEX(i)), pattern)) then ! fpp
                        lenPOP = lenPOP + 1_IK
                        if (lenPOP > lenPOPMax) exit loopFindPOP ! This condition is crucial when `maxInstance < lenPOPMax`.
                        !print *, "POP", POP
                        POP(lenPOP) = i
                        i = i + lenPattern
                    else
                        i = i + 1_IK
                    end if
                    if (i > iLast) exit loopFindPOP
                end do loopFindPOP
            else blockFullEmptyPattern
                lenPOP = lenArray + 1_IK
#if             CusIns_ENABLED
                maxInstance = maxval(instance)
                if (minval(instance) >= 0_IK .and. maxInstance < lenPOP) lenPOP = maxInstance
#endif
                allocate(POP(lenPOP))
                do i = 1, lenPOP
                    POP(i) = i + offset
                end do
            end if blockFullEmptyPattern
            ! Replace all requested instances of pattern.
            blockInstanceExists: if (lenPOP > 0_IK) then
#if             CusIns_ENABLED
                ! Convert all negative and positive instances to counts from the beginning within the possible range [1, lenPOP].
                !lenInstance = size(instance, kind = IK) ! this is now moved up to quit if zero-length instance is encountered.
                allocate(instanceNew(lenInstance))
                lenInstanceNew = 0_IK
                i = 0_IK
                ! This loop requires lenInstance to be at least 1, which is guaranteed by the condition after `lenInstance` definition in the above.
                do
                    i = i + 1_IK
                    if (instance(i) > 0_IK .and. instance(i) <= lenPOP) then
                        lenInstanceNew = lenInstanceNew  + 1_IK
                        instanceNew(lenInstanceNew) = instance(i)
                    elseif (instance(i) < 0_IK .and. instance(i) + lenPOP + 1_IK > 0_IK) then
                        lenInstanceNew = lenInstanceNew  + 1_IK
                        instanceNew(lenInstanceNew) = instance(i) + lenPOP + 1_IK
                    end if
                    if (i == lenInstance) exit
                end do
                sorted_def = .false._LK
                if (present(sorted)) sorted_def = sorted
                if (.not. sorted_def) call setSorted(instanceNew(1:lenInstanceNew))
                unique_def = .false._LK
                if (present(unique)) unique_def = unique
                if (unique_def) then
                    lenPOP = lenInstanceNew
                else
                    instanceNew = getUnique(instanceNew(1:lenInstanceNew))
                    lenPOP = size(instanceNew, kind = IK)
                end if
                if (lenPOP == 0_IK) then ! instance is empty, return the input array, untouched.
#if                 getReplaced_ENABLED
                    arrayNew = array
#endif
                    ! The following deallocations are essential since gfortran, as of version 10.3, cannot automatically deallocate array upon return.
#if                 CusIns_ENABLED
                    deallocate(instanceNew)
#endif
                    deallocate(POP)
                    return
                end if
#define         INSTANCENEW(i) instanceNew(i)
#elif           DefIns_ENABLED
#define         INSTANCENEW(i) i
#else
#error          "Unrecognized Interface."
#endif
                lenDiff = lenReplacement - lenPattern
                lenArrayNew = lenArray + lenPOP * lenDiff
#if             SK_ENABLED && D0_D0_D0_ENABLED
                allocate(character(lenArrayNew,SKG) :: arrayNew)
                !>  \bug
                !>  This string vector allocation must be separated from the following because of a bug in Intel ifort 2021.5.
                !>  The bug is related to the separation of module interface from implementation.
#elif           SK_ENABLED && getReplaced_ENABLED && (D1_D0_D0_ENABLED || D1_D0_D1_ENABLED || D1_D1_D0_ENABLED || D1_D1_D1_ENABLED)
                allocate(character(len(array,IK),SKG) :: arrayNew(1_IK + offset : lenArrayNew + offset))
#else
                allocate(arrayNew(1_IK + offset : lenArrayNew + offset))
#endif
!#if             getReplacedDefComCusIns_D1_D0_D1_IK_ENABLED || getReplacedDefComCusIns_D1_D1_D1_IK_ENABLED
                !print *, "size(replacement)", size(replacement)
                !print *, "instanceNew", instanceNew
                !print *, "INSTANCENEW(1_IK)", INSTANCENEW(1_IK)
                !print *, "offset", offset
                !print *, "POP", POP
!#endif
                newPOP = POP(INSTANCENEW(1_IK))
                arrayNew(1_IK + offset : POP(INSTANCENEW(1_IK)) - 1_IK) = array(1_IK + offset : POP(INSTANCENEW(1_IK)) - 1_IK)
                do i = 1_IK, lenPOP - 1_IK
                    arrayNew(newPOP : newPOP + lenReplacement - 1_IK) = replacement
                    newPOPNext = POP(INSTANCENEW(i + 1_IK)) + i * lenDiff
                    arrayNew(newPOP + lenReplacement : newPOPNext - 1_IK) = array(POP(INSTANCENEW(i)) + lenPattern : POP(INSTANCENEW(i+1_IK)) - 1_IK)
                    newPOP = newPOPNext
                end do
                arrayNew(newPOP : newPOP + lenReplacement - 1_IK) = replacement
                arrayNew(newPOP + lenReplacement : lenArrayNew + offset) = array(POP(INSTANCENEW(i)) + lenPattern : lenArray + offset)
#if             CusIns_ENABLED
                deallocate(instanceNew) ! This is essential since gfortran, as of version 10.3, cannot automatically deallocate array upon return.
#endif
#if             setReplaced_ENABLED
                call move_alloc(from = arrayNew, to = array)
#elif           getReplaced_ENABLED
            else blockInstanceExists
                arrayNew = array
#else
#error          "Unrecognized interface."
#endif
            end if blockInstanceExists
            deallocate(POP)
        elseif (lenArray == lenPattern) then
            if (ISEQ(array(GET_INDEX(1_IK + offset)), pattern)) then
#if             CusIns_ENABLED
                !   \bug
                !   Bizarrely, if this condition is merged with the above, then both ifort and gfortran occasionally
                !   (but in different situations yield .true., even when expression is `.true. and .false.`.
                if (any(abs(instance) == 1_IK)) then
#endif
#if                 setReplaced_ENABLED && D0_D0_D0_ENABLED && SK_ENABLED
                    array = replacement
#elif               setReplaced_ENABLED
                    deallocate(array)
                    allocate(array(1_IK + offset : lenReplacement + offset), source = replacement)
#elif               getReplaced_ENABLED
#if                 D0_D0_D0_ENABLED && SK_ENABLED
                    arrayNew = replacement
#else
                    allocate(arrayNew(1_IK + offset : lenReplacement + offset), source = replacement)
#endif
#else
#error              "Unrecognized Interface."
#endif
                    return
#if             CusIns_ENABLED
                end if
#endif
            end if
#if         getReplaced_ENABLED
            arrayNew = array
        else ! array is smaller than pattern.
            arrayNew = array
#endif
        end if
#undef  lenReplacement_ENABLED
#undef  lenPattern_ENABLED
#undef  INSTANCENEW
#undef  GET_INDEX
#undef  GET_SIZE
#undef  IS_EQUAL
#undef  ISEQ
#undef  ANY
#undef  ALL