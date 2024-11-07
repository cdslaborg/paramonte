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
!>  This file contains procedure implementations of [pm_arraySort](@ref pm_arraySort).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the custom comparison macro for recursive sorting.
#if     DefCom_ENABLED
#define ISSORTED
#elif   CusCom_ENABLED
#define ISSORTED , isSorted
#else
#error  "Unrecognized interface."
#endif
        ! Set the sorting rules.
#if     CusCom_ENABLED
#define IS_SORTED(i,j) isSorted(i,j)
#elif   LK_ENABLED && DefCom_ENABLED
#define IS_SORTED(i,j) j .and. .not. i
#elif   CK_ENABLED && DefCom_ENABLED
#define IS_SORTED(i,j) i%re < j%re
#elif   (PSSK_ENABLED || BSSK_ENABLED) && DefCom_ENABLED
#define IS_SORTED(i,j) i%val < j%val
#elif   (SK_ENABLED || IK_ENABLED || RK_ENABLED) && DefCom_ENABLED
#define IS_SORTED(i,j) i < j
#else
#error  "Unrecognized interface."
#endif
        ! Set the indexing rules.
#if     D0_ENABLED && SK_ENABLED
#define GET_INDEX(i) i:i
#define GET_SIZE len
#elif   D1_ENABLED
#define GET_INDEX(i) i
#define GET_SIZE size
#else
#error  "Unrecognized interface."
#endif
        ! Set the types and kinds.
#if     SK_ENABLED && D0_ENABLED
#define TYPE_KIND character(1,SKG)
#elif   SK_ENABLED && D1_ENABLED
#define TYPE_KIND character(len(array,IK),SKG)
#elif   IK_ENABLED && D1_ENABLED
#define TYPE_KIND integer(IKG)
#elif   LK_ENABLED && D1_ENABLED
#define TYPE_KIND logical(LKG)
#elif   CK_ENABLED && D1_ENABLED
#define TYPE_KIND complex(CKG)
#elif   RK_ENABLED && D1_ENABLED
#define TYPE_KIND real(RKG)
#elif   PSSK_ENABLED && D1_ENABLED
#define TYPE_KIND type(css_pdt(SKG))
#elif   BSSK_ENABLED && D1_ENABLED
#define TYPE_KIND type(css_type)
#else
#error  "Unrecognized interface."
#endif
        ! Set the len type parameter for the string kind.
#if     SK_ENABLED
#define TYPE_KIND_LEN(LEN) character(LEN,SKG)
#else
#define TYPE_KIND_LEN(LEN) TYPE_KIND
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getSorted_ENABLED && Ind_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     DefCom_ENABLED
        call setSorted(array, sorting)
#elif   CusCom_ENABLED
        call setSorted(array, sorting, isSorted)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getSorted_ENABLED && Arr_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     DefCom_ENABLED
#define CALL_SETSORTED(METHOD) call setSorted(sorting, METHOD)
#elif   CusCom_ENABLED
#define CALL_SETSORTED(METHOD) call setSorted(sorting, isSorted, METHOD)
#else
#error  "Unrecognized interface."
#endif
        sorting = array
        blockMethod: if (present(method)) then
            select type (method)
            type is (qsorti_type)
                exit blockMethod
            type is (qsortr_type)
                CALL_SETSORTED(method)
            type is (qsortrdp_type)
                CALL_SETSORTED(method)
            type is (bubble_type)
                CALL_SETSORTED(method)
            type is (heapi_type)
                CALL_SETSORTED(method)
            type is (heapr_type)
                CALL_SETSORTED(method)
            type is (insertionl_type)
                CALL_SETSORTED(method)
            type is (insertionb_type)
                CALL_SETSORTED(method)
            type is (merger_type)
                CALL_SETSORTED(method)
            type is (selection_type)
                CALL_SETSORTED(method)
            type is (shell_type)
                CALL_SETSORTED(method)
            class default
                error stop MODULE_NAME//SK_"@getSorted(): Unrecognized `method`."
            end select
            return
        end if blockMethod
        CALL_SETSORTED(qsorti)
#undef  CALL_SETSORTED

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setSorted_ENABLED && Ind_ENABLED && (Qsorti_ENABLED || Def_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_KIND :: temp
        integer(IK), parameter :: MIN_ARRAY_SIZE = 15_IK
        integer(IK) :: i, j, swapIndex, tempIndex, low, high, mid, stackCounter, stack(2_IK * bit_size(0_IK)) ! storage_size(0_IK)
        high = GET_SIZE(array, kind = IK)
        do concurrent(j = 1 : high); index(j) = j; end do
        CHECK_ASSERTION(__LINE__, high == size(index, 1, IK), SK_"@setSorted(): The condition `size(array) == size(index)` must hold. size(array), size(index) = "//getStr([high, size(index, 1, IK)])) ! fpp
        !if (high == 0_IK) return
        stackCounter = 0_IK
        low = 1_IK
        loopMain: do
            if (high - low < MIN_ARRAY_SIZE) then
                do j = low + 1_IK, high
                    tempIndex = index(j)
                    temp = array(GET_INDEX(tempIndex))
                    do i = j - 1_IK, low, -1_IK
                        if (IS_SORTED(temp, array(GET_INDEX(index(i))))) then ! fpp
                            index(i + 1_IK) = index(i)
                            cycle
                        end if
                        exit
                        !if (array(GET_INDEX(index(i))) <= temp) exit ! fpp
                        !index(i + 1_IK) = index(i)
                    end do
                    index(i + 1_IK) = tempIndex
                end do
                if (stackCounter == 0_IK) exit loopMain
                high = stack(stackCounter)
                low = stack(stackCounter - 1_IK)
                stackCounter = stackCounter - 2_IK
            else
                mid = (low + high) / 2_IK
                tempIndex = index(mid)
                index(mid) = index(low + 1_IK)
                index(low + 1_IK) = tempIndex
                if (IS_SORTED(array(GET_INDEX(index(high))) , array(GET_INDEX(index(low))))) then ! fpp
                    swapIndex = index(low)
                    index(low) = index(high)
                    index(high) = swapIndex
                end if
                if (IS_SORTED(array(GET_INDEX(index(high))) , array(GET_INDEX(index(low + 1_IK))))) then ! fpp
                    swapIndex = index(low + 1_IK)
                    index(low + 1_IK) = index(high)
                    index(high) = swapIndex
                end if
                if (IS_SORTED(array(GET_INDEX(index(low + 1_IK))) , array(GET_INDEX(index(low))))) then ! fpp
                    swapIndex = index(low)
                    index(low) = index(low + 1_IK)
                    index(low + 1_IK) = swapIndex
                end if
                !call exchangeIndex(index(low),index(high))
                !call exchangeIndex(index(low + 1_IK),index(high))
                !call exchangeIndex(index(low),index(low + 1_IK))
                i = low + 1_IK
                j = high
                tempIndex = index(low + 1_IK)
                temp = array(GET_INDEX(tempIndex)) ! fpp
                do
                    do
                        i = i + 1_IK
                        !if (temp <= array(GET_INDEX(index(i))) ) exit ! fpp
                        if (IS_SORTED(array(GET_INDEX(index(i))) , temp ) ) cycle ! fpp
                        exit
                    end do
                    do
                        j = j - 1_IK
                        !if (array(GET_INDEX(index(j))) <= temp) exit
                        if (IS_SORTED(temp, array(GET_INDEX(index(j))))) cycle ! fpp
                        exit
                    end do
                    if (j < i) exit
                    swapIndex = index(i)
                    index(i) = index(j)
                    index(j) = swapIndex
                end do
                index(low + 1_IK) = index(j)
                index(j) = tempIndex
                stackCounter = stackCounter + 2_IK
                CHECK_ASSERTION(__LINE__, size(stack, kind = IK) > stackCounter, SK_"@setSorted(): The stack size exceeded. This is highly unusual. size(stack), stackCounter = "//getStr([size(stack, kind = IK), stackCounter])) ! fpp
                if (j - low <= high - i + 1_IK) then
                    stack(stackCounter) = high
                    stack(stackCounter - 1_IK) = i
                    high = j - 1_IK
                else
                    stack(stackCounter) = j - 1_IK
                    stack(stackCounter - 1_IK) = low
                    low = i
                end if
            end if
        end do loopMain
!    contains
!        pure subroutine exchangeIndex(i,j)
!            integer(IK), intent(inout) :: i,j
!            integer(IK)                :: swp
!            if (array(GET_INDEX(j)) < array(GET_INDEX(i))) then
!                swp = i
!                i = j
!                j = swp
!            end if
!        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setSorted_ENABLED && Arr_ENABLED && (Qsorti_ENABLED || Def_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_KIND :: temp, pivot
        integer :: stackCounter
        integer(IK) :: i, j, left, right, low, high, mid
        integer(IK) :: stack(2_IK * bit_size(0_IK)) ! storage_size(i)
#define UseSelection_ENABLED 1
#if     UseSelection_ENABLED
        integer(IK) :: k
#endif
        low = 1_IK
        high = GET_SIZE(array, kind = IK)
        stackCounter = 1_IK
        do
            if (high - low < 30_IK) then
#if             UseSelection_ENABLED
                ! use setSortedSelection on small Arrays. Works better than setSortedInsertion.
                loopSelectionSort: do i = low, high - 1_IK
#if                 1
                    temp = array(GET_INDEX(i))
                    j = i
                    do k = i + 1_IK, high ! minloc
                        if (IS_SORTED(array(GET_INDEX(k)), temp)) then
                            temp = array(GET_INDEX(k))
                            j = k
                        end if
                    end do
#else
                    j = minloc(array(i:high), dim = 1, kind = IK) + i - 1_IK
#endif
                    temp = array(GET_INDEX(i))
                    array(GET_INDEX(i)) = array(GET_INDEX(j))
                    array(GET_INDEX(j)) = temp
                end do loopSelectionSort
#else
                ! use setSortedInsertion on small Arrays.
                loopInsertionSort: do i = low + 1_IK, high
                    temp = array(GET_INDEX(i))
                    do j = i - 1_IK, low, -1_IK
                        if (IS_SORTED(temp, array(GET_INDEX(j)))) then
                            array(GET_INDEX(j + 1_IK)) = array(GET_INDEX(j))
                            cycle
                        end if
                        exit
                    end do
                    array(GET_INDEX(j + 1_IK)) = temp
                end do loopInsertionSort
#endif
                ! pop from stack.
                if (stackCounter == 1) return
                stackCounter = stackCounter - 2
                high = stack(stackCounter+1)
                low = stack(stackCounter)
                cycle
            end if
            ! Find median of three pivot and place sentinels at first and last elements.
            mid = (low + high) / 2_IK
            left = low + 1_IK
            temp = array(GET_INDEX(mid))
            array(GET_INDEX(mid)) = array(GET_INDEX(left))
            if (IS_SORTED(array(GET_INDEX(high)), temp)) then
                array(GET_INDEX(left)) = array(GET_INDEX(high))
                array(GET_INDEX(high)) = temp
            else
                array(GET_INDEX(left)) = temp
            end if
            if (IS_SORTED(array(GET_INDEX(high)) , array(GET_INDEX(low)))) then
                temp = array(GET_INDEX(low))
                array(GET_INDEX(low)) = array(GET_INDEX(high))
                array(GET_INDEX(high)) = temp
            end if
            if (IS_SORTED(array(GET_INDEX(left)) , array(GET_INDEX(low)))) then
                temp = array(GET_INDEX(low))
                array(GET_INDEX(low)) = array(GET_INDEX(left))
                array(GET_INDEX(left)) = temp
            end if
            pivot = array(GET_INDEX(left))
            right = high - 1_IK
            left = left + 1_IK
            do
                if (IS_SORTED(array(GET_INDEX(left)), pivot))then
                    do
                        left = left + 1_IK
                        if (IS_SORTED(array(GET_INDEX(left)), pivot))cycle
                        exit
                    end do
                end if
                if (IS_SORTED(pivot , array(GET_INDEX(right)))) then
                    do
                        right = right - 1_IK
                        if (IS_SORTED(pivot , array(GET_INDEX(right)))) cycle
                        exit
                    end do
                end if
                if (left >= right) exit
                temp = array(GET_INDEX(left))
                array(GET_INDEX(left)) = array(GET_INDEX(right))
                array(GET_INDEX(right)) = temp
                right = right - 1_IK
                left = left + 1_IK
            end do
            if (left == right) left = left + 1_IK
            if (left < mid) then
                stack(stackCounter) = left
                stack(stackCounter + 1) = high
                stackCounter = stackCounter + 2
                high = left - 1_IK
            else
                stack(stackCounter) = low
                stack(stackCounter + 1) = left - 1_IK
                stackCounter = stackCounter + 2
                low = left
            end if
        end do
#undef  UseSelection_ENABLED

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setSorted_ENABLED && Arr_ENABLED && Qsortr_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_KIND :: temp, pivot
        integer(IK) :: i, j, left, lenArray
        lenArray = GET_SIZE(array, kind = IK)
        if(lenArray > 15_IK) then
            ! partition.
            i = 0_IK
            pivot = array(GET_INDEX(1_IK))
            j = GET_SIZE(array, kind = IK) + 1_IK
            do
                j = j - 1_IK
                do
                    if (IS_SORTED(pivot, array(GET_INDEX(j)))) then
                        j = j - 1_IK
                        cycle
                    end if
                    exit
                end do
                i = i + 1_IK
                do
                    if (IS_SORTED(array(GET_INDEX(i)), pivot)) then
                        i = i + 1_IK
                        cycle
                    end if
                    exit
                end do
                if (i < j) then ! exchange array(i) and array(j)
                    temp = array(GET_INDEX(i))
                    array(GET_INDEX(i)) = array(GET_INDEX(j))
                    array(GET_INDEX(j)) = temp
                elseif (i == j) then
                    left = i + 1_IK
                    exit
                else
                    left = i
                    exit
                endif
            end do
            call setSorted(array(1 : left - 1_IK)ISSORTED, qsortr)
            call setSorted(array(left : lenArray)ISSORTED, qsortr)
        elseif (lenArray > 1_IK) then ! use insertion sort on small Arrays
            do i = 2_IK, lenArray
                temp = array(GET_INDEX(i))
                do j = i - 1_IK, 1_IK, -1_IK
                    if (IS_SORTED(temp, array(GET_INDEX(j)))) then
                        array(GET_INDEX(j + 1_IK)) = array(GET_INDEX(j))
                        cycle
                    end if
                    exit
                end do
                array(GET_INDEX(j + 1_IK)) = temp
            end do
        endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setSorted_ENABLED && Arr_ENABLED && Qsortrdp_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_KIND :: temp, pivot1, pivot2
        integer(IK) :: i, j, last, l, k, g
        last = GET_SIZE(array, kind = IK)
        if (last < 15_IK) then ! use insertion sort on small Arrays
            do i = 2_IK, last
                temp = array(GET_INDEX(i))
                do j = i - 1_IK, 1_IK, -1_IK
                    if (IS_SORTED(temp, array(GET_INDEX(j)))) then
                        array(GET_INDEX(j + 1)) = array(GET_INDEX(j))
                        cycle
                    end if
                    exit
                end do
                array(GET_INDEX(j + 1)) = temp
            end do
            return
        end if
        pivot1 = array(GET_INDEX(last / 3_IK))
        pivot2 = array(GET_INDEX(2_IK * last / 3_IK))
        if (IS_SORTED(pivot2, pivot1)) then
            temp = pivot1
            pivot1 = pivot2
            pivot2 = temp
        end if
        array(GET_INDEX(last / 3_IK)) = array(GET_INDEX(1_IK))
        array(GET_INDEX(1_IK)) = pivot1
        array(GET_INDEX(2_IK * last / 3_IK)) = array(GET_INDEX(last))
        array(GET_INDEX(last)) = pivot2
        g = last
        l = 2_IK
        do while (IS_SORTED(array(GET_INDEX(l)), pivot1))
            l = l + 1_IK
        end do
        k = l
        do while(k < g)
            temp = array(GET_INDEX(k))
            if (IS_SORTED(temp, pivot1)) then
                array(GET_INDEX(k)) = array(GET_INDEX(l))
                array(GET_INDEX(l)) = temp
                l = l + 1_IK
            elseif (IS_SORTED(pivot2, temp)) then
                do while(IS_SORTED(pivot2, array(GET_INDEX(g - 1_IK))))
                    g = g - 1_IK
                end do
                if (k >= g) exit
                g = g - 1_IK
                if (IS_SORTED(array(GET_INDEX(g)), pivot1)) then
                    array(GET_INDEX(k)) = array(GET_INDEX(l))
                    array(GET_INDEX(l)) = array(GET_INDEX(g))
                    array(GET_INDEX(g)) = temp
                    l = l + 1_IK
                else
                    array(GET_INDEX(k)) = array(GET_INDEX(g))
                    array(GET_INDEX(g)) = temp
                end if
            end if
            k = k + 1_IK
        end do
        if (l > 2_IK) then
            array(GET_INDEX(1)) = array(GET_INDEX(l - 1_IK))
            array(GET_INDEX(l - 1_IK)) = pivot1
            call setSorted(array(1_IK : l - 2_IK)ISSORTED, qsortrdp) ! fpp
        end if
        call setSorted(array(l : g - 1_IK)ISSORTED, qsortrdp) ! fpp
        if (g < last) then
            array(GET_INDEX(last)) = array(GET_INDEX(g))
            array(GET_INDEX(g)) = pivot2
            call setSorted(array(g + 1_IK : last)ISSORTED, qsortrdp) ! fpp
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setSorted_ENABLED && Arr_ENABLED && Bubble_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_KIND :: temp
        integer(IK) :: i, j
        do i = GET_SIZE(array, kind = IK) - 1_IK, 1_IK, -1_IK
            do j = 1_IK, i
                if (IS_SORTED(array(GET_INDEX(j + 1_IK)), array(GET_INDEX(j)))) then
                    temp = array(GET_INDEX(j + 1_IK))
                    array(GET_INDEX(j + 1_IK)) = array(GET_INDEX(j))
                    array(GET_INDEX(j)) = temp
                end if
            end do
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setSorted_ENABLED && Arr_ENABLED && (Heapi_ENABLED || Heapr_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_KIND :: temp
        integer(IK) :: i, last
        last = GET_SIZE(array, kind = IK)
        ! Build heap
        do i = last / 2_IK, 1_IK, -1_IK
            call heapify(array, i)
        end do
        ! Unpick heap
        do i = last, 2_IK, -1_IK
            temp = array(GET_INDEX(1))
            array(GET_INDEX(1)) = array(GET_INDEX(i))
            array(GET_INDEX(i)) = temp
            call heapify(array(1:i-1), 1_IK)
        end do
    contains
#if     DefCom_ENABLED
        pure &
#endif
#if     Heapr_ENABLED
        recursive &
#endif
        subroutine heapify(array, i)
            integer(IK), intent(in) :: i
#if         SK_ENABLED && D0_ENABLED
            character(*,SKG), intent(inout):: array
#elif       SK_ENABLED && D1_ENABLED
            character(*,SKG), intent(inout), contiguous :: array(:)
#elif       D1_ENABLED
            TYPE_KIND, intent(inout) :: array(:)
#else
#error      "Unrecognized interface."
#endif
            TYPE_KIND :: temp
            integer(IK) :: left, right, root, last, largest
            last = GET_SIZE(array, kind = IK)
            root = i
            largest = root
            left = 2_IK * root
#if         Heapi_ENABLED
            temp = array(GET_INDEX(root))
            do while(left <= last)
                right = left + 1_IK
                if (left <= last) then
                    if(IS_SORTED(array(GET_INDEX(largest)), array(GET_INDEX(left)))) largest = left
                end if
                if (right <= last) then
                    if(IS_SORTED(array(GET_INDEX(largest)), array(GET_INDEX(right)))) largest = right
                end if
                if (largest == root) exit
                array(GET_INDEX(root)) = array(GET_INDEX(largest))
                array(GET_INDEX(largest)) = temp
                root = largest
                left = 2_IK * root
            end do
#elif       Heapr_ENABLED
            right = left + 1_IK
            if (left <= last) then
                if(IS_SORTED(array(GET_INDEX(largest)), array(GET_INDEX(left)))) largest = left
            end if
            if (right <= last) then
                if(IS_SORTED(array(GET_INDEX(largest)), array(GET_INDEX(right)))) largest = right
            end if
            if (largest /= root) then
                temp = array(GET_INDEX(root))
                array(GET_INDEX(root)) = array(GET_INDEX(largest))
                array(GET_INDEX(largest)) = temp
                call heapify(array, largest)
            end if
#else
#error      "Unrecognized interface."
#endif
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setSorted_ENABLED && Arr_ENABLED && Insertionl_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_KIND :: temp
        integer(IK) :: i, pos
        do i = 2_IK, GET_SIZE(array, kind = IK)
            temp = array(GET_INDEX(i))
            ! Do linear search
            do pos = i - 1_IK, 1_IK, -1_IK
                if (IS_SORTED(temp, array(GET_INDEX(pos)))) then
                    array(GET_INDEX(pos + 1_IK)) = array(GET_INDEX(pos))
                    cycle
                end if
                exit
            end do
            ! Insert
            array(GET_INDEX(pos + 1_IK)) = temp
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setSorted_ENABLED && Arr_ENABLED && Insertionb_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_KIND :: temp
        integer(IK) :: i, j, a, b, mid, pos
        do i = 2_IK, GET_SIZE(array, kind = IK)
            temp = array(GET_INDEX(i))
            ! Do binary search
            a = 1_IK
            b = i - 1_IK
            loopBinarySearch: do
                if (a == b) then
                    if (IS_SORTED(temp, array(GET_INDEX(a)))) then
                        pos = a
                    else
                        pos = a + 1_IK
                    end if
                    exit loopBinarySearch
                end if
                if (a > b) then
                    pos = a
                    exit loopBinarySearch
                end if
                mid = (a + b) / 2_IK
                if (IS_SORTED(array(GET_INDEX(mid)), temp)) then
                    a = mid + 1_IK
                elseif (IS_SORTED(temp, array(GET_INDEX(mid)))) then
                    b = mid - 1_IK
                else
                    pos = mid
                    exit loopBinarySearch
                end if
            end do loopBinarySearch
            do j = i, pos + 1_IK, -1_IK
                array(GET_INDEX(j)) = array(GET_INDEX(j - 1_IK))
            end do
            ! Insert
            array(GET_INDEX(pos)) = temp
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setSorted_ENABLED && Arr_ENABLED && Merger_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK_ENABLED && D0_ENABLED
        character(len(array, kind = IK),SKG) :: temp
#elif   D1_ENABLED
        TYPE_KIND :: temp(size(array, kind = IK))
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: last
        last = GET_SIZE(array, kind = IK)
        if (last < 1_IK) return
        if (last < 15_IK) then
            call setSorted(array ISSORTED, insertionl)
            return
        end if
        call setSorted(array(1 : last / 2)ISSORTED, merger) ! fpp
        call setSorted(array(last / 2 + 1 : last)ISSORTED, merger) ! fpp
        call setMerged(temp, array(1 : last / 2), array(last / 2 + 1 : last)ISSORTED) ! fpp
        array = temp

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setSorted_ENABLED && Arr_ENABLED && Selection_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_KIND :: temp
        integer(IK) :: i, j, k, lenArray
        lenArray = GET_SIZE(array, kind = IK)
        do i = 1, lenArray - 1_IK
            ! minloc
            temp = array(GET_INDEX(i))
            j = i
            do k = i + 1_IK, lenArray
                if (IS_SORTED(array(GET_INDEX(k)), temp)) then
                    temp = array(GET_INDEX(k))
                    j = k
                end if
            end do
            !j = minloc(array(i:), 1) + i - 1_IK
            temp = array(GET_INDEX(i))
            array(GET_INDEX(i)) = array(GET_INDEX(j))
            array(GET_INDEX(j)) = temp
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setSorted_ENABLED && Arr_ENABLED && Shell_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK_ENABLED && D0_ENABLED
        character(len(array, kind = IK),SKG) :: temp
#elif   D1_ENABLED
        TYPE_KIND :: temp(size(array, kind = IK))
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: i, j, counter, interval, lenArray, oneThirdLenArray
        interval = 1_IK
        lenArray = GET_SIZE(array, kind = IK)
        oneThirdLenArray = lenArray / 3_IK
        do while(interval < oneThirdLenArray)
            interval = 3_IK * interval + 1_IK
        end do
        do while(interval > 0_IK)
            do i = 1_IK, interval
                ! copy array(i), array(i+interval), array(i+2*interval), ...
                counter = 0_IK
                do j = i, lenArray, interval
                    counter = counter + 1_IK
                    temp(GET_INDEX(counter)) = array(GET_INDEX(j))
                end do
                call setSorted(temp(1:counter)ISSORTED, insertionl)
                counter = 0_IK
                do j = i, lenArray, interval
                    counter = counter + 1_IK
                    array(GET_INDEX(j)) = temp(GET_INDEX(counter))
                end do
            end do
            interval = (interval - 1_IK) / 3_IK
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   isSorted_ENABLED || isAscending_ENABLED || isAscendingAll_ENABLED || isDescending_ENABLED || isDescendingAll_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the sorting component.
#if     PSSK_ENABLED || BSSK_ENABLED
#define COMPONENT %val
#elif   CK_ENABLED
#define COMPONENT %re
#else
#define COMPONENT
#endif
        ! Define the comparison operator.
#if     DefCom_ENABLED && isAscending_ENABLED
#define NOT_COMPARES_WITH >
#define SORTED ascending
#elif   DefCom_ENABLED && isDescending_ENABLED
#define NOT_COMPARES_WITH <
#define SORTED descending
#elif   DefCom_ENABLED && isAscendingAll_ENABLED
#define NOT_COMPARES_WITH >=
#define SORTED ascendingAll
#elif   DefCom_ENABLED && isDescendingAll_ENABLED
#define NOT_COMPARES_WITH <=
#define SORTED descendingAll
#elif   !isSorted_ENABLED
#error  "Unrecognized interface."
#endif
        ! Define the logical comparison operator.
#if     LK_ENABLED && isAscending_ENABLED
        use pm_logicalCompare, only: operator(>)
#elif   LK_ENABLED && isDescending_ENABLED
        use pm_logicalCompare, only: operator(<)
#elif   LK_ENABLED && isAscendingAll_ENABLED
        use pm_logicalCompare, only: operator(>=)
#elif   LK_ENABLED && isDescendingAll_ENABLED
        use pm_logicalCompare, only: operator(<=)
#elif   LK_ENABLED && DefCom_ENABLED && isSorted_ENABLED
        use pm_logicalCompare, only: operator(<=), operator(>=)
#endif
        integer(IK) :: i
#if     DefCom_ENABLED && (isAscending_ENABLED || isAscendingAll_ENABLED || isDescending_ENABLED || isDescendingAll_ENABLED)
        SORTED = .false._LK
        do i = 1_IK, GET_SIZE(array, kind = IK) - 1_IK
            if (array(GET_INDEX(i))COMPONENT NOT_COMPARES_WITH array(GET_INDEX(i + 1))COMPONENT) return
        end do
        SORTED = .true._LK
#elif   DefCom_ENABLED && isSorted_ENABLED
        logical(LK) :: isAscending, isDescending
        sorted = .true._LK
        isAscending = .true._LK
        isDescending = .true._LK
        do i = 1_IK, GET_SIZE(array, kind = IK) - 1_IK
            isAscending  = isAscending  .and. array(GET_INDEX(i))COMPONENT <= array(GET_INDEX(i+1_IK))COMPONENT
            isDescending = isDescending .and. array(GET_INDEX(i))COMPONENT >= array(GET_INDEX(i+1_IK))COMPONENT
            sorted = isAscending .or. isDescending
            if (.not. sorted) exit
        end do
#elif   CusCom_ENABLED && isSorted_ENABLED
        sorted = .true._LK
        do i = 1_IK, GET_SIZE(array, kind = IK) - 1_IK
            sorted = sorted .and. isSorted(array(GET_INDEX(i))COMPONENT, array(GET_INDEX(i+1_IK))COMPONENT)
            if (.not. sorted) exit
        end do
#else
#error  "Unrecognized interface."
#endif
#undef  NOT_COMPARES_WITH
#undef  COMPONENT
#undef  SORTED

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef TYPE_KIND_LEN
#undef TYPE_KIND
#undef IS_SORTED
#undef GET_INDEX
#undef GET_SIZE
#undef ISSORTED