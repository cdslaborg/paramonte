! Test the performance of different sorting algorithms in pm_arraySort.
program benchmark

    use pm_kind, only: IK, LK, RK, SK
    use pm_bench, only: bench_type

    implicit none

    logical(LK)                         :: issorted                 !<  The logical flag indicating whether the sorting is being on pre-sorted arrays.
    integer(IK)                         :: i                        !<  The sorting procedure benchmark counter.
    integer(IK)                         :: iarr                     !<  The array size counter.
    integer(IK)                         :: fileUnitRandom           !<  The output file unit for sorting benchmark with random input data.
    integer(IK)                         :: fileUnitSorted           !<  The output file unit for sorted sorting benchmark with sorted input data.
    integer(IK) , parameter             :: NARR = 12_IK             !<  The number of benchmark array sizes.
    integer(IK) , parameter             :: NBENCH = 11_IK           !<  The number of sorting procedures.
    integer(IK)                         :: arraySize(NARR)          !<  The benchmarking array sizes.
    real(RK)  , allocatable             :: array(:)                 !<  The arrays to be sorted.
    type(bench_type)                    :: bench(2,NBENCH)          !<  The Benchmark array.
                                                                    !!  The first dimension holds benchmarks of a given sorting algorithm for random and sorted array.
                                                                    !!  The second dimension holds benchmarks of different sorting algorithms.

    bench(1:2, 1) = bench_type(name = SK_"setSortedQsorti       ", exec = setSortedQsorti       , overhead = setOverhead)
    bench(1:2, 2) = bench_type(name = SK_"setSortedQsortr       ", exec = setSortedQsortr       , overhead = setOverhead)
    bench(1:2, 4) = bench_type(name = SK_"setSortedQsortrdp     ", exec = setSortedQsortrdp     , overhead = setOverhead)
    bench(1:2, 3) = bench_type(name = SK_"setSortedBubble       ", exec = setSortedBubble       , overhead = setOverhead)
    bench(1:2, 5) = bench_type(name = SK_"setSortedHeapi        ", exec = setSortedHeapi        , overhead = setOverhead)
    bench(1:2, 6) = bench_type(name = SK_"setSortedHeapr        ", exec = setSortedHeapr        , overhead = setOverhead)
    bench(1:2, 7) = bench_type(name = SK_"setSortedInsertionl   ", exec = setSortedInsertionl   , overhead = setOverhead)
    bench(1:2, 8) = bench_type(name = SK_"setSortedInsertionb   ", exec = setSortedInsertionb   , overhead = setOverhead)
    bench(1:2, 9) = bench_type(name = SK_"setSortedMerge        ", exec = setSortedMerge        , overhead = setOverhead)
    bench(1:2,10) = bench_type(name = SK_"setSortedSelection    ", exec = setSortedSelection    , overhead = setOverhead)
    bench(1:2,11) = bench_type(name = SK_"setSortedShell        ", exec = setSortedShell        , overhead = setOverhead)

    arraySize = [( 2_IK**iarr, iarr = 3_IK, NARR + 2_IK )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "setSorted() vs. others."
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnitRandom, file = "main.random.out", status = "replace")
    open(newunit = fileUnitSorted, file = "main.sorted.out", status = "replace")

        write(fileUnitRandom       ,"(*(g0,:,','))") "arraySize", (bench(1,i)%name, i = 1, NBENCH)
        write(fileUnitSorted,"(*(g0,:,','))") "arraySize", (bench(2,i)%name, i = 1, NBENCH)

        loopOverArraySize: do iarr = 1, NARR

            allocate(array(arraySize(iarr)))
            write(*,"(*(g0,:,' '))") "Benchmarking sorting algorithms with array size", arraySize(iarr)

            ! Time sorting algorithms against the default setSorted() procedures.

            issorted = .false._LK
            do i = 1, NBENCH
                bench(1,i)%timing = bench(1,i)%getTiming()
            end do
            write(fileUnitRandom,"(*(g0,:,','))") arraySize(iarr), (bench(1,i)%timing%mean, i = 1, NBENCH)

            ! Time sorting algorithms with sorted input arrays against unordered input arrays.

            issorted = .true._LK
            do i = 1, NBENCH
                bench(2,i)%timing = bench(2,i)%getTiming()
            end do
            write(fileUnitSorted,"(*(g0,:,','))") arraySize(iarr), (bench(2,i)%timing%mean, i = 1, NBENCH)

            deallocate(array)

        end do loopOverArraySize
        write(*,"(*(g0,:,' '))")

    close(fileUnitSorted)
    close(fileUnitRandom)

contains

    ! sorting procedure wrappers.

    subroutine setOverhead()
        call setArray()
    end subroutine

    subroutine setArray()
        use pm_arraySpace, only: getLinSpace
        if (issorted) then
            array(:) = getLinSpace(0._RK, 1._RK, count = size(array, kind = IK))
        else
            call random_number(array)
        end if
    end subroutine

    subroutine setSortedQsorti()
        call setArray()
        block; use pm_arraySort, only: setSorted, qsorti; call setSorted(array, qsorti); end block
    end subroutine

    subroutine setSortedQsortr()
        call setArray()
        block; use pm_arraySort, only: setSorted, qsortr; call setSorted(array, qsortr); end block
    end subroutine

    subroutine setSortedQsortrdp()
        call setArray()
        block; use pm_arraySort, only: setSorted, qsortrdp; call setSorted(array, qsortrdp); end block
    end subroutine

    subroutine setSortedBubble()
        call setArray()
        block; use pm_arraySort, only: setSorted, bubble; call setSorted(array, bubble); end block
    end subroutine

    subroutine setSortedHeapi()
        call setArray()
        block; use pm_arraySort, only: setSorted, heapi; call setSorted(array, heapi); end block
    end subroutine

    subroutine setSortedHeapr()
        call setArray()
        block; use pm_arraySort, only: setSorted, heapr; call setSorted(array, heapr); end block
    end subroutine

    subroutine setSortedInsertionl()
        call setArray()
        block; use pm_arraySort, only: setSorted, insertionl; call setSorted(array, insertionl); end block
    end subroutine

    subroutine setSortedInsertionb()
        call setArray()
        block; use pm_arraySort, only: setSorted, insertionb; call setSorted(array, insertionb); end block
    end subroutine

    subroutine setSortedMerge()
        call setArray()
        block; use pm_arraySort, only: setSorted, merger; call setSorted(array, merger); end block
    end subroutine

    subroutine setSortedSelection()
        call setArray()
        block; use pm_arraySort, only: setSorted, selection; call setSorted(array, selection); end block
    end subroutine

    subroutine setSortedShell()
        call setArray()
        block; use pm_arraySort, only: setSorted, shell; call setSorted(array, shell); end block
    end subroutine

end program benchmark