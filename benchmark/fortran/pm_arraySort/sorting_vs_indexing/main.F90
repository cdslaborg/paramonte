program benchmark

    use pm_kind, only: IK, LK, RKG => RK
    use pm_arrayResize, only: setResized
    use pm_distUnif, only: getUnifRand
    use pm_arraySort, only: setSorted
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: i,iarr                   !< The counters.
    integer(IK)                         :: fileUnitRandom           !< The output file unit for random data benchmark.
    integer(IK)                         :: fileUnitSorted           !< The output file unit for ordered data benchmark.
    integer(IK) , parameter             :: NARR = 12_IK             !< The number of benchmark array sizes.
    integer(IK) , parameter             :: NBENCH = 3_IK            !< The number of benchmark wrappers.
    integer(IK)                         :: arraySize(NARR)          !< The benchmarking array sizes.
    real(RKG)   , allocatable           :: array(:)                 !< The arrays to be sorted.
    integer(IK) , allocatable           :: index(:)                 !< The arrays to be sorted.
    logical(LK)                         :: issorted
    type(bench_type)                    :: bench(2,NBENCH)

    bench(:,1) = bench_type("setSortedArray", setSortedArray, setOverhead)
    bench(:,2) = bench_type("setSortedIndex", setSortedIndex, setOverhead)
    bench(:,3) = bench_type("setSortedindexWithSorting", setSortedindexWithSorting, setOverhead)

    arraySize = [( 2_IK**i, i = 3_IK, NARR + 2_IK )]

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Time the setSortedArray() against setSortedIndex() for different array sizes.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "setSortedArray() vs. setSortedIndex()."
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnitRandom, file = "main.random.out", status = "replace")
    open(newunit = fileUnitSorted, file = "main.sorted.out", status = "replace")
    write(fileUnitRandom,"(*(g0,:,','))") "arraySize", (bench(1,i)%name, i = 1, NBENCH)
    write(fileUnitSorted,"(*(g0,:,','))") "arraySize", (bench(2,i)%name, i = 1, NBENCH)

    loopOverarraySize: do iarr = 1, NARR

        write(*,"(*(g0,:,' '))") "Benchmarking with array size", arraySize(iarr)

        array = getUnifRand(1._RKG, 2._RKG, arraySize(iarr))
        call setResized(index, arraySize(iarr))

        issorted = .false._LK
        do i = 1, NBENCH
            bench(1,i)%timing = bench(1,i)%getTiming()
        end do

        write(fileUnitRandom,"(*(g0,:,','))") arraySize(iarr), (bench(1,i)%timing%mean, i = 1, NBENCH)

        issorted = .true._LK
        do i = 1, NBENCH
            bench(2,i)%timing = bench(2,i)%getTiming()
        end do

        write(fileUnitSorted,"(*(g0,:,','))") arraySize(iarr), (bench(2,i)%timing%mean, i = 1, NBENCH)

    end do loopOverarraySize

    close(fileUnitRandom)
    close(fileUnitSorted)

contains

    subroutine setarray()
        use pm_arraySpace, only: getLinSpace
        if (issorted) then
            array(:) = getLinSpace(0._RKG, 1._RKG, count = size(array, kind = IK))
        else
            call random_number(array)
        end if
    end subroutine

    subroutine setOverhead()
        call setarray()
    end subroutine

    subroutine setSortedArray()
        call setSorted(array)
    end subroutine

    subroutine setSortedIndex()
        call setSorted(array, index)
    end subroutine

    subroutine setSortedindexWithSorting()
        call setSorted(array, index)
        array(:) = array(index)
    end subroutine

end program benchmark