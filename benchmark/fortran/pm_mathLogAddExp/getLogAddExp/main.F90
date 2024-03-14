! Test the performance of `getLogAddExp()` with and without the `control` argument.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_bench, only: bench_type
    use pm_kind, only: IK, LK, RK, SK

    implicit none

    integer(IK)                         :: i                !<  The procedure benchmark counter.
    integer(IK)                         :: iarr             !<  The array size counter.
    integer(IK)                         :: fileUnitN        !<  The output file unit for benchmark results when the input array is normal.
    integer(IK)                         :: fileUnitU        !<  The output file unit for benchmark results when the input array has many cases of underflow.
    integer(IK) , parameter             :: NARR = 11_IK     !<  The number of benchmark array sizes.
    integer(IK) , parameter             :: NBENCH = 2_IK    !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(NARR)  !<  The sizes of the benchmark array.
    real(RK)                            :: dummySum = 0._RK !<  The dummy computation to prevent the compiler from aggressive optimizations.
    real(RK)                            :: maxArray         !<  The maximum value of the benchmark array.
    real(RK)    , allocatable           :: Larger(:)        !<  The benchmark array.
    real(RK)    , allocatable           :: Smaller(:)       !<  The benchmark array.
    real(RK)    , allocatable           :: LogAddExp(:)     !<  The result.
    type(bench_type)                    :: bench(2,NBENCH)  !<  The Benchmark array.
    logical(LK)                         :: underflowEnabled !<  The logical flag indicating whether an array with many instances of underflow should be generated.

    bench(1:2,1) = bench_type(name = SK_"getLogAddExpSequence", exec = getLogAddExpSequence , overhead = setOverhead)
    bench(1:2,2) = bench_type(name = SK_"getLogAddExpSelection", exec = getLogAddExpSelection , overhead = setOverhead)

    arraySize = [( 2_IK**iarr, iarr = 1_IK, NARR )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "getLogAddExp(...) vs. getLogAddExp(..., control) vs. direct method."
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnitN, file = "main.normal.out", status = "replace")
    open(newunit = fileUnitU, file = "main.underflow.out", status = "replace")

        write(fileUnitN, "(*(g0,:,','))") "arraySize", (bench(1,i)%name, i = 1, NBENCH)
        write(fileUnitU, "(*(g0,:,','))") "arraySize", (bench(2,i)%name, i = 1, NBENCH)

        loopOverArraySize: do iarr = 1, NARR

            allocate(Larger(arraySize(iarr)))
            allocate(Smaller(arraySize(iarr)))
            allocate(LogAddExp(arraySize(iarr)))
            write(*,"(*(g0,:,' '))") "Benchmarking with array size", arraySize(iarr)

            underflowEnabled = .false._LK
            do i = 1, NBENCH
                bench(1,i)%timing = bench(1,i)%getTiming(minsec = 0.07_RK)
            end do
            write(fileUnitN,"(*(g0,:,','))") arraySize(iarr), (bench(1,i)%timing%mean, i = 1, NBENCH)

            underflowEnabled = .true._LK
            do i = 1, NBENCH
                bench(2,i)%timing = bench(2,i)%getTiming(minsec = 0.07_RK)
            end do
            write(fileUnitU,"(*(g0,:,','))") arraySize(iarr), (bench(2,i)%timing%mean, i = 1, NBENCH)

            deallocate(LogAddExp, Smaller, Larger)

        end do loopOverArraySize
        write(*,"(*(g0,:,' '))") dummySum
        write(*,"(*(g0,:,' '))")

    close(fileUnitN)
    close(fileUnitU)

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! procedure wrappers.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setOverhead()
        call setArray()
        call getDummy()
    end subroutine

    subroutine setArray()
        use pm_mathMinMax, only: setMinMax
        call random_number(Larger)
        call random_number(Smaller)
        if (underflowEnabled) then
            Larger = Larger * (maxexponent(0._RK) - minexponent(0._RK)) + minexponent(0._RK)
            Smaller = Smaller * (maxexponent(0._RK) - minexponent(0._RK)) + minexponent(0._RK)
        end if
        call setMinMax(Smaller, Larger)
    end subroutine

    subroutine getDummy()
        dummySum = dummySum + LogAddExp(1)
    end subroutine

    subroutine getLogAddExpSequence()
        use pm_mathLogAddExp, only: getLogAddExp
        call setArray()
        LogAddExp(:) = getLogAddExp(Smaller, Larger)
        call getDummy()
    end subroutine

    subroutine getLogAddExpSelection()
        use pm_mathLogAddExp, only: getLogAddExp, selection
        call setArray()
        LogAddExp(:) = getLogAddExp(Smaller, Larger, selection)
        call getDummy()
    end subroutine

end program benchmark