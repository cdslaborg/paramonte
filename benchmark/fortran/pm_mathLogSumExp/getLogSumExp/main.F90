! Test the performance of `getLogSumExp()` with and without underflow check vs. direct computation of `logSumExp`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_kind, only: IK, LK, RK, SK
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: i                !<  The procedure benchmark counter.
    integer(IK)                         :: iarr             !<  The array size counter.
    integer(IK)                         :: fileUnitN        !<  The output file unit for benchmark results when the input array is normal.
    integer(IK)                         :: fileUnitU        !<  The output file unit for benchmark results when the input array has many cases of underflow.
    integer(IK) , parameter             :: NARR = 11_IK     !<  The number of benchmark array sizes.
    integer(IK) , parameter             :: NBENCH = 3_IK    !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(NARR)  !<  The sizes of the benchmark array.
    real(RK)                            :: dummySum = 0._RK !<  The dummy computation to prevent the compiler from aggressive optimizations.
    real(RK)                            :: maxArray         !<  The maximum value of the benchmark array.
    real(RK)    , allocatable           :: array(:)         !<  The benchmark array.
    real(RK)                            :: logSumExp        !<  The logarithm of the sum of the exponential of the `array`.
    type(bench_type)                    :: bench(2,NBENCH)  !<  The Benchmark array.
    logical(LK)                         :: underflowEnabled !<  The logical flag indicating whether an array with many instances of underflow should be generated.

    bench(1:2,1) = bench_type(name = SK_"getLogSumExpMaxSequence", exec = getLogSumExpMaxSequence , overhead = setOverhead)
    bench(1:2,2) = bench_type(name = SK_"getLogSumExpMaxSelection", exec = getLogSumExpMaxSelection , overhead = setOverhead)
    bench(1:2,3) = bench_type(name = SK_"getLogSumExpDirect", exec = getLogSumExpDirect, overhead = setOverhead)

    arraySize = [( 2_IK**iarr, iarr = 1_IK, NARR )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "getLogSumExp(...) vs. getLogSumExp(..., control = selection) vs. direct method."
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnitN, file = "main.normal.out", status = "replace")
    open(newunit = fileUnitU, file = "main.underflow.out", status = "replace")

        write(fileUnitN, "(*(g0,:,','))") "arraySize", (bench(1,i)%name, i = 1, NBENCH)
        write(fileUnitU, "(*(g0,:,','))") "arraySize", (bench(2,i)%name, i = 1, NBENCH)

        loopOverArraySize: do iarr = 1, NARR

            allocate(array(arraySize(iarr)))
            write(*,"(*(g0,:,' '))") "Benchmarking with array size", arraySize(iarr)

            underflowEnabled = .false._LK
            do i = 1, NBENCH
            !do i = NBENCH, 1, -1
                bench(1,i)%timing = bench(1,i)%getTiming(minsec = 0.07_RK)
            end do
            write(fileUnitN,"(*(g0,:,','))") arraySize(iarr), (bench(1,i)%timing%mean, i = 1, NBENCH)

            underflowEnabled = .true._LK
            do i = 1, NBENCH
                bench(2,i)%timing = bench(2,i)%getTiming(minsec = 0.07_RK)
            end do
            write(fileUnitU,"(*(g0,:,','))") arraySize(iarr), (bench(2,i)%timing%mean, i = 1, NBENCH)

            deallocate(array)

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
        call random_number(array)
        if (underflowEnabled) array = array * (maxexponent(0._RK) - minexponent(0._RK)) + minexponent(0._RK)
        maxArray = maxval(array)
    end subroutine

    subroutine getDummy()
        dummySum = dummySum + logSumExp
    end subroutine

    subroutine getLogSumExpMaxSequence()
        use pm_mathLogSumExp, only: getLogSumExp
        call setArray()
        logSumExp = getLogSumExp(array, maxArray)
        call getDummy()
    end subroutine

    subroutine getLogSumExpMaxSelection()
        use pm_mathLogSumExp, only: getLogSumExp, selection
        call setArray()
        logSumExp = getLogSumExp(array, maxArray, selection)
        call getDummy()
    end subroutine

    subroutine getLogSumExpDirect()
        call setArray()
        logSumExp = maxArray + log(sum(exp(array - maxArray)))
        call getDummy()
    end subroutine

end program benchmark