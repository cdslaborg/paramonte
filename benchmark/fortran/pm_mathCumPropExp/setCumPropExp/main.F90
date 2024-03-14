! Test the performance of `setCumPropExp()` with and without the `cenabled` argument.
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
    integer(IK) , parameter             :: NBENCH = 3_IK    !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(NARR)  !<  The sizes of the benchmark array.
    real(RK)                            :: dummySum = 0._RK !<  The dummy computation to prevent the compiler from aggressive optimizations.
    real(RK)                            :: maxArray         !<  The maximum value of the benchmark array.
    real(RK)    , allocatable           :: array(:)         !<  The benchmark array.
    real(RK)    , allocatable           :: cumPropExp(:)    !<  The cumulative proportion of the exponential of the `array`.
    type(bench_type)                    :: bench(2,NBENCH)  !<  The Benchmark array.
    logical(LK)                         :: underflowEnabled !<  The logical flag indicating whether an array with many instances of underflow should be generated.

    bench(1:2,1) = bench_type(name = SK_"setCumPropExpSequence", exec = setCumPropExpSequence, overhead = setOverhead)
    bench(1:2,2) = bench_type(name = SK_"setCumPropExpSelection", exec = setCumPropExpSelection, overhead = setOverhead)
    bench(1:2,3) = bench_type(name = SK_"setCumPropExpDirect", exec = setCumPropExpDirect, overhead = setOverhead)

    arraySize = [( 2_IK**iarr, iarr = 1_IK, NARR )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "setCumPropExp(..., sequence) vs. setCumPropExp(..., selection) vs. direct method."
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnitN, file = "main.normal.out", status = "replace")
    open(newunit = fileUnitU, file = "main.underflow.out", status = "replace")

        write(fileUnitN, "(*(g0,:,','))") "arraySize", (bench(1,i)%name, i = 1, NBENCH)
        write(fileUnitU, "(*(g0,:,','))") "arraySize", (bench(2,i)%name, i = 1, NBENCH)

        loopOverArraySize: do iarr = 1, NARR

            allocate(array(arraySize(iarr)))
            allocate(cumPropExp(arraySize(iarr)), source = 0._RK)
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

            deallocate(array, cumPropExp)

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
        call getArray()
        call getDummy()
    end subroutine

    subroutine getArray()
        call random_number(array)
        if (underflowEnabled) array = array * (maxexponent(0._RK) - minexponent(0._RK)) + minexponent(0._RK)
        maxArray = maxval(array)
    end subroutine

    subroutine getDummy()
        dummySum = dummySum + cumPropExp(1) + array(1)
    end subroutine

    subroutine setCumPropExpSequence()
        use pm_mathCumPropExp, only: setCumPropExp, sequence
        call getArray()
        call setCumPropExp(cumPropExp, array, maxArray, sequence)
        call getDummy()
    end subroutine

    subroutine setCumPropExpSelection()
        use pm_mathCumPropExp, only: setCumPropExp, selection
        call getArray()
        call setCumPropExp(cumPropExp, array, maxArray, selection)
        call getDummy()
    end subroutine

    subroutine setCumPropExpDirect()
        use pm_mathCumSum, only: setCumSum
        call getArray()
        call setCumSum(cumPropExp, exp(array - maxArray))
        cumPropExp = cumPropExp / cumPropExp(size(cumPropExp))
        call getDummy()
    end subroutine

end program benchmark