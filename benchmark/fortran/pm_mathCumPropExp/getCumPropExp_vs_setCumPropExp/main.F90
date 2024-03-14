! Test the performance of `getCumPropExp()` vs. `setCumPropExp()`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_bench, only: bench_type
    use pm_kind, only: IK, RK, SK

    implicit none

    integer(IK)                         :: i                !<  The procedure benchmark counter.
    integer(IK)                         :: iarr             !<  The array size counter.
    integer(IK)                         :: fileUnit         !<  The output file unit for benchmark results.
    integer(IK)     , parameter         :: NARR = 11_IK     !<  The number of benchmark array sizes.
    integer(IK)                         :: arraySize(NARR)  !<  The sizes of the benchmark array.
    real(RK)                            :: dummySum = 0._RK !<  The dummy computation to prevent the compiler from aggressive optimizations.
    real(RK)                            :: maxArray         !<  The maximum value of the benchmark array.
    real(RK)        , allocatable       :: array(:)         !<  The benchmark array.
    real(RK)        , allocatable       :: cumPropExp(:)    !<  The cumulative proportion of the exponential of the `array`.
    type(bench_type), allocatable       :: bench(:)         !<  The Benchmark array.

    bench = [ bench_type(name = SK_"setCumPropExp", exec = setCumPropExp , overhead = setOverhead) &
            , bench_type(name = SK_"getCumPropExp", exec = getCumPropExp , overhead = setOverhead) &
            ]

    arraySize = [( 2_IK**iarr, iarr = 1_IK, NARR )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "getCumPropExp() vs. setCumPropExp()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1, size(bench))

        loopOverArraySize: do iarr = 1, NARR

            allocate(array(arraySize(iarr)))
            allocate(cumPropExp(arraySize(iarr)), source = 0._RK)
            write(*,"(*(g0,:,' '))") "Benchmarking with array size", arraySize(iarr)

            do i = 1, size(bench)
                bench(i)%timing = bench(i)%getTiming(minsec = 0.1_RK)
            end do
            write(fileUnit,"(*(g0,:,','))") arraySize(iarr), (bench(i)%timing%mean, i = 1, size(bench))

            deallocate(array, cumPropExp)

        end do loopOverArraySize
        write(*,"(*(g0,:,' '))") dummySum
        write(*,"(*(g0,:,' '))")

    close(fileUnit)

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
        maxArray = maxval(array)
    end subroutine

    subroutine getDummy()
        dummySum = dummySum + cumPropExp(1)
    end subroutine

    subroutine getCumPropExp()
        block
            use pm_mathCumPropExp, only: getCumPropExp
            call getArray()
            cumPropExp = getCumPropExp(array, maxArray)
            call getDummy()
        end block
    end subroutine

    subroutine setCumPropExp()
        block
            use pm_mathCumPropExp, only: setCumPropExp, sequence
            call getArray()
            call setCumPropExp(cumPropExp, array, maxArray, sequence)
            call getDummy()
        end block
    end subroutine

end program benchmark