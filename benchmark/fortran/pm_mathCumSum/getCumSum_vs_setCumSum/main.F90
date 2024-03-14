! Test the performance of `getCumSum()` vs. `setCumSum()`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_kind, only: IK, RK, SK
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: i                !<  The procedure benchmark counter.
    integer(IK)                         :: iarr             !<  The array size counter.
    integer(IK)                         :: fileUnit         !<  The output file unit for benchmark results.
    integer(IK)     , parameter         :: NARR = 11_IK     !<  The number of benchmark array sizes.
    integer(IK)                         :: arraySize(NARR)  !<  The sizes of the benchmark array.
    real(RK)                            :: dummySum = 0._RK !<  The dummy computation to prevent the compiler from aggressive optimizations.
    real(RK)        , allocatable       :: array(:)         !<  The benchmark array.
    real(RK)        , allocatable       :: cumsum(:)        !<  The cumulative proportion of the exponential of the `array`.
    type(bench_type), allocatable       :: bench(:)         !<  The Benchmark array.

    bench = [ bench_type(name = SK_"setCumSum", exec = setCumSum, overhead = setOverhead) &
            , bench_type(name = SK_"getCumSum", exec = getCumSum, overhead = setOverhead) &
            , bench_type(name = SK_"setCumSum_overwrite", exec = setCumSum_overwrite, overhead = setOverhead) &
            , bench_type(name = SK_"getCumSum_withBounds", exec = getCumSum_withBounds, overhead = setOverhead) &
            ]

    arraySize = [( 2_IK**iarr, iarr = 1_IK, NARR )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "setCumSum() vs. getCumSum() vs. getCumSum_withBounds()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1, size(bench))

        loopOverArraySize: do iarr = 1, NARR

            allocate(array(arraySize(iarr)))
            allocate(cumsum(arraySize(iarr)), source = 0._RK)
            write(*,"(*(g0,:,' '))") "Benchmarking with array size", arraySize(iarr)

            do i = 1, size(bench)
                bench(i)%timing = bench(i)%getTiming() !, minsec = 0.1_RK)
            end do
            write(fileUnit,"(*(g0,:,','))") arraySize(iarr), (bench(i)%timing%mean, i = 1, size(bench))

            deallocate(array, cumsum)

        end do loopOverArraySize
        write(*,"(*(g0,:,' '))") dummySum
        write(*,"(*(g0,:,' '))")

    close(fileUnit)

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
    end subroutine

    subroutine getDummy()
        dummySum = dummySum + cumsum(1) + array(1)
    end subroutine

    subroutine setCumSum()
        block
            use pm_mathCumSum, only: setCumSum
            call setArray()
            call setCumSum(cumsum, array)
            call getDummy()
        end block
    end subroutine

    subroutine setCumSum_overwrite()
        block
            use pm_mathCumSum, only: setCumSum
            call setArray()
            call setCumSum(array)
            call getDummy()
        end block
    end subroutine

    subroutine getCumSum()
        block
            use pm_mathCumSum, only: getCumSum
            call setArray()
            cumsum = getCumSum(array)
            call getDummy()
        end block
    end subroutine

    subroutine getCumSum_withBounds()
        block
            use pm_mathCumSum, only: getCumSum
            call setArray()
            cumsum(:) = getCumSum(array)
            call getDummy()
        end block
    end subroutine

end program benchmark