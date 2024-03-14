! Test the performance of `getReversed()` vs. `setReversed()`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_bench, only: bench_type
    use pm_kind, only: IK, RK, SK

    implicit none

    integer(IK)                         :: i                            !<  The procedure benchmark counter.
    integer(IK)                         :: iSIZE                        !<  The Array size counter.
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK) , parameter             :: NSIZE = 11_IK                !<  The number of benchmark array sizes.
    integer(IK) , parameter             :: NBENCH = 3_IK                !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(NSIZE)             !<  The sizes of the benchmark Array.
    real(RK)                            :: dummySum = 0._RK             !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    real(RK)    , allocatable           :: Array(:), ArrayReverse(:)    !<  The benchmark arrays.
    type(bench_type)                    :: bench(NBENCH)                !<  The Benchmark object.

    bench(1) = bench_type(name = SK_"getReversed", exec = getReversed, overhead = setOverhead)
    bench(2) = bench_type(name = SK_"setReversed", exec = setReversed, overhead = setOverhead)
    bench(3) = bench_type(name = SK_"setReversed_overwrite", exec = setReversed_overwrite , overhead = setOverhead)

    arraySize = [( 2_IK**iSIZE, iSIZE = 1_IK, NSIZE )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "getReversed() vs. setReversed() vs. setReversed_overwrite()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1, NBENCH)

        loopOverArraySize: do iSIZE = 1, NSIZE

            allocate(Array(arraySize(iSIZE)), ArrayReverse(arraySize(iSIZE)))
            write(*,"(*(g0,:,' '))") "Benchmarking with rank", arraySize(iSIZE)

            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming(minsec = 0.05_RK)
            end do

            write(fileUnit,"(*(g0,:,','))") arraySize(iSIZE), (bench(i)%timing%mean, i = 1, NBENCH)
            deallocate(Array, ArrayReverse)

        end do loopOverArraySize
        write(*,"(*(g0,:,' '))") dummySum
        write(*,"(*(g0,:,' '))")

    close(fileUnit)

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! procedure wrappers.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setOverhead()
        call getDummy()
    end subroutine

    subroutine getDummy()
        dummySum = dummySum + sum(Array)
    end subroutine

    subroutine getReversed()
        block
            use pm_arrayReverse, only: getReversed
            ArrayReverse = getReversed(Array)
            call getDummy()
        end block
    end subroutine

    subroutine setReversed()
        block
            use pm_arrayReverse, only: setReversed
            call setReversed(Array, ArrayReverse)
            call getDummy()
        end block
    end subroutine

    subroutine setReversed_overwrite()
        block
            use pm_arrayReverse, only: setReversed
            call setReversed(Array)
            call getDummy()
        end block
    end subroutine

end program benchmark