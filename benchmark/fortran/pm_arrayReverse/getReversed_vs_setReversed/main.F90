! Test the performance of `getReversed()` vs. `setReversed()`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_arrayResize, only: setResized
    use pm_distUnif, only: getUnifRand
    use pm_bench, only: bench_type
    use pm_kind, only: IK, RK, SK

    implicit none

    integer(IK)                         :: i                            !<  The procedure benchmark counter.
    integer(IK)                         :: isize                        !<  The array size counter.
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK) , parameter             :: nsize = 11_IK                !<  The number of benchmark array sizes.
    integer(IK) , parameter             :: NBENCH = 3_IK                !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(nsize)             !<  The sizes of the benchmark array.
    real(RK)                            :: dummySum = 0._RK             !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    real(RK)    , allocatable           :: array(:), arrayReverse(:)    !<  The benchmark arrays.
    type(bench_type)                    :: bench(NBENCH)                !<  The Benchmark object.

    bench(1) = bench_type(name = SK_"getReversed", exec = getReversed, overhead = setOverhead)
    bench(2) = bench_type(name = SK_"setReversed", exec = setReversed, overhead = setOverhead)
    bench(3) = bench_type(name = SK_"setReversed_overwrite", exec = setReversed_overwrite , overhead = setOverhead)

    arraySize = [( 2_IK**isize, isize = 1_IK, nsize )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "getReversed() vs. setReversed() vs. setReversed_overwrite()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1, NBENCH)

        loopOverArraySize: do isize = 1, nsize

            array = getUnifRand(1_IK, 9_IK, arraySize(isize))
            call setResized(arrayReverse, arraySize(isize))
            write(*,"(*(g0,:,' '))") "Benchmarking with rank", arraySize(isize)

            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming(minsec = 0.05_RK)
            end do

            write(fileUnit,"(*(g0,:,','))") arraySize(isize), (bench(i)%timing%mean, i = 1, NBENCH)

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
        dummySum = dummySum + sum(array)
    end subroutine

    subroutine getReversed()
        block
            use pm_arrayReverse, only: getReversed
            arrayReverse = getReversed(array)
            call getDummy()
        end block
    end subroutine

    subroutine setReversed()
        block
            use pm_arrayReverse, only: setReversed
            call setReversed(array, arrayReverse)
            call getDummy()
        end block
    end subroutine

    subroutine setReversed_overwrite()
        block
            use pm_arrayReverse, only: setReversed
            call setReversed(array)
            call getDummy()
        end block
    end subroutine

end program benchmark