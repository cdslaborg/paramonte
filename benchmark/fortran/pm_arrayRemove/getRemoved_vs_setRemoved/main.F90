! Test the performance of `getRemoved()` vs. `setRemoved()`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_kind, only: IK, LK, RK, SK
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: i                            !<  The procedure benchmark counter.
    integer(IK)                         :: isize                        !<  The array size counter.
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK) , parameter             :: NSIZE = 12_IK                !<  The number of benchmark ranks.
    integer(IK)                         :: arraySize(NSIZE)             !<  The sizes of the benchmark array.
    logical(LK)                         :: dummy = .true._LK            !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    character(:, SK), allocatable       :: array                        !<  The array whose elements will have to be removed.
    character(*, SK), parameter         :: pattern = "a"                !<  The pattern to be removed from array.
    type(bench_type), allocatable       :: bench(:)                     !<  The Benchmark array.

    bench = [ bench_type(name = SK_"setRemoved", exec = setRemoved, overhead = setOverhead) &
            , bench_type(name = SK_"getRemoved", exec = getRemoved, overhead = setOverhead) &
            ]
    arraySize = [( 2_IK**isize, isize = 1_IK, NSIZE )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "setRemoved() vs. getRemoved()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1, size(bench, 1, IK))

        loopOverArraySize: do isize = 1, NSIZE

            write(*,"(*(g0,:,' '))") "Benchmarking with size", arraySize(isize)

            do i = 1, size(bench, 1, IK)
                bench(i)%timing = bench(i)%getTiming(minsec = 0.05_RK)
            end do

            write(fileUnit,"(*(g0,:,','))") arraySize(isize), (bench(i)%timing%mean, i = 1, size(bench, 1, IK))

        end do loopOverArraySize
        write(*,"(*(g0,:,' '))") dummy
        write(*,"(*(g0,:,' '))")

    close(fileUnit)

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! procedure wrappers.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setOverhead()
        call initialize()
        call finalize()
    end subroutine

    subroutine initialize()
        array = repeat(pattern, arraySize(isize))
    end subroutine

    subroutine finalize()
        dummy = dummy .and. 0_IK < len(array, IK)
        deallocate(array)
    end subroutine

    subroutine setRemoved()
        block
            use pm_arrayRemove, only: setRemoved
            call initialize()
            call setRemoved(array, pattern)
            call finalize()
        end block
    end subroutine

    subroutine getRemoved()
        block
            use pm_arrayRemove, only: getRemoved
            call initialize()
            array = getRemoved(array, pattern)
            call finalize()
        end block
    end subroutine

end program benchmark