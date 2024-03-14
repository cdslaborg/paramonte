! Test the performance of setRemoved() with a vector `pattern` vs. scalar `pattern`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_bench, only: bench_type
    use pm_kind, only: IK, LK, RK, SK

    implicit none

    integer(IK)                         :: i                            !<  The procedure benchmark counter.
    integer(IK)                         :: isize                        !<  The Array size counter.
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK) , parameter             :: NSIZE = 15_IK                !<  The number of benchmark ranks.
    integer(IK) , parameter             :: NBENCH = 2_IK                !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(NSIZE)             !<  The sizes of the benchmark Array.
    logical(LK)                         :: dummy = .true._LK            !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    integer(IK) , allocatable           :: Array(:)                     !<  The Array to remove pattern from.
    integer(IK) , parameter             :: pattern(1) = 0_IK            !<  The pattern to remove.
    type(bench_type)                    :: bench(NBENCH)                !<  The Benchmark array.

    bench(1) = bench_type(name = SK_"scalarPattern", exec = scalarPattern , overhead = setOverhead)
    bench(2) = bench_type(name = SK_"vectorPattern", exec = vectorPattern , overhead = setOverhead)

    arraySize = [( 2_IK**isize, isize = 1_IK, NSIZE )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "scalarPattern() vs. vectorPattern()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1_IK, NBENCH)

        loopOverArraySize: do isize = 1, NSIZE

            write(*,"(*(g0,:,' '))") "Benchmarking with size", arraySize(isize)
            allocate(Array(arraySize(isize)))

            do i = 1_IK, NBENCH
                bench(i)%timing = bench(i)%getTiming(minsec = 0.05_RK)
            end do

            deallocate(Array)
            write(fileUnit,"(*(g0,:,','))") arraySize(isize), (bench(i)%timing%mean, i = 1_IK, NBENCH)

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
        Array(:) = 1_IK
    end subroutine

    subroutine finalize()
        dummy = dummy .and. size(Array, kind = IK) < 1_IK
    end subroutine

    subroutine scalarPattern()
        use pm_arrayRemove, only: setRemoved
        call initialize()
        call setRemoved(Array, pattern(1))
        call finalize()
    end subroutine

    subroutine vectorPattern()
        block
            use pm_arrayRemove, only: setRemoved
            call initialize()
            call setRemoved(Array, pattern)
            call finalize()
        end block
    end subroutine

end program benchmark