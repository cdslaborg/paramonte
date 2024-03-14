! Test the performance of `getReplaced()` vs. `setReplaced()`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_kind, only: IK, LK, RK, SK
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: i                            !<  The procedure benchmark counter.
    integer(IK)                         :: isize                        !<  The array size counter.
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK) , parameter             :: NSIZE = 12_IK                !<  The number of benchmark ranks.
    integer(IK) , parameter             :: NBENCH = 3_IK                !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(NSIZE)             !<  The sizes of the benchmark array.
    logical(LK)                         :: dummy = .true._LK            !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    character(:, SK), allocatable       :: array                        !<  The array whose elements will have to be replaced.
    character(*, SK), parameter         :: pattern = "a"                !<  The pattern to be replaced in array.
    character(*, SK), parameter         :: Replacement = "**"           !<  The replacement for pattern in array.
    type(bench_type)                    :: bench(NBENCH)                !<  The Benchmark array.

    bench(1) = bench_type(name = SK_"setReplaced", exec = setReplaced , overhead = setOverhead)
    bench(2) = bench_type(name = SK_"getReplaced", exec = getReplaced , overhead = setOverhead)
    bench(3) = bench_type(name = SK_"getReplaced_recurs_alloc", exec = getReplaced_recurs_alloc , overhead = setOverhead)

    arraySize = [( 2_IK**isize, isize = 1_IK, NSIZE )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "setReplaced() vs. getReplaced() vs. getReplaced_recurs_alloc()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1, NBENCH)

        loopOverArraySize: do isize = 1, NSIZE

            write(*,"(*(g0,:,' '))") "Benchmarking with size", arraySize(isize)

            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming(minsec = 0.05_RK)
            end do

            write(fileUnit,"(*(g0,:,','))") arraySize(isize), (bench(i)%timing%mean, i = 1, NBENCH)

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
        allocate( character(arraySize(isize)) :: array )
        !block
        !    use pm_distUnif, only: setUnifRand
        !    call setUnifRand(array)
        !end block
        block
            array(:) = repeat(pattern, len(array))
        end block
    end subroutine

    subroutine finalize()
        dummy = dummy .and. array(1:1) == pattern
        deallocate(array)
    end subroutine

    subroutine setReplaced()
        block
            use pm_arrayReplace, only: setReplaced
            call initialize()
            call setReplaced(array, pattern, Replacement)
            call finalize()
        end block
    end subroutine

    subroutine getReplaced()
        block
            use pm_arrayReplace, only: getReplaced
            call initialize()
            array = getReplaced(array, pattern, Replacement)
            call finalize()
        end block
    end subroutine

    subroutine getReplaced_recurs_alloc()
        block
            use pm_arrayReplace, only: getReplaced_recurs_alloc
            call initialize()
            array = getReplaced_recurs_alloc(array, pattern, Replacement)
            call finalize()
        end block
    end subroutine

end program benchmark