! Test the performance of `getLoc()` vs. `setLoc()`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_kind, only: IK, LK, RK, SK
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: i                            !<  The procedure benchmark counter.
    integer(IK)                         :: nloc                         !<  The length of `loc`.
    integer(IK)                         :: isize                        !<  The array size counter.
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK)     , parameter         :: blindness = 1                !<  The search blindness.
    integer(IK)     , parameter         :: NSIZE = 18_IK                !<  The number of benchmark ranks.
    integer(IK)     , parameter         :: NBENCH = 2_IK                !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(NSIZE)             !<  The sizes of the benchmark array.
    logical(LK)                         :: dummy = .true._LK            !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    integer(IK)     , allocatable       :: loc(:)                       !<  The array of indices of instances of pattern in array.
    character(:, SK), allocatable       :: array                        !<  The array whose elements will have to be found.
    character(1,SK)                     :: pattern                      !<  The pattern to be found from array.
    type(bench_type)                    :: bench(NBENCH)                !<  The Benchmark array.

    bench(1) = bench_type(name = SK_"setLoc", exec = setLoc , overhead = setOverhead)
    bench(2) = bench_type(name = SK_"getLoc", exec = getLoc , overhead = setOverhead)

    arraySize = [( 2_IK**isize, isize = 1_IK, NSIZE )]
    loc = [integer(IK) ::] ! essential for `setLoc()`.
    nloc = 0_IK

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "setLoc() vs. getLoc()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1, NBENCH)

        loopOverArraySize: do isize = 1, NSIZE

            write(*,"(*(g0,:,' '))") "Benchmarking with size", arraySize(isize)

            allocate(character(arraySize(isize)) :: array)
            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming(minsec = 0.05_RK)
            end do
            deallocate(array)

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
        use pm_distUnif, only: setUnifRand
        !allocate( character(arraySize(isize)) :: array )
        call setUnifRand(pattern)
        array(:) = repeat(pattern, len(array))
    end subroutine

    subroutine finalize()
        dummy = dummy .and. size(loc, 1, IK) == nloc
        !deallocate(array)
    end subroutine

    subroutine setLoc()
        block
            use pm_arrayFind, only: setLoc
            call initialize()
            call setLoc(loc, nloc, array, pattern, blindness)
            call finalize()
        end block
    end subroutine

    subroutine getLoc()
        block
            use pm_arrayFind, only: getLoc
            call initialize()
            loc = getLoc(array, pattern)
            call finalize()
        end block
    end subroutine

end program benchmark