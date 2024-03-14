! Test the performance of `getInserted()` vs. `setInserted()`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_arrayRange, only: getRange
    use pm_kind, only: IK, LK, RK, SK
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: i                            !<  The procedure benchmark counter.
    integer(IK)                         :: isize                        !<  The array size counter.
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK) , parameter             :: NSIZE = 15_IK                !<  The number of benchmark ranks.
    integer(IK) , parameter             :: NBENCH = 2_IK                !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(NSIZE)             !<  The sizes of the benchmark array.
    logical(LK)                         :: dummy = .true._LK            !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    real(RK)                            :: insertion = 1._RK            !<  The insertion array.
    real(RK)    , allocatable           :: array(:), arrayNew(:)        !<  The array within which insertion is to be made.
    integer(IK) , allocatable           :: index(:)                     !<  The indices of positions of insertions.
    type(bench_type)                    :: bench(NBENCH)                !<  The Benchmark array.

    bench(1) = bench_type(name = SK_"setInserted", exec = setInserted , overhead = setOverhead)
    bench(2) = bench_type(name = SK_"getInserted", exec = getInserted , overhead = setOverhead)

    arraySize = [( 2_IK**isize, isize = 1_IK, NSIZE )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "setInserted() vs. getInserted()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1, NBENCH)

        loopOverArraySize: do isize = 1, NSIZE

            write(*,"(*(g0,:,' '))") "Benchmarking with size", arraySize(isize)

            index = getRange(1_IK, arraySize(isize))
            allocate(arrayNew(arraySize(isize)+size(index)))
            allocate(array(arraySize(isize)), source = 1._RK)
            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming(minsec = 0.05_RK)
            end do
            deallocate(array, arrayNew)

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
        call random_number(insertion)
    end subroutine

    subroutine finalize()
        dummy = dummy .and. arrayNew(1) == 0.5_RK
    end subroutine

    subroutine setInserted()
        block
            use pm_arrayInsert, only: setInserted
            call initialize()
            call setInserted(arrayNew, array, insertion, index, positive = .true._LK, sorted = .true._LK)
            call finalize()
        end block
    end subroutine

    subroutine getInserted()
        block
            use pm_arrayInsert, only: getInserted
            call initialize()
            arrayNew(:) = getInserted(array, insertion, index, positive = .true._LK, sorted = .true._LK)
            call finalize()
        end block
    end subroutine

end program benchmark