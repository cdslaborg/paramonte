! Test the performance of `getSelected()` vs. `setSelected()`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_kind, only: IK, LK, RK, SK
    use pm_arrayRange, only: getRange
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: i                            !<  The procedure benchmark counter.
    integer(IK)                         :: isize                        !<  The array size counter.
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK) , parameter             :: NSIZE = 20_IK                !<  The number of benchmark ranks.
    integer(IK) , parameter             :: NBENCH = 2_IK                !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(NSIZE)             !<  The sizes of the benchmark array.
    logical(LK)                         :: dummy = .true._LK            !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    real(RK)                            :: selection = 1._RK            !<  The selection array.
    real(RK)    , allocatable           :: array(:)                     !<  The array whose `rank`th smallest is to be selected.
    integer(IK)                         :: rank                         !<  The rank of the selected element of the array.
    type(bench_type)                    :: bench(NBENCH)                !<  The Benchmark array.

    bench(1) = bench_type(name = SK_"setSelected", exec = setSelected , overhead = setOverhead)
    bench(2) = bench_type(name = SK_"getSelected", exec = getSelected , overhead = setOverhead)

    arraySize = [( 2_IK**isize, isize = 1_IK, NSIZE )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "setSelected() vs. getSelected()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1, NBENCH)

        loopOverArraySize: do isize = 1, NSIZE

            write(*,"(*(g0,:,' '))") "Benchmarking with size", arraySize(isize)

            allocate(array(arraySize(isize)))
            rank = 1_IK ! arraySize(isize)
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
        call random_number(array)
    end subroutine

    subroutine finalize()
        dummy = dummy .and. selection == 0.5_RK
    end subroutine

    subroutine setSelected()
        block
            use pm_arraySelect, only: setSelected
            call initialize()
            call setSelected(selection, array, rank)
            call finalize()
        end block
    end subroutine

    subroutine getSelected()
        block
            use pm_arraySelect, only: getSelected
            call initialize()
            selection = getSelected(array, rank)
            call finalize()
        end block
    end subroutine

end program benchmark