program benchmark

    use pm_kind, only: IK, LK, RK, SK
    use iso_fortran_env, only: error_unit
    use pm_arrayShuffle, only: setShuffled
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: i                    !<  The procedure benchmark counter.
    integer(IK)                         :: isize                !<  The array size counter.
    integer(IK)                         :: fileUnit             !<  The output file unit for benchmark results.
    integer(IK) , parameter             :: NSIZE = 15_IK        !<  The number of benchmark sizes.
    integer(IK) , parameter             :: NBENCH = 3_IK        !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(NSIZE)     !<  The sizes of the benchmark array.
    real(RK)                            :: dummy = 0._RK        !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    real(RK)    , allocatable           :: array(:)             !<  The array whose elements will have to be remapped.
    integer(IK) , allocatable           :: index(:)             !<  The array of indices.
    type(bench_type)                    :: bench(NBENCH)        !<  The Benchmark array.

    bench(1) = bench_type(name = SK_"setRemapped", exec = setRemapped, overhead = setOverhead)
    bench(2) = bench_type(name = SK_"getRemapped", exec = getRemapped, overhead = setOverhead)
    bench(3) = bench_type(name = SK_"direct", exec = direct , overhead = setOverhead)

    arraySize = [( 2_IK**isize, isize = 1_IK, NSIZE )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "setRemapped() vs. getRemapped() vs. direct()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1, NBENCH)

        loopOverArraySize: do isize = 1, NSIZE

            write(*,"(*(g0,:,' '))") "Benchmarking with size", arraySize(isize)

            index = [( i, i = 1, arraySize(isize) )]
            allocate(array(arraySize(isize)))
            call random_number(array)
            call setShuffled(index)
            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming(minsec = 0.07_RK)
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
        call finalize()
    end subroutine

    subroutine finalize()
        dummy = dummy + array(1)
    end subroutine

    subroutine setRemapped()
        block
            use pm_arrayRemap, only: setRemapped, reverse
            call setRemapped(array, index, reverse)
            call finalize()
        end block
    end subroutine

    subroutine getRemapped()
        block
            use pm_arrayRemap, only: getRemapped, reverse
            array = getRemapped(array, index, reverse)
            call finalize()
        end block
    end subroutine

    subroutine direct()
        array = array(index(arraySize(isize):1:-1))
        call finalize()
    end subroutine

end program benchmark
