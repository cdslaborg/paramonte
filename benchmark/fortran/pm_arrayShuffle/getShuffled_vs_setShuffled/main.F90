program benchmark

    use pm_kind, only: IK, LK, RK, SK
    use iso_fortran_env, only: error_unit
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: i                    !<  The procedure benchmark counter.
    integer(IK)                         :: isize                !<  The array size counter.
    integer(IK)                         :: fileUnit             !<  The output file unit for benchmark results.
    integer(IK) , parameter             :: NSIZE = 15_IK        !<  The number of benchmark sizes.
    integer(IK) , parameter             :: NBENCH = 2_IK        !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(NSIZE)     !<  The sizes of the benchmark array.
    real(RK)                            :: dummy = 0._RK        !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    real(RK)    , allocatable           :: Array(:)             !<  The array whose elements will have to be remapped.
    type(bench_type)                    :: bench(NBENCH)        !<  The Benchmark array.

    bench(1) = bench_type(name = SK_"setShuffled", exec = setShuffled, overhead = setOverhead)
    bench(2) = bench_type(name = SK_"getShuffled", exec = getShuffled, overhead = setOverhead)

    arraySize = [( 2_IK**isize, isize = 1_IK, NSIZE )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "setShuffled() vs. getShuffled()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1, NBENCH)

        loopOverArraySize: do isize = 1, NSIZE

            write(*,"(*(g0,:,' '))") "Benchmarking with size", arraySize(isize)

            allocate(Array(arraySize(isize)))
            call random_number(Array)
            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming(minsec = 0.07_RK)
            end do
            deallocate(Array)

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
        dummy = dummy + Array(1)
    end subroutine

    subroutine setShuffled()
        block
            use pm_arrayShuffle, only: setShuffled
            call setShuffled(Array)
            call finalize()
        end block
    end subroutine

    subroutine getShuffled()
        block
            use pm_arrayShuffle, only: getShuffled
            Array = getShuffled(Array)
            call finalize()
        end block
    end subroutine

end program benchmark
