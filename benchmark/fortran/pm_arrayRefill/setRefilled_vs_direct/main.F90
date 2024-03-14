program benchmark

    use iso_fortran_env, only: error_unit
    use pm_bench, only: bench_type
    use pm_kind, only: IK, RK, SK

    implicit none

    integer(IK)                         :: i                    !<  The procedure benchmark counter.
    integer(IK)                         :: isize                !<  The array size counter.
    integer(IK)                         :: fileUnit             !<  The output file unit for benchmark results.
    integer(IK) , parameter             :: NSIZE = 21_IK        !<  The number of benchmark sizes.
    integer(IK) , parameter             :: NBENCH = 2_IK        !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(NSIZE)     !<  The sizes of the benchmark array.
    real(RK)                            :: dummy = 0._RK        !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    real(RK)    , allocatable           :: array(:)             !<  The array whose elements will have to be replaced.
    type(bench_type)                    :: bench(NBENCH)        !<  The Benchmark array.

    bench(1) = bench_type(name = SK_"setRefilled", exec = setRefilled , overhead = setOverhead)
    bench(2) = bench_type(name = SK_"direct", exec = direct , overhead = setOverhead)

    arraySize = [( 2_IK**isize, isize = 1_IK, NSIZE )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "setRefilled() vs. direct()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1, NBENCH)

        loopOverArraySize: do isize = 1, NSIZE - 1

            write(*,"(*(g0,:,' '))") "Benchmarking with size", arraySize(isize)

            allocate(array(arraySize(isize)))
            call random_number(array)
            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming(minsec = 0.1_RK, miniter = 20_IK)
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
        use pm_arrayShuffle, only: setShuffled
        call setShuffled(array(1:arraySize(isize)))
    end subroutine

    subroutine finalize()
        dummy = dummy + array(1)
    end subroutine

    subroutine setRefilled()
        block
            use pm_arrayRefill, only: setRefilled
            call initialize()
                call setRefilled(array, 0._RK, arraySize(isize+1))
            call finalize()
        end block
    end subroutine

    subroutine direct()
        real(RK), allocatable :: Temp(:)
        call initialize()
            Temp = array(1:arraySize(isize))
            deallocate(array)
            allocate(array(arraySize(isize+1)))
            array(1:arraySize(isize)) = Temp
            array(arraySize(isize)+1:arraySize(isize+1)) = 0._RK
        call finalize()
        deallocate(Temp)
    end subroutine

end program benchmark