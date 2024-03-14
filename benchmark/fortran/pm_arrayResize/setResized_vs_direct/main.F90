program benchmark

    use iso_fortran_env, only: error_unit
    use pm_kind, only: IK, RK, SK
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: i                    !<  The procedure benchmark counter.
    integer(IK)                         :: isize                !<  The array size counter.
    integer(IK)                         :: fileUnit             !<  The output file unit for benchmark results.
    integer(IK)         , parameter     :: NSIZE = 21_IK        !<  The number of benchmark sizes.
    integer(IK)                         :: arraySize(NSIZE)     !<  The sizes of the benchmark array.
    real(RK)                            :: dummy = 0._RK        !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    real(RK)            , allocatable   :: array(:)             !<  The array whose elements will have to be replaced.
    type(bench_type)    , allocatable   :: bench(:)             !<  The Benchmark array.

    bench = [ bench_type(name = SK_"setResized", exec = setResized , overhead = setOverhead) &
            , bench_type(name = SK_"direct", exec = direct , overhead = setOverhead) &
            ]

    arraySize = [( 2_IK**isize, isize = 1_IK, NSIZE )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "setResized() vs. direct()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1, size(bench, 1, IK))

        loopOverArraySize: do isize = 1, NSIZE - 1

            write(*,"(*(g0,:,' '))") "Benchmarking with size", arraySize(isize)

            allocate(array(arraySize(isize)))
            call random_number(array)
            do i = 1, size(bench, 1, IK)
                bench(i)%timing = bench(i)%getTiming(minsec = 0.1_RK)!, miniter = 20_IK)
            end do
            deallocate(array)

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
        use pm_arrayShuffle, only: setShuffled
        call setShuffled(array(1:arraySize(isize)))
    end subroutine

    subroutine finalize()
        dummy = dummy + array(1)
    end subroutine

    subroutine setResized()
        block
            use pm_arrayResize, only: setResized
            call initialize()
                call setResized(array, arraySize(isize+1))
            call finalize()
            call initialize()
                call setResized(array, arraySize(isize))
            call finalize()
        end block
    end subroutine

    subroutine direct()
        real(RK), allocatable :: temp(:)
        call initialize()
            temp = array
            deallocate(array)
            allocate(array(arraySize(isize+1)))
            array(1:arraySize(isize)) = temp
        call finalize()
        call initialize()
            temp = array(1:arraySize(isize))
            array = temp
        call finalize()
        deallocate(temp)
    end subroutine

end program benchmark