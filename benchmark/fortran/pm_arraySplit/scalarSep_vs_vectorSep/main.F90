! Test the performance of split with `vector_sep` vs. `scalar_sep`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_container, only: cvi_type
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
    integer(IK) , allocatable           :: array(:)                     !<  The array to split.
    integer(IK) , parameter             :: sep(1) = 0_IK                !<  The sepiter.
    type(bench_type)                    :: bench(NBENCH)                !<  The Benchmark array.
    type(cvi_type), allocatable   :: ArraySplit(:)                !<  The array parts after split.

    bench(1) = bench_type(name = SK_"scalar_sep", exec = scalar_sep , overhead = setOverhead)
    bench(2) = bench_type(name = SK_"vector_sep", exec = vector_sep , overhead = setOverhead)

    arraySize = [( 2_IK**isize, isize = 1_IK, NSIZE )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "scalar_sep() vs. vector_sep()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1, NBENCH)

        loopOverArraySize: do isize = 1, NSIZE

            write(*,"(*(g0,:,' '))") "Benchmarking with size", arraySize(isize)
            allocate(array(arraySize(isize)))

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
        array(:) = 1_IK
    end subroutine

    subroutine finalize()
        dummy = dummy .and. size(ArraySplit, kind = IK) == 1_IK
    end subroutine

    subroutine scalar_sep()
        use pm_arraySplit, only: setSplit
        call initialize()
        call setSplit(ArraySplit, array, sep(1))
        call finalize()
    end subroutine

    subroutine vector_sep()
        block
            use pm_arraySplit, only: setSplit
            call initialize()
            call setSplit(ArraySplit, array, sep)
            call finalize()
        end block
    end subroutine

end program benchmark