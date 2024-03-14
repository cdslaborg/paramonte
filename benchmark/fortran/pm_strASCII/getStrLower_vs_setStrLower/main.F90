! Test the performance of `getStrLower()` vs. `setStrLower()`.
program benchmark

    use pm_bench, only: bench_type
    use pm_kind, only: IK, LK, RK, SK
    use pm_distUnif, only: setUnifRand
    use iso_fortran_env, only: error_unit

    implicit none

    integer(IK)                         :: i                            !<  The procedure benchmark counter.
    integer(IK)                         :: isize                        !<  The Array size counter.
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK)     , parameter         :: NSIZE = 18_IK                !<  The number of benchmark ranks.
    integer(IK)     , parameter         :: NBENCH = 2_IK                !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(NSIZE)             !<  The sizes of the benchmark Array.
    logical(LK)                         :: dummy = .true._LK            !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    character(:, SK), allocatable       :: Array, arrayNew              !<  The Array whose elements will have to be found.
    type(bench_type)                    :: bench(NBENCH)                !<  The Benchmark array.
    real(RK)                            :: TimingMean(2) = 0._RK

    bench(1) = bench_type(name = SK_"setStrLower", exec = setStrLower , overhead = setOverhead)
    bench(2) = bench_type(name = SK_"getStrLower", exec = getStrLower , overhead = setOverhead)

    arraySize = [( 2_IK**isize, isize = 1_IK, NSIZE )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "setStrLower() vs. getStrLower()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1, NBENCH)

        loopOverArraySize: do isize = 1, NSIZE

            write(*,"(*(g0,:,' '))") "Benchmarking with size", arraySize(isize)
            allocate(character(arraySize(isize)) :: Array, arrayNew)
            call setUnifRand(Array)

            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming(minsec = 0.025_RK)
                TimingMean(i) = TimingMean(i) + bench(i)%timing%mean
            end do

            do i = NBENCH, 1, -1
                bench(i)%timing = bench(i)%getTiming(minsec = 0.025_RK)
                TimingMean(i) = TimingMean(i) + bench(i)%timing%mean
            end do

            deallocate(Array, arrayNew)
            write(fileUnit,"(*(g0,:,','))") arraySize(isize), TimingMean/2

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
        !allocate( character(arraySize(isize)) :: Array )
        arrayNew = Array
    end subroutine

    subroutine finalize()
        dummy = dummy .and. arrayNew(1:1) == "a"
        !deallocate(Array)
    end subroutine

    subroutine setStrLower()
        block
            use pm_strASCII, only: setStrLower!, setStrUpper
            call initialize()
            call setStrLower(arrayNew)
           !call setStrUpper(arrayNew)
           !call setStrLower(arrayNew)
            call finalize()
        end block
    end subroutine

    subroutine getStrLower()
        block
            use pm_strASCII, only: getStrLower!, getStrUpper
            call initialize()
            arrayNew = getStrLower(arrayNew)
           !arrayNew = getStrUpper(arrayNew)
           !arrayNew = getStrLower(arrayNew)
            call finalize()
        end block
    end subroutine

end program benchmark