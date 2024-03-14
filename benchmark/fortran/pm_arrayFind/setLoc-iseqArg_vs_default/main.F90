! Test the performance of setLoc() with and without the optional external function `iseq` input argument.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_kind, only: IK, LK, RK, SK
    use pm_arrayRange, only: getRange
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: i                            !<  The procedure benchmark counter.
    integer(IK)                         :: nloc                         !<  The length of `loc`.
    integer(IK)                         :: isize                        !<  The array size counter.
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK)     , parameter         :: blindness = 1                !<  The search blindness.
    integer(IK)     , parameter         :: NSIZE = 20_IK                !<  The number of benchmark ranks.
    integer(IK)     , parameter         :: NBENCH = 2_IK                !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(NSIZE)             !<  The sizes of the benchmark array.
    logical(LK)                         :: dummy = .true._LK            !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    integer(IK)     , allocatable       :: instance(:)                  !<  The array of indices of instances of pattern to find in array.
    integer(IK)     , allocatable       :: loc(:)                       !<  The array of indices of instances of pattern in array.
    integer(IK)     , allocatable       :: array(:)                     !<  The array to find instances of pattern within.
    integer(IK)                         :: pattern                      !<  The pattern to find.
    type(bench_type)                    :: bench(NBENCH)                !<  The Benchmark array.
    real(RK)                            :: unifrnd, MeanTime(2) = 0._RK

    bench(1) = bench_type(name = SK_"default", exec = default , overhead = setOverhead)
    bench(2) = bench_type(name = SK_"iseqArg", exec = iseqArg , overhead = setOverhead)

    arraySize = [( 2_IK**isize, isize = 1_IK, NSIZE )]
    loc = [integer(IK) ::] ! essential for `setLoc()`.
    nloc = 0_IK

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "iseqArg vs. default"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1_IK, NBENCH)

        loopOverArraySize: do isize = 1_IK, NSIZE

            write(*,"(*(g0,:,' '))") "Benchmarking with size", arraySize(isize)
            instance = getRange(1_IK, arraySize(isize))
            allocate(array(arraySize(isize)))

            do i = 1_IK, NBENCH
                bench(i)%timing = bench(i)%getTiming(minsec = 0.03_RK)
                MeanTime(i) = MeanTime(i) + bench(i)%timing%mean
            end do

            do i = NBENCH, 1_IK, -1_IK
                bench(i)%timing = bench(i)%getTiming(minsec = 0.03_RK)
                MeanTime(i) = MeanTime(i) + bench(i)%timing%mean
            end do

            deallocate(array, instance)
            write(fileUnit,"(*(g0,:,','))") arraySize(isize), MeanTime / 2_IK

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
        call random_number(unifrnd)
        pattern = nint(unifrnd, kind = IK)
        array(:) = pattern + 1_IK
    end subroutine

    subroutine finalize()
        dummy = dummy .and. size(loc, 1, IK) == nloc
    end subroutine

    subroutine iseqArg()
        use pm_arrayFind, only: setLoc
        call initialize()
        call setLoc(loc, nloc, array, pattern, iseq, 1_IK)
        call finalize()
    end subroutine

    subroutine default()
        block
            use pm_arrayFind, only: setLoc
            call initialize()
            call setLoc(loc, nloc, array, pattern, 1_IK)
            call finalize()
        end block
    end subroutine

    pure function iseq(arraysegment, pattern) result(equivalent)
        use pm_kind, only: IK, LK
        integer(IK) , intent(in)    :: pattern, arraySegment
        logical(LK)                 :: equivalent
        equivalent = arraySegment == pattern
    end function

end program benchmark