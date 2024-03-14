! Test the performance of getBin() with and without the optional external function `isLess` input argument.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_kind, only: IK, LK, SK, RK, RKC => RK32
    use pm_distUnif, only: setUnifRand
    use pm_arraySpace, only: getLinSpace
    use pm_mathMinMax, only: getMinMax
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: i                            !<  The procedure benchmark counter.
    integer(IK)                         :: isize                        !<  The Array size counter.
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK)     , parameter         :: NSIZE = 20_IK                !<  The number of benchmark ranks.
    integer(IK)     , parameter         :: NBENCH = 2_IK                !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(NSIZE)             !<  The sizes of the benchmark Array.
    logical(LK)                         :: dummy = .true._LK            !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    real(RKC)         , allocatable     :: Array(:)                     !<  The Array to find instances of value within.
    real(RKC)                           :: value                        !<  The value to find.
    real(RKC)                           :: Limit(2)                     !<  The array limits.
    type(bench_type)                    :: bench(NBENCH)                !<  The Benchmark array.
    integer(IK)                         :: index, bin
    real(RK)                            :: MeanTime(2) = 0._RK

    bench(1) = bench_type(name = SK_"default"  , exec = default    , overhead = setOverhead)
    bench(2) = bench_type(name = SK_"isLessArg", exec = isLessArg  , overhead = setOverhead)

    arraySize = [( 2_IK**(isize+1), isize = 1_IK, NSIZE )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "isLessArg vs. default"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,', '))") "arraySize", (bench(i)%name, i = 1, NBENCH)

        loopOverArraySize: do isize = 1, NSIZE

            write(*,"(*(g0,:,' '))") "Benchmarking with size", arraySize(isize)

            call random_number(Limit)
            Limit = getMinMax(Limit(1), Limit(2))
            Array = getLinSpace(Limit(1), Limit(2), count = arraySize(isize))

            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming(minsec = 0.03_RK)
                MeanTime(i) = MeanTime(i) + bench(i)%timing%mean
            end do

            do i = NBENCH, 1, -1
                bench(i)%timing = bench(i)%getTiming(minsec = 0.03_RK)
                MeanTime(i) = MeanTime(i) + bench(i)%timing%mean
            end do

            write(fileUnit,"(*(g0,:,', '))") arraySize(isize), MeanTime / 2

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
        call setUnifRand(index, 1_IK, arraySize(isize))
        value = Array(index)
    end subroutine

    subroutine finalize()
        dummy = dummy .and. bin == index
    end subroutine

    subroutine isLessArg()
        use pm_arraySearch, only: getBin
        call initialize()
        bin = getBin(Array, value, isLess = isLess)
        call finalize()
    end subroutine

    subroutine default()
        use pm_arraySearch, only: getBin
        call initialize()
        bin = getBin(Array, value)
        call finalize()
    end subroutine

    pure function isLess(value, segment) result(less)
        use pm_kind, only: IK, LK
        integer(IK) , intent(in)    :: value, segment
        logical(LK)                 :: less
        less = value < segment
    end function

end program benchmark