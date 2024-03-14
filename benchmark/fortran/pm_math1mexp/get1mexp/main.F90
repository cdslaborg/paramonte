! Test the performance of `get1mexp()` with and without the selection `control` argument.
program benchmark

    use pm_bench, only: bench_type
    use pm_arraySpace, only: getLinSpace
    use pm_kind, only: IK, LK, RKC => RK64, RK, SK
    use iso_fortran_env, only: error_unit

    implicit none

    integer(IK)                         :: i                    !<  The procedure benchmark counter.
    integer(IK)                         :: iarr                 !<  The array size counter.
    integer(IK)                         :: fileUnit             !<  The output file unit for benchmark results.
    integer(IK) , parameter             :: NBENCH = 2_IK        !<  The number of benchmark procedures.
    real(RKC)                           :: dummySum = 0._RKC    !<  The dummy computation to prevent the compiler from aggressive optimizations.
    real(RKC)   , allocatable           :: X(:)                 !<  The benchmark value.
    real(RKC)                           :: onemexp              !<  The cumulative proportion of the exponential of the `X(:)`.
    type(bench_type)                    :: bench(NBENCH)        !<  The Benchmark array.
    logical(LK)                         :: underflowEnabled     !<  The logical flag indicating whether an array with many instances of underflow should be generated.

    bench(1) = bench_type(name = SK_"get1mexpSelection", exec = get1mexpSelection, overhead = setOverhead)
    bench(2) = bench_type(name = SK_"get1mexpSequence", exec = get1mexpSequence, overhead = setOverhead)

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "get1mexp(...) vs. get1mexp(..., control = selection)"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "X", (bench(i)%name, i = 1, NBENCH)
        X = getLinSpace(x1 = -2 * log(huge(onemexp)), x2 = log(epsilon(onemexp)), count = 20_IK)
        loopOverArraySize: do iarr = 1, size(X)

            write(*,"(*(g0,:,' '))") "Benchmarking with X", X(iarr)
            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming(minsec = 0.07_RK)
            end do
            write(fileUnit,"(*(g0,:,','))") X(iarr), (bench(i)%timing%mean, i = 1, NBENCH)

        end do loopOverArraySize
        write(*,"(*(g0,:,' '))") dummySum
        write(*,"(*(g0,:,' '))")

    close(fileUnit)

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! procedure wrappers.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setOverhead()
        call getDummy()
    end subroutine

    subroutine getDummy()
        dummySum = dummySum + onemexp
    end subroutine

    subroutine get1mexpSelection()
        use pm_math1mexp, only: get1mexp, selection
        onemexp = get1mexp(X(iarr), selection)
        call getDummy()
    end subroutine

    subroutine get1mexpSequence()
        use pm_math1mexp, only: get1mexp
        onemexp = get1mexp(X(iarr))
        call getDummy()
    end subroutine

end program benchmark