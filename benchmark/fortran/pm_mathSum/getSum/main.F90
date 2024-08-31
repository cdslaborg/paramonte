! Test the performance of `getSum()` with and without the selection `control` argument.
program benchmark

    use pm_bench, only: bench_type
    use pm_distUnif, only: setUnifRand
    use pm_mathCumSum, only: setCumSum
    use pm_arrayResize, only: setResized
    use pm_kind, only: SK, IK, LK, RKH, RK, RKG => RKD
    use iso_fortran_env, only: error_unit

    implicit none

    integer(IK)                         :: ibench               !<  The procedure benchmark counter.
    integer(IK)                         :: iarr                 !<  The array size counter.
    integer(IK)                         :: arrlen               !<  The array size.
    integer(IK)                         :: fileUnit             !<  The output file unit for benchmark results.
    real(RKG)                           :: dumsum = 0._RKG      !<  The dummy computation to prevent the compiler from aggressive optimizations.
    real(RKG)                           :: array(10**8)         !<  The benchmark values.
    real(RKG)                           :: sumres               !<  The summation.
    type(bench_type)    , allocatable   :: bench(:)             !<  The Benchmark array.

    bench = [ bench_type(name = SK_"sum()", exec = getSumFortran, overhead = setOverhead) &
            , bench_type(name = SK_"fablocked", exec = getSumFAB, overhead = setOverhead) &
            , bench_type(name = SK_"nablocked", exec = getSumNAB, overhead = setOverhead) &
            , bench_type(name = SK_"kahanbabu", exec = getSumKAB, overhead = setOverhead) &
            , bench_type(name = SK_"iteration", exec = getSumIte, overhead = setOverhead) &
            , bench_type(name = SK_"recursion", exec = getSumRec, overhead = setOverhead) &
            ]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "sum() vs. getSum()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        call setUnifRand(array)
        write(fileUnit, "(*(g0,:,','))") "Array Size", (bench(ibench)%name, ibench = 1, size(bench))
        loopOverArraySize: do iarr = 2, 26, 2

            arrlen = 2**iarr
            write(*,"(*(g0,:,' '))") "Benchmarking with array size", arrlen
            do ibench = 1, size(bench)
                bench(ibench)%timing = bench(ibench)%getTiming(minsec = 0.07_RK)
            end do
            write(fileUnit,"(*(g0,:,','))") arrlen, (bench(ibench)%timing%mean, ibench = 1, size(bench))

        end do loopOverArraySize
        write(*,"(*(g0,:,' '))") dumsum
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
        dumsum = dumsum + sumres
    end subroutine

    subroutine getSumFortran()
        sumres = sum(array(1 : arrlen))
        call getDummy()
    end subroutine

    subroutine getSumFAB()
        use pm_mathSum, only: getSum!, fablocked
        sumres = getSum(array(1 : arrlen))!, fablocked)
        call getDummy()
    end subroutine

    subroutine getSumNAB()
        use pm_mathSum, only: getSum, nablocked
        sumres = getSum(array(1 : arrlen), nablocked)
        call getDummy()
    end subroutine

    subroutine getSumKAB()
        use pm_mathSum, only: getSum, kahanbabu
        sumres = getSum(array(1 : arrlen), kahanbabu)
        call getDummy()
    end subroutine

    subroutine getSumIte()
        use pm_mathSum, only: getSum, iteration
        sumres = getSum(array(1 : arrlen), iteration)
        call getDummy()
    end subroutine

    subroutine getSumRec()
        use pm_mathSum, only: getSum, recursion
        sumres = getSum(array(1 : arrlen), recursion)
        call getDummy()
    end subroutine

end program benchmark