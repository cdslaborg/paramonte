! Test the performance of `getDot()` with and without the selection `control` argument.
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
    real(RKG)                           :: dotres               !<  The summation.
    type(bench_type)    , allocatable   :: bench(:)             !<  The Benchmark array.
    logical(LK)                         :: underflowEnabled     !<  The logical flag indicating whether an array with many instances of underflow should be generated.

    bench = [ bench_type(name = SK_"dot_product()", exec = getDotFortran, overhead = setOverhead) &
            , bench_type(name = SK_"fablocked", exec = getDotFAB, overhead = setOverhead) &
            , bench_type(name = SK_"nablocked", exec = getDotNAB, overhead = setOverhead) &
            , bench_type(name = SK_"kahanbabu", exec = getDotKAB, overhead = setOverhead) &
            , bench_type(name = SK_"iteration", exec = getDotIte, overhead = setOverhead) &
            , bench_type(name = SK_"recursion", exec = getDotRec, overhead = setOverhead) &
            ]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "dot_product() vs. getDot()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        call setUnifRand(array)
        !truth = array; call setCumSum(truth)
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
        dumsum = dumsum + dotres
    end subroutine

    subroutine getDotFortran()
        dotres = dot_product(array(1 : arrlen), array(1 : arrlen))
        call getDummy()
    end subroutine

    subroutine getDotFAB()
        use pm_mathSum, only: getDot!, fablocked
        dotres = getDot(array(1 : arrlen), array(1 : arrlen))!, fablocked)
        call getDummy()
    end subroutine

    subroutine getDotNAB()
        use pm_mathSum, only: getDot, nablocked
        dotres = getDot(array(1 : arrlen), array(1 : arrlen), nablocked)
        call getDummy()
    end subroutine

    subroutine getDotKAB()
        use pm_mathSum, only: getDot, kahanbabu
        dotres = getDot(array(1 : arrlen), array(1 : arrlen), kahanbabu)
        call getDummy()
    end subroutine

    subroutine getDotIte()
        use pm_mathSum, only: getDot, iteration
        dotres = getDot(array(1 : arrlen), array(1 : arrlen), iteration)
        call getDummy()
    end subroutine

    subroutine getDotRec()
        use pm_mathSum, only: getDot, recursion
        dotres = getDot(array(1 : arrlen), array(1 : arrlen), recursion)
        call getDummy()
    end subroutine

end program benchmark