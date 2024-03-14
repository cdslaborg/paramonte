! Test the overhead of calling `setUnifRand()` vs. Fortran intrinsic procedure `random_number()`.
program benchmark

    use pm_kind, only: IK, RK, SK
    use pm_bench, only: bench_type
    use pm_distUnif, only: setUnifRand, xoshiro256ssw_type

    implicit none

    type(xoshiro256ssw_type)            :: rngx
    integer(IK)                         :: i                !<  The procedure benchmark counter.
    integer(IK)                         :: iarr             !<  The array size counter.
    integer(IK)                         :: fileUnit         !<  The output file unit for benchmark results.
    integer(IK)     , parameter         :: NARR = 11_IK     !<  The number of benchmark array sizes.
    integer(IK)                         :: arraySize(NARR)  !<  The benchmarking array sizes.
    real(RK)        , allocatable       :: rand(:)          !<  The Random vector.
    real(RK)                            :: dummy = 0._RK    !<  The dummy computation to prevent aggressive optimizations.
    type(bench_type), allocatable       :: bench(:)         !<  The Benchmark array.

    rngx = xoshiro256ssw_type()
    bench = [ bench_type(name = SK_"random_number ", exec = random_number, overhead = setOverhead) &
            , bench_type(name = SK_"setUnifRandRNGD", exec = setUnifRandRNGD, overhead = setOverhead) &
            , bench_type(name = SK_"setUnifRandRNGX", exec = setUnifRandRNGX, overhead = setOverhead) &
            ]

    arraySize = [( 2_IK**iarr, iarr = 1_IK, NARR )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "setUnifRand() vs. random_number()."
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1, size(bench))

        loopOverMatrixSize: do iarr = 1, NARR

            allocate(rand(arraySize(iarr)))
            write(*,"(*(g0,:,' '))") "Benchmarking setUnifRand() vs. random_number() with array size", arraySize(iarr)

            do i = 1, size(bench)
                bench(i)%timing = bench(i)%getTiming()
            end do

            write(fileUnit,"(*(g0,:,','))") arraySize(iarr), (bench(i)%timing%mean, i = 1, size(bench))
            deallocate(rand)

        end do loopOverMatrixSize
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
        dummy = dummy + rand(1)
    end subroutine

    subroutine setUnifRandRNGD()
        block
            call setUnifRand(rand)
            call getDummy()
        end block
    end subroutine

    subroutine setUnifRandRNGX()
        block
            call setUnifRand(rngx, rand)
            call getDummy()
        end block
    end subroutine

    subroutine random_number()
        block
            intrinsic :: random_number
            call random_number(rand)
            call getDummy()
        end block
    end subroutine

end program benchmark