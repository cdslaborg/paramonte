program benchmark

    use pm_kind, only: IK, RK, SK
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: i                !<  The procedure benchmark counter.
    integer(IK)                         :: iarr             !<  The array size counter.
    integer(IK)                         :: fileUnit         !<  The output file unit for benchmark results.
    integer(IK) , parameter             :: NARR = 16_IK     !<  The number of benchmark array sizes.
    integer(IK) , parameter             :: NBENCH = 2_IK    !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(NARR)  !<  The benchmarking array sizes.
    real(RK)    , parameter             :: LOWER = -2._RK   !<  The lower limit of the CDF.
    real(RK)    , parameter             :: UPPER = +2._RK   !<  The upper limit of the CDF.
    real(RK)    , allocatable           :: cdf(:)           !<  The CDF vector.
    real(RK)    , allocatable           :: Point(:)         !<  The x vector.
    real(RK)                            :: dummy = 0._RK    !<  The dummy computation to prevent aggressive optimizations.
    type(bench_type)                    :: bench(NBENCH)    !<  The Benchmark array.

    bench(1) = bench_type(name = SK_"getUnifCDF", exec = getUnifCDF, overhead = setOverhead)
    bench(2) = bench_type(name = SK_"setUnifCDF", exec = setUnifCDF, overhead = setOverhead)

    arraySize = [( 2_IK**iarr, iarr = 1_IK, NARR )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "getUnifCDF() vs. setUnifCDF()."
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1, NBENCH)

        loopOverMatrixSize: do iarr = 1, NARR

            write(*,"(*(g0,:,' '))") "Benchmarking getUnifCDF() vs. setUnifCDF() with array size", arraySize(iarr)
            allocate(cdf(arraySize(iarr)), Point(arraySize(iarr)))
            call random_number(Point)

            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming()
            end do

            deallocate(cdf, Point)
            write(fileUnit,"(*(g0,:,','))") arraySize(iarr), (bench(i)%timing%mean, i = 1, NBENCH)

        end do loopOverMatrixSize
        write(*,"(*(g0,:,' '))")

    close(fileUnit)

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! procedure wrappers.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setOverhead()
        call finalize()
    end subroutine

    subroutine finalize()
        dummy = dummy + sum(cdf)
    end subroutine

    subroutine getUnifCDF()
        block
            use pm_distUnif, only: getUnifCDF
            cdf(:) = getUnifCDF(Point, LOWER, UPPER)
            call finalize()
        end block
    end subroutine

    subroutine setUnifCDF()
        block
            use pm_distUnif, only: setUnifCDF
            call setUnifCDF(cdf, Point, LOWER, UPPER)
            call finalize()
        end block
    end subroutine

end program benchmark