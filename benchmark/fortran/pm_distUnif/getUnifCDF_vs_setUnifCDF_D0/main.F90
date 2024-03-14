program benchmark

    use pm_kind, only: IK, RK, SK
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: i,j                  !<  The procedure benchmark counter.
    integer(IK)                         :: iarr                 !<  The array size counter.
    integer(IK)                         :: fileUnit             !<  The output file unit for benchmark results.
    integer(IK) , parameter             :: MINITER = 10**5_IK   !<  The number of benchmark procedures.
    integer(IK) , parameter             :: NBENCH = 2_IK        !<  The number of benchmark procedures.
    real(RK)                            :: cdf                  !<  The CDF scalar.
    real(RK)                            :: point                !<  The x scalar.
    real(RK)                            :: dummy = 0._RK        !<  The dummy computation to prevent aggressive optimizations.
    type(bench_type)                    :: bench(NBENCH)        !<  The Benchmark array.

    bench(1) = bench_type(name = SK_"getUnifCDF", exec = getUnifCDF, overhead = setOverhead)
    bench(2) = bench_type(name = SK_"setUnifCDF", exec = setUnifCDF, overhead = setOverhead)

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "getUnifCDF() vs. setUnifCDF()."
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") (bench(i)%name, i = 1, NBENCH)

        call random_number(point)
        do i = 1, NBENCH
            bench(i)%timing = bench(i)%getTiming(miniter = MINITER)
        end do

        do j = 1, MINITER
            write(fileUnit,"(*(g0,:,','))") (max(epsilon(0._RK),bench(i)%timing%values(j)), i = 1, NBENCH)
        end do

        write(*,"(*(g0,:,' '))") dummy
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
        dummy = dummy + cdf
    end subroutine

    subroutine getUnifCDF()
        block
            use pm_distUnif, only: getUnifCDF
            cdf = getUnifCDF(point)
            call finalize()
        end block
    end subroutine

    subroutine setUnifCDF()
        block
            use pm_distUnif, only: setUnifCDF
            call setUnifCDF(cdf, point)
            call finalize()
        end block
    end subroutine

end program benchmark