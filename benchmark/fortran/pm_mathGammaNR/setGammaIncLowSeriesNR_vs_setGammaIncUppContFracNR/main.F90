! Test the performance of `setGammaIncLowSeriesNR()` vs. `setGammaIncUppContFracNR()`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_kind, only: IK, LK, RKG => RK, RK, SK
    use pm_arraySpace, only: setLogSpace
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: i                                !<  The procedure benchmark counter.
    integer(IK)                         :: ipnt                             !<  The Point size counter.
    integer(IK)                         :: fileUnit                         !<  The output file unit for benchmark results.
    integer(IK)     , parameter         :: NPNT = 20_IK                     !<  The number of benchmark points.
    integer(IK)     , parameter         :: NBENCH = 2_IK                    !<  The number of benchmark procedures.
    real(RKG)       , parameter         :: kappa = 1._RKG                   !<  The kappa parameter of the Gamma function.
    real(RKG)       , parameter         :: logGammaKappa = log_gamma(kappa) !<  The log_gamma(kappa).
   !real(RKG)       , parameter         :: TOL = 1000 * epsilon(0._RKG)     !<  The tolerance.
    real(RKG)                           :: Point(NPNT)                      !<  The benchmark array.
    real(RKG)                           :: dummy = 0._RKG                   !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    real(RKG)                           :: gammaInc                         !<  The incomplete gamma function.
    integer(IK)                         :: info                             !<  The logical convergence failure check.
    type(bench_type)                    :: bench(NBENCH)                    !<  The Benchmark array.

    bench(1) = bench_type(name = SK_"setGammaIncLowSeriesNR", exec = setGammaIncLowSeriesNR, overhead = setOverhead)
    bench(2) = bench_type(name = SK_"setGammaIncUppContFracNR", exec = setGammaIncUppContFracNR, overhead = setOverhead)

    call setLogSpace(Point, logx1 = log(0.04_RKG), logx2 = log(100._RKG))


    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' vs. '))") (bench(i)%name, i = 1, NBENCH)
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "x / (kappa + 1)", (bench(i)%name, i = 1, NBENCH)

        loopOverPoint: do ipnt = 1, NPNT

            write(*,"(*(g0,:,' '))") "Benchmarking with point", Point(ipnt)

            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming(minsec = 0.05_RK)
            end do

            write(fileUnit,"(*(g0,:,','))") Point(ipnt) / (kappa + 1), (bench(i)%timing%mean, i = 1, NBENCH)

        end do loopOverPoint

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
        if (info < 0_IK) error stop
        dummy = dummy + gammaInc
    end subroutine

    subroutine setGammaIncLowSeriesNR()
        block
            use pm_mathGammaNR, only: setGammaIncLowSeriesNR
            call setGammaIncLowSeriesNR ( gammaIncLow = gammaInc &
                                        , x = Point(ipnt) &
                                        , logGammaKappa = logGammaKappa &
                                        , kappa = kappa &
                                        , info = info &
                                       !, tol = TOL &
                                        )
            call finalize()
        end block
    end subroutine

    subroutine setGammaIncUppContFracNR()
        block
            use pm_mathGammaNR, only: setGammaIncUppContFracNR
            call setGammaIncUppContFracNR( gammaIncUpp = gammaInc &
                                        , x = Point(ipnt) &
                                        , logGammaKappa = logGammaKappa &
                                        , kappa = kappa &
                                        , info = info &
                                       !, tol = TOL &
                                        )
            call finalize()
        end block
    end subroutine

end program benchmark