! Test the performance of `getFactorial()` vs. `getLogFactorial()`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_kind, only: IK, LK, SK, RKD, IKD
    use pm_arrayRange, only: getRange
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: i                                !<  The procedure benchmark counter.
    integer(IK)                         :: ipnt                             !<  The Point size counter.
    integer(IK)                         :: fileUnit                         !<  The output file unit for benchmark results.
    integer(IK)     , parameter         :: NPNT = 10_IK                     !<  The number of benchmark points.
    integer(IK)     , parameter         :: NBENCH = 2_IK                    !<  The number of benchmark procedures.
    real(RKD)                           :: Point_RKD(NPNT)                  !<  The benchmark array.
    real(RKD)                           :: dummy = 0._RKD                   !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    real(RKD)                           :: factorial_RKD = 0._RKD           !<  The factorial.
    integer(IKD)                        :: factorial_IKD = 0_IKD            !<  The factorial.
    integer(IKD)                        :: point_IKD(NPNT)                  !<  The benchmark array.
    type(bench_type)                    :: bench(NBENCH)                    !<  The Benchmark array.

    bench(1) = bench_type(name = SK_"getFactorial", exec = getFactorial , overhead = setOverhead)
    bench(2) = bench_type(name = SK_"getLogFactorial", exec = getLogFactorial , overhead = setOverhead)

    point_IKD = getRange(start = 2_IKD, stop = 20_IKD, step = 2_IKD)
    Point_RKD = real(point_IKD, RKD)
    

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' vs. '))") (bench(i)%name, i = 1, NBENCH)
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "Point", (bench(i)%name, i = 1, NBENCH)

        loopOverPoint: do ipnt = 1, NPNT

            write(*,"(*(g0,:,' '))") "Benchmarking with point", point_IKD(ipnt)

            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming()
            end do

            write(fileUnit,"(*(g0,:,','))") point_IKD(ipnt), (bench(i)%timing%mean, i = 1, NBENCH)

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
        dummy = dummy + factorial_IKD + factorial_RKD
    end subroutine

    subroutine getFactorial()
        block
            use pm_mathFactorial, only: getFactorial
            factorial_IKD = getFactorial(point_IKD(ipnt))
            call finalize()
        end block
    end subroutine

    subroutine getLogFactorial()
        block
            use pm_mathFactorial, only: getLogFactorial
            factorial_RKD = log(getLogFactorial(Point_RKD(ipnt)))
            call finalize()
        end block
    end subroutine

end program benchmark