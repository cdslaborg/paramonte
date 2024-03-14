! Test the performance of `getFactorial()` vs. `getLogFactorial()`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_kind, only: IK, LK, SK, RK64, IK64
    use pm_arrayRange, only: getRange
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: i                                !<  The procedure benchmark counter.
    integer(IK)                         :: ipnt                             !<  The Point size counter.
    integer(IK)                         :: fileUnit                         !<  The output file unit for benchmark results.
    integer(IK)     , parameter         :: NPNT = 10_IK                     !<  The number of benchmark points.
    integer(IK)     , parameter         :: NBENCH = 2_IK                    !<  The number of benchmark procedures.
    real(RK64)                          :: Point_RK64(NPNT)                 !<  The benchmark array.
    real(RK64)                          :: dummy = 0._RK64                  !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    real(RK64)                          :: factorial_RK64 = 0._RK64         !<  The factorial.
    integer(IK64)                       :: factorial_IK64 = 0_IK64          !<  The factorial.
    integer(IK64)                       :: Point_IK64(NPNT)                 !<  The benchmark array.
    type(bench_type)                    :: bench(NBENCH)                    !<  The Benchmark array.

    bench(1) = bench_type(name = SK_"getFactorial", exec = getFactorial , overhead = setOverhead)
    bench(2) = bench_type(name = SK_"getLogFactorial", exec = getLogFactorial , overhead = setOverhead)

    Point_IK64 = getRange(start = 2_IK64, stop = 20_IK64, step = 2_IK64)
    Point_RK64 = real(Point_IK64, RK64)
    

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' vs. '))") (bench(i)%name, i = 1, NBENCH)
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "Point", (bench(i)%name, i = 1, NBENCH)

        loopOverPoint: do ipnt = 1, NPNT

            write(*,"(*(g0,:,' '))") "Benchmarking with point", Point_IK64(ipnt)

            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming()
            end do

            write(fileUnit,"(*(g0,:,','))") Point_IK64(ipnt), (bench(i)%timing%mean, i = 1, NBENCH)

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
        dummy = dummy + factorial_IK64 + factorial_RK64
    end subroutine

    subroutine getFactorial()
        block
            use pm_mathFactorial, only: getFactorial
            factorial_IK64 = getFactorial(Point_IK64(ipnt))
            call finalize()
        end block
    end subroutine

    subroutine getLogFactorial()
        block
            use pm_mathFactorial, only: getLogFactorial
            factorial_RK64 = log(getLogFactorial(Point_RK64(ipnt)))
            call finalize()
        end block
    end subroutine

end program benchmark