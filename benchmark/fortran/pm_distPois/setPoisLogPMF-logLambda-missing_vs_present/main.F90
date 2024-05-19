program benchmark

    use iso_fortran_env, only: error_unit
    use pm_bench, only: bench_type
    use pm_distUnif, only: xoshiro256ssw_type
    use pm_distUnif, only: setUnifRand
    use pm_kind, only: SK, IK, RK, RKG => RK

    implicit none

    integer(IK)                         :: i                            !<  The procedure benchmark counter.
    integer(IK)                         :: isize                        !<  The logPMF size counter.
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK)     , parameter         :: NSIZE = 18_IK                !<  The number of benchmark ranks.
    integer(IK)     , parameter         :: NBENCH = 2_IK                !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(0:NSIZE)           !<  The sizes of the benchmark logPMF.
    integer(IK)     , allocatable       :: count(:)                     !<  The benchmark point.
    real(RKG)       , allocatable       :: logPMF(:)                    !<  The benchmark array.
    real(RKG)       , allocatable       :: logLambda(:), lambda(:)      !<  The distribution parameters.
    real(RKG)                           :: dummy = 0._RKG               !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    type(bench_type)                    :: bench(NBENCH)                !<  The Benchmark array.
    type(xoshiro256ssw_type)            :: rng

    rng = xoshiro256ssw_type()

    bench(1) = bench_type(name = SK_"logLambdaMissing", exec = logLambdaMissing , overhead = setOverhead)
    bench(2) = bench_type(name = SK_"logLambdaPresent", exec = logLambdaPresent , overhead = setOverhead)

    arraySize = [( 2_IK**isize, isize = 0_IK, NSIZE )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' vs. '))") (bench(i)%name, i = 1, NBENCH)
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1, NBENCH)

        loopOverarraySize: do isize = 0, NSIZE

            write(*,"(*(g0,:,' '))") "Benchmarking with size", arraySize(isize)

            allocate(logPMF(arraySize(isize)), count(arraySize(isize)), lambda(arraySize(isize)))
            call setUnifRand(rng, count, 0_IK, 1023_IK)
            call random_number(lambda)
            lambda = 1._RKG - lambda
            logLambda = log(lambda)
            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming(minsec = 0.05_RK)
            end do
            deallocate(logPMF, count, lambda)

            write(fileUnit,"(*(g0,:,','))") arraySize(isize), (bench(i)%timing%mean, i = 1, NBENCH)

        end do loopOverarraySize

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
        !call random_number(lambda)
        !logLambda = log(lambda)
    end subroutine

    subroutine finalize()
        dummy = dummy + logPMF(1)
    end subroutine

    subroutine logLambdaPresent()
        use pm_distPois, only: setPoisLogPMF
        call initialize()
        if (arraySize(isize) > 1_IK) then
            call setPoisLogPMF(logPMF, count, lambda, logLambda)
        else
            call setPoisLogPMF(logPMF(1), count(1), lambda(1), logLambda(1))
        end if
        call finalize()
    end subroutine

    subroutine logLambdaMissing()
        use pm_distPois, only: getPoisLogPMF
        call initialize()
        if (arraySize(isize) > 1_IK) then
            logPMF = getPoisLogPMF(count, lambda = lambda)
        else
            logPMF(1) = getPoisLogPMF(count(1), lambda = lambda(1))
        end if
        call finalize()
    end subroutine

end program benchmark