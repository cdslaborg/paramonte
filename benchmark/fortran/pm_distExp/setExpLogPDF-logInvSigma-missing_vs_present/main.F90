program benchmark

    use iso_fortran_env, only: error_unit
    use pm_bench, only: bench_type
    use pm_kind, only: IK, RK, SK

    implicit none

    integer(IK)                         :: i                            !<  The procedure benchmark counter.
    integer(IK)                         :: isize                        !<  The logPDF size counter.
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK)     , parameter         :: NSIZE = 18_IK                !<  The number of benchmark ranks.
    integer(IK)     , parameter         :: NBENCH = 2_IK                !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(0:NSIZE)           !<  The sizes of the benchmark logPDF.
    real(RK)        , allocatable       :: logPDF(:), point(:)           !<  The benchmark array.
    real(RK)        , allocatable       :: logInvSigma(:), invSigma(:)      !<  The distribution parameters.
    real(RK)                            :: dummy = 0._RK                !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    type(bench_type)                    :: bench(NBENCH)                !<  The Benchmark array.

    bench(1) = bench_type(name = SK_"logInvSigmaMissing", exec = logInvSigmaMissing , overhead = setOverhead)
    bench(2) = bench_type(name = SK_"logInvSigmaPresent", exec = logInvSigmaPresent , overhead = setOverhead)

    arraySize = [( 2_IK**isize, isize = 0_IK, NSIZE )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' vs. '))") (bench(i)%name, i = 1, NBENCH)
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1, NBENCH)

        loopOverArraySize: do isize = 0, NSIZE

            write(*,"(*(g0,:,' '))") "Benchmarking with size", arraySize(isize)

            allocate(logPDF(arraySize(isize)), point(arraySize(isize)), invSigma(arraySize(isize)))
            call random_number(point)
            call random_number(invSigma)
            logInvSigma = log(invSigma)
            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming(minsec = 0.05_RK)
            end do
            deallocate(logPDF, point, invSigma)

            write(fileUnit,"(*(g0,:,','))") arraySize(isize), (bench(i)%timing%mean, i = 1, NBENCH)

        end do loopOverArraySize

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
        !call random_number(invSigma)
        !logInvSigma = log(invSigma)
    end subroutine

    subroutine finalize()
        dummy = dummy + logPDF(1)
    end subroutine

    subroutine logInvSigmaPresent()
        use pm_distExp, only: setExpLogPDF
        call initialize()
        if (arraySize(isize) > 1_IK) then
            call setExpLogPDF(logPDF, point, invSigma, logInvSigma)
        else
            call setExpLogPDF(logPDF(1), point(1), invSigma(1), logInvSigma(1))
        end if
        call finalize()
    end subroutine

    subroutine logInvSigmaMissing()
        use pm_distExp, only: getExpLogPDF
        call initialize()
        if (arraySize(isize) > 1_IK) then
            logPDF = getExpLogPDF(point, invSigma = invSigma)
        else
            logPDF(1) = getExpLogPDF(point(1), invSigma = invSigma(1))
        end if
        call finalize()
    end subroutine

end program benchmark
