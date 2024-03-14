! Test the performance of `getGenExpGammaLogPDF()` vs. `setGenExpGammaLogPDF()`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_bench, only: bench_type
    use pm_kind, only: IK, RK, SK

    implicit none

    integer(IK)                         :: i                            !<  The procedure benchmark counter.
    integer(IK)                         :: isize                        !<  The Array size counter.
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK)     , parameter         :: NSIZE = 18_IK                !<  The number of benchmark ranks.
    integer(IK)     , parameter         :: NBENCH = 2_IK                !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(NSIZE)             !<  The sizes of the benchmark Array.
    real(RK)        , allocatable       :: Array(:), Point(:)           !<  The benchmark array.
    real(RK)                            :: dummy = 0._RK                !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    type(bench_type)                    :: bench(NBENCH)                !<  The Benchmark array.

    bench(1) = bench_type(name = SK_"setGenExpGammaLogPDF", exec = setGenExpGammaLogPDF , overhead = setOverhead)
    bench(2) = bench_type(name = SK_"getGenExpGammaLogPDF", exec = getGenExpGammaLogPDF , overhead = setOverhead)

    arraySize = [( 2_IK**isize, isize = 1_IK, NSIZE )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' vs. '))") (bench(i)%name, i = 1, NBENCH)
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1, NBENCH)

        loopOverArraySize: do isize = 1, NSIZE

            write(*,"(*(g0,:,' '))") "Benchmarking with size", arraySize(isize)

            allocate(Array(arraySize(isize)), Point(arraySize(isize)))
            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming(minsec = 0.05_RK)
            end do
            deallocate(Array, Point)

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
        call random_number(Point)
    end subroutine

    subroutine finalize()
        dummy = dummy + Array(1)
    end subroutine

    subroutine setGenExpGammaLogPDF()
        block
            use pm_distGenExpGamma, only: setGenExpGammaLogPDF
            call initialize()
            call setGenExpGammaLogPDF(Array, Point)
            call finalize()
        end block
    end subroutine

    subroutine getGenExpGammaLogPDF()
        block
            use pm_distGenExpGamma, only: getGenExpGammaLogPDF
            call initialize()
            Array = getGenExpGammaLogPDF(Point)
            call finalize()
        end block
    end subroutine

end program benchmark