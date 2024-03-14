! Test the performance of `getPoisLogPMF()` vs. `setPoisLogPMF()`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_kind, only: SK, IK, RK, RKC => RK
    use pm_distUnif, only: xoshiro256ssw_type
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: i                            !<  The procedure benchmark counter.
    integer(IK)                         :: isize                        !<  The array size counter.
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK)     , parameter         :: NSIZE = 18_IK                !<  The number of benchmark ranks.
    integer(IK)     , parameter         :: NBENCH = 2_IK                !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(NSIZE)             !<  The sizes of the benchmark array.
    integer(IK)     , allocatable       :: count(:)                     !<  The benchmark point.
    real(RKC)       , allocatable       :: array(:)                     !<  The benchmark array.
    real(RKC)                           :: dummy = 0._RKC               !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    type(bench_type)                    :: bench(NBENCH)                !<  The Benchmark array.
    type(xoshiro256ssw_type)            :: rng

    rng = xoshiro256ssw_type()

    bench(1) = bench_type(name = SK_"getPoisLogPMF", exec = getPoisLogPMF , overhead = setOverhead)
    bench(2) = bench_type(name = SK_"setPoisLogPMF", exec = setPoisLogPMF , overhead = setOverhead)

    arraySize = [( 2_IK**isize, isize = 1_IK, NSIZE )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' vs. '))") (bench(i)%name, i = 1, NBENCH)
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1, NBENCH)

        loopOverarraySize: do isize = 1, NSIZE

            write(*,"(*(g0,:,' '))") "Benchmarking with size", arraySize(isize)

            allocate(array(arraySize(isize)), count(arraySize(isize)))
            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming(minsec = 0.05_RK)
            end do
            deallocate(array, count)

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
        use pm_distUnif, only: setUnifRand
        call setUnifRand(rng, count, 0_IK, 1023_IK)
    end subroutine

    subroutine finalize()
        dummy = dummy + array(1)
    end subroutine

    subroutine setPoisLogPMF()
        block
            use pm_distPois, only: setPoisLogPMF
            call initialize()
            call setPoisLogPMF(array, count, lambda = 1._RKC)
            call finalize()
        end block
    end subroutine

    subroutine getPoisLogPMF()
        block
            use pm_distPois, only: getPoisLogPMF
            call initialize()
            array = getPoisLogPMF(count, lambda = 1._RKC)
            call finalize()
        end block
    end subroutine

end program benchmark