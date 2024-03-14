program benchmark

    use iso_fortran_env, only: error_unit
    use pm_bench, only: bench_type
    use pm_kind, only: IK, RK, SK

    implicit none

    integer(IK)                         :: i
    integer(IK)                         :: isize
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK)     , parameter         :: NSIZE = 15_IK                !<  The number of benchmark procedures.
    integer(IK)     , parameter         :: NBENCH = 3_IK                !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(NSIZE)             !<  The number of benchmark iterations.
    real(RK)        , allocatable       :: InvScale(:)                  !<  The input argument.
    real(RK)        , allocatable       :: shape(:)                     !<  The input argument.
    real(RK)                            :: dummy = 0._RK                !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    type(bench_type)                    :: Bench(NBENCH)                !<  The Benchmark array.

    Bench(1) = bench_type(name = SK_"scalarAddition", exec = scalarAddition, overhead = setOverhead)
    Bench(2) = bench_type(name = SK_"logNormFacWithShape", exec = logNormFacWithShape, overhead = setOverhead)
    Bench(3) = bench_type(name = SK_"logNormFacWithShapeScale", exec = logNormFacWithShapeScale, overhead = setOverhead)

    arraySize = [( 2**isize, isize = 1, NSIZE )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' vs. '))") (Bench(i)%name, i = 1, NBENCH)
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (Bench(i)%name, i = 1, NBENCH)

        loopOverArraySize: do isize = 1, NSIZE

            write(*,"(*(g0,:,' '))") "Benchmarking with rank", arraySize(isize)

            allocate(shape(arraySize(isize)), InvScale(arraySize(isize)))
            do i = 1, NBENCH
                Bench(i)%timing = Bench(i)%getTiming(minsec = 0.05_RK)
            end do
            deallocate(shape, InvScale)

            write(fileUnit,"(*(g0,:,','))") arraySize(isize), (Bench(i)%timing%mean, i = 1, NBENCH)

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
    end subroutine

    subroutine initialize()
        call random_number(shape)
        shape = shape + epsilon(0._RK)
        InvScale = shape + 1._RK
    end subroutine

    subroutine scalarAddition()
        call initialize()
        dummy = dummy + sum(scalarAdditionCallee(shape,InvScale))
    end subroutine

    impure elemental function scalarAdditionCallee(shape, invScale) result(summation)
        real(RK), intent(in) :: shape, invScale
        real(RK) :: summation
        summation = shape + invScale
    end function

    subroutine logNormFacWithShape()
        use pm_distGenExpGamma, only: getGenExpGammaLogPDFNF
        call initialize()
        dummy = dummy + sum(getGenExpGammaLogPDFNF(shape))
    end subroutine

    subroutine logNormFacWithShapeScale()
        use pm_distGenExpGamma, only: getGenExpGammaLogPDFNF
        call initialize()
        dummy = dummy + sum(getGenExpGammaLogPDFNF(shape, InvScale))
    end subroutine

end program benchmark
