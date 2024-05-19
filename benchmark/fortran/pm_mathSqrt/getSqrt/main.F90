! Test the performance of `getSqrt()` with and without the selection `control` argument.
program benchmark

    use pm_bench, only: bench_type
    use pm_arrayUnique, only: getUnique
    use pm_arraySpace, only: getLogSpace
    use pm_kind, only: IK, LK, IKG => IK, RK, SK
    use iso_fortran_env, only: error_unit

    implicit none

    integer(IK)                         :: ibench
    integer(IK)                         :: iposint
    integer(IK)                         :: fileUnit
    integer(IKG)                        :: intSqrt, dumm
    integer(IKG)        , allocatable   :: posint(:)
    type(bench_type)    , allocatable   :: bench(:)

    bench = [ bench_type(name = SK_"getSqrtBinary", exec = getSqrtBinary, overhead = setOverhead) &
            , bench_type(name = SK_"getSqrtLinear", exec = getSqrtLinear, overhead = setOverhead) &
            , bench_type(name = SK_"floor_sqrt", exec = floor_sqrt, overhead = setOverhead) &
            ]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' vs. '))") (bench(ibench)%name, ibench = 1, size(bench))
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "Integer", (bench(ibench)%name, ibench = 1, size(bench))
        posint = getUnique(int(getLogSpace(0._RK, log(real(huge(0_IK), RK)), count = 50_IK), IKG))
        loopOverArraySize: do iposint = 1, size(posint)

            write(*,"(*(g0,:,' '))") "Benchmarking with posint = ", posint(iposint)
            do ibench = 1, size(bench)
                bench(ibench)%timing = bench(ibench)%getTiming(minsec = 0.04_RK)
            end do
            write(fileUnit,"(*(g0,:,','))") posint(iposint), (bench(ibench)%timing%mean, ibench = 1, size(bench))

        end do loopOverArraySize
        write(*,"(*(g0,:,' '))") dumm
        write(*,"(*(g0,:,' '))")

    close(fileUnit)

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! procedure wrappers.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setOverhead()
        call getDummy()
    end subroutine

    subroutine getDummy()
        dumm = intSqrt - dumm
    end subroutine

    subroutine getSqrtBinary()
        use pm_mathSqrt, only: getSqrt, binary
        intSqrt = getSqrt(posint(iposint), binary)
        call getDummy()
    end subroutine

    subroutine getSqrtLinear()
        use pm_mathSqrt, only: getSqrt, linear
        intSqrt = getSqrt(posint(iposint), linear)
        call getDummy()
    end subroutine

    subroutine floor_sqrt()
        intSqrt = floor(sqrt(real(posint(iposint), RK)), IKG)
        call getDummy()
    end subroutine

end program benchmark