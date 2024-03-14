! Test the overhead of calling `setNormRand()` vs. Fortran intrinsic procedure `random_number()`.
program benchmark

    use pm_kind, only: IK, RK, RKC => RK, SK
    use pm_distNorm, only: xoshiro256ssw_type
    use pm_distNorm, only: setNormRand
    use pm_distNorm, only: getZigNorm
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: ibench
    integer(IK)                         :: nlay
    integer(IK)                         :: iset             !<  The ziggurat set counter.
    integer(IK)                         :: fileUnit         !<  The output file unit for benchmark results.
    integer(IK)         , parameter     :: NSET = 13_IK     !<  The number of ziggurat sets.
    integer(IK)         , parameter     :: NSIM = 10000_IK  !<  The number of ziggurat sets.
    real(RKC)           , allocatable   :: zig(:,:)         !<  The ziggurat set.
    real(RKC)                           :: rand(NSIM)       !<  The Random vector.
    real(RKC)                           :: dummy = 0._RKC   !<  The dummy computation to prevent aggressive optimizations.
    type(bench_type)    , allocatable   :: bench(:)
    type(xoshiro256ssw_type)            :: rng

    rng = xoshiro256ssw_type()
    bench = [ bench_type(name = SK_"setNormRandX256", exec = setNormRandX256, overhead = setOverhead) &
            , bench_type(name = SK_"setNormRandFRNG", exec = setNormRandFRNG, overhead = setOverhead) &
            ]

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "ZigguratLayerCount", (bench(ibench)%name, ibench = 1, size(bench))

        loopOverZigSets: do iset = 1, NSET

            nlay = 2**iset
            zig = getZigNorm(nlay, dummy)
            write(*,"(*(g0,:,' '))") "Benchmarking setNormRand() with ziggurat set size and abserr", nlay, dummy
            do ibench = 1, size(bench)
                bench(ibench)%timing = bench(ibench)%getTiming(miniter = 10_IK)
            end do
            write(fileUnit, "(*(g0,:,','))") nlay, (bench(ibench)%timing%mean / NSIM, ibench = 1, size(bench))
            write(*,"(*(g0,:,' '))") dummy

        end do loopOverZigSets

    close(fileUnit)

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! procedure wrappers.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setOverhead()
        call getDummy()
    end subroutine

    subroutine getDummy()
        dummy = dummy + sum(rand)
    end subroutine

    subroutine setNormRandFRNG()
        call setNormRand(rand, zig)
        call getDummy()
    end subroutine

    subroutine setNormRandX256()
        call setNormRand(rng, rand, zig)
        call getDummy()
    end subroutine

end program benchmark