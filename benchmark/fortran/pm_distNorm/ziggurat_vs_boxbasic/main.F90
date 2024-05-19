module zig_mod
    use pm_distNorm, only: xoshiro256ssw_type
    use pm_distNorm, only: getZigNorm
    use pm_kind, only: RKG => RK
    implicit none
    real(RKG) :: abserr
    real(RKG), allocatable :: zig(:,:)
end module zig_mod

! Test the performance of `setNormRandZiggurat()` vs. `setNormRandBoxBasic()`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_bench, only: bench_type, benchMulti_type
    use pm_io, only: display_type
    use pm_kind, only: SK, IK, LK, RK
    use zig_mod

    implicit none

    integer(IK)                         :: isim
    integer(IK)                         :: itime
    integer(IK)                         :: ibench
    integer(IK)         , parameter     :: NSIM = 10000_IK ! must be even number.
    real(RKG)                           :: rand(NSIM) = 0._RKG
    real(RKG)                           :: dummy = 0._RKG
    type(benchMulti_type)               :: bench
    type(display_type)                  :: disp
    integer(IK)                         :: miniter = 10_IK
    type(xoshiro256ssw_type)            :: rng

    zig = getZigNorm(256_IK, abserr)
    rng = xoshiro256ssw_type()
    bench = benchMulti_type([ bench_type(name = SK_"setNormRandZiggurat", exec = setNormRandZiggurat, overhead = setOverhead, minsec = 0._RK, miniter = miniter) &
                            , bench_type(name = SK_"setNormRandZigX256S", exec = setNormRandZigX256S, overhead = setOverhead, minsec = 0._RK, miniter = miniter) &
                            , bench_type(name = SK_"setNormRandBoxBasic", exec = setNormRandBoxBasic, overhead = setOverhead, minsec = 0._RK, miniter = miniter) &
                            ], sorted = .true._LK, repeat = 1_IK)

    disp = display_type()
    call disp%show(bench%name, tmsize = 1_IK, bmsize = 1_IK)
    disp = display_type(file = "main.out")

    write(disp%unit, "(*(g0,:,','))") (bench%case(ibench)%name, ibench = 1, size(bench%case))
    do itime = 1, size(bench%case(1)%timing%values)
        call disp%show([(bench%case(ibench)%timing%values(itime) / NSIM, ibench = 1, size(bench%case))])
    end do

    write(*,"(*(g0,:,' '))") dummy
    write(*,"(*(g0,:,' '))")

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! procedure wrappers.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setOverhead()
        call finalize()
    end subroutine

    subroutine finalize()
        dummy = dummy + sum(rand)
    end subroutine

    subroutine setNormRandZiggurat()
        use pm_distNorm, only: setNormRand
        call setNormRand(rand)!, zig)
        call finalize()
    end subroutine

    subroutine setNormRandZigX256S()
        use pm_distNorm, only: setNormRand
        call setNormRand(rng, rand)!, zig)
        call finalize()
    end subroutine

    subroutine setNormRandBoxBasic()
        use pm_distNorm, only: setNormRandBox
        call random_number(rand)
        call setNormRandBox(rand(1:NSIM-1:2), rand(2:NSIM:2))
        call finalize()
    end subroutine

end program benchmark