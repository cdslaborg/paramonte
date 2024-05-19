module internal
    use pm_kind, only: RKG => RK
    implicit none
contains
    impure subroutine setNormRandBox(rand1, rand2)
        real(RKG)   , intent(out)   :: rand1, rand2
        real(RKG)                   :: factor, rsq
        do
            call random_number(rand1)
            call random_number(rand2)
            rand1 = 2._RKG * rand1 - 1._RKG
            rand2 = 2._RKG * rand2 - 1._RKG
            rsq = rand1**2 + rand2**2
            if (0._RKG < rsq .and. rsq < 1._RKG) exit
        end do
        factor = sqrt(-2._RKG * log(rsq) / rsq)
        rand1 = rand1 * factor
        rand2 = rand2 * factor
    end subroutine
end module


! Test the performance of `setNormRandBoxBasic()` vs. `setNormRandBoxPolar()`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_bench, only: bench_type, benchMulti_type
    use pm_io, only: display_type
    use pm_kind, only: SK, IK, LK, RK
    use internal, only: RKG

    implicit none

    integer(IK)                         :: itime
    integer(IK)                         :: ibench
    real(RKG)                           :: rand1 = 0._RKG, rand2 = 0._RKG
    real(RKG)                           :: dummy = 0._RKG
    type(benchMulti_type)               :: bench
    type(display_type)                  :: disp
    integer(IK)                         :: miniter = 10000_IK

    bench = benchMulti_type([ bench_type(name = SK_"setNormRandBoxBasic", exec = setNormRandBoxBasic, overhead = setOverhead, minsec = 0._RK, miniter = miniter) &
                            , bench_type(name = SK_"setNormRandBoxPolar", exec = setNormRandBoxPolar, overhead = setOverhead, minsec = 0._RK, miniter = miniter) &
                            , bench_type(name = SK_"setNormRandBoxPolarImpure", exec = setNormRandBoxPolarImpure, overhead = setOverhead, minsec = 0._RK, miniter = miniter) &
                           !, bench_type(name = SK_"getNormRandBox", exec = getNormRandBox, overhead = setOverhead, minsec = 0._RK, miniter = miniter) &
                            ], sorted = .true._LK, repeat = 1_IK)

    disp = display_type()
    call disp%show(bench%name, tmsize = 1_IK, bmsize = 1_IK)
    disp = display_type(file = "main.out")

    write(disp%unit, "(*(g0,:,','))") (bench%case(ibench)%name, ibench = 1, size(bench%case))
    do itime = 1, size(bench%case(1)%timing%values)
        call disp%show([(bench%case(ibench)%timing%values(itime), ibench = 1, size(bench%case))])
    end do

    write(*,"(*(g0,:,' '))") dummy
    write(*,"(*(g0,:,' '))")

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! procedure wrappers.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setOverhead()
        call initialize()
        call finalize()
    end subroutine

    subroutine initialize()
        call random_number(rand1)
        call random_number(rand2)
    end subroutine

    subroutine finalize()
        dummy = dummy + rand1 + rand2
    end subroutine

    subroutine setNormRandBoxBasic()
        use pm_distNorm, only: setNormRandBox
        call initialize()
        call setNormRandBox(rand1, rand2)
        call finalize()
    end subroutine

    subroutine setNormRandBoxPolar()
        use pm_distNorm, only: setNormRandBox
        logical(LK) :: failed
        call initialize()
        do
            call setNormRandBox(rand1, rand2, failed)
            if (.not. failed) exit
            call random_number(rand1)
            call random_number(rand2)
        end do
        call finalize()
    end subroutine

    !subroutine getNormRandBox()
    !    block
    !        use pm_distNorm, only: getNormRandBox
    !        rand1 = getNormRandBox(mean = 0._RKG)
    !        rand2 = getNormRandBox(mean = 0._RKG)
    !        call finalize()
    !    end block
    !end subroutine

    subroutine setNormRandBoxPolarImpure()
        use internal, only: setNormRandBox
        call setNormRandBox(rand1, rand2)
        call finalize()
    end subroutine

end program benchmark