program benchmark

    use pm_kind, only: SK, IK, LK, RKG => RK, IKG => IK, LKG => LK
    use pm_distUnif, only: xoshiro256ssg_type
    use pm_distUnif, only: xoshiro256ssw_type
    use pm_distUnif, only: setUnifRand
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: i, j, fileUnit
    integer(IK)         , parameter     :: NSIM = 100000_IK
    logical(LKG)                        :: dumm_LK = .false._LKG    !<  The dummy value to prevent aggressive optimization.
    logical(LKG)                        :: rand_LK(NSIM)            !<  The Random vector.
    integer(IKG)                        :: rand_IK(NSIM)            !<  The Random vector.
    real(RKG)                           :: rand_RK(NSIM)            !<  The Random vector.
    type(bench_type)    , allocatable   :: bench(:)                 !<  The Benchmark array.
    type(xoshiro256ssg_type)            :: xoshiro256ssg
    type(xoshiro256ssw_type)            :: xoshiro256ssw
    xoshiro256ssg = xoshiro256ssg_type()
    xoshiro256ssw = xoshiro256ssw_type()

    bench = [ bench_type(name = SK_"random_number_LK", exec = random_number_LK, overhead = setOverhead_LK) &
            , bench_type(name = SK_"xoshiro256ssg_type_LK", exec = xoshiro256ssg_type_LK, overhead = setOverhead_LK) &
            , bench_type(name = SK_"xoshiro256ssw_type_LK", exec = xoshiro256ssw_type_LK, overhead = setOverhead_LK) &
            , bench_type(name = SK_"random_number_IK", exec = random_number_IK, overhead = setOverhead_LK) &
            , bench_type(name = SK_"xoshiro256ssg_type_IK", exec = xoshiro256ssg_type_IK, overhead = setOverhead_IK) &
            , bench_type(name = SK_"xoshiro256ssw_type_IK", exec = xoshiro256ssw_type_IK, overhead = setOverhead_IK) &
            , bench_type(name = SK_"random_number_RK", exec = random_number_RK, overhead = setOverhead_RK) &
            , bench_type(name = SK_"xoshiro256ssg_type_RK", exec = xoshiro256ssg_type_RK, overhead = setOverhead_RK) &
            , bench_type(name = SK_"xoshiro256ssw_type_RK", exec = xoshiro256ssw_type_RK, overhead = setOverhead_RK) &
            ]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "xoshiro256ssg_type() vs. xoshiro256ssw_type()."
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") (bench(i)%name, i = 1, size(bench))
        do i = 1, size(bench)
            bench(i)%timing = bench(i)%getTiming()
        end do
        do j = 1, minval([(size(bench(i)%timing%values), i = 1, size(bench))])
            write(fileUnit,"(*(g0,:,','))") (bench(i)%timing%values(j) / NSIM, i = 1, size(bench))
        end do
        write(*,"(*(g0,:,' '))")

    close(fileUnit)

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! procedure wrappers.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setOverhead_LK()
        call getDummy_LK()
    end subroutine

    subroutine setOverhead_IK()
        call getDummy_IK()
    end subroutine

    subroutine setOverhead_RK()
        call getDummy_RK()
    end subroutine

    subroutine getDummy_LK()
        dumm_LK = dumm_LK .or. count(rand_LK) == 0
    end subroutine

    subroutine getDummy_IK()
        dumm_LK = dumm_LK .or. any(rand_IK == 0_IKG)
    end subroutine

    subroutine getDummy_RK()
        dumm_LK = dumm_LK .or. any(rand_RK == 0._RKG)
    end subroutine

    subroutine random_number_LK()
        call setUnifRand(rand_LK)
        !block
        !    real :: rand
        !    call random_number(rand)
        !    rand_LK = logical(rand < 0.5, LKG)
        !    call getDummy_LK()
        !end block
    end subroutine

    subroutine xoshiro256ssg_type_LK()
        call setUnifRand(xoshiro256ssg, rand_LK)
        call getDummy_LK()
    end subroutine

    subroutine xoshiro256ssw_type_LK()
        call setUnifRand(xoshiro256ssw, rand_LK)
        call getDummy_LK()
    end subroutine


    subroutine random_number_IK()
        call setUnifRand(rand_IK)
        call getDummy_IK()
    end subroutine

    subroutine xoshiro256ssg_type_IK()
        call setUnifRand(xoshiro256ssg, rand_IK)
        call getDummy_IK()
    end subroutine

    subroutine xoshiro256ssw_type_IK()
        call setUnifRand(xoshiro256ssw, rand_IK)
        call getDummy_IK()
    end subroutine


    subroutine random_number_RK()
        call setUnifRand(rand_RK)
        call getDummy_RK()
    end subroutine

    subroutine xoshiro256ssg_type_RK()
        call setUnifRand(xoshiro256ssg, rand_RK)
        call getDummy_RK()
    end subroutine

    subroutine xoshiro256ssw_type_RK()
        call setUnifRand(xoshiro256ssw, rand_RK)
        call getDummy_RK()
    end subroutine

end program benchmark