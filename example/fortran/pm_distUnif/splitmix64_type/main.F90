program example

    use pm_kind, only: SK, IK, LK
    use pm_err, only: setAsserted
    use pm_io, only: display_type
    use pm_io, only: getFormat
    use pm_distUnif, only: getUnifRand
    use pm_distUnif, only: setUnifRand
    use pm_distUnif, only: splitmix64_type

    implicit none

    integer :: itry, ntry = 5
    type(splitmix64_type) :: rng
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("rng = splitmix64_type()")
                    rng = splitmix64_type()

    block
        use pm_kind, only: IKD
        integer :: rand, lb, ub
        do itry = 1, ntry
        call disp%show("lb = -3_IKD; ub = 5_IKD")
                        lb = -3_IKD; ub = 5_IKD
        call disp%show("rand = getUnifRand(rng, lb, ub)")
                        rand = getUnifRand(rng, lb, ub)
        call disp%show("[lb, rand, ub]")
        call disp%show( [lb, rand, ub] )
        call disp%show("call setAsserted(lb <= rand .and. rand <= ub)")
                        call setAsserted(lb <= rand .and. rand <= ub)
        end do
        call disp%skip()
    end block

    block
        use pm_kind, only: IKC => IKS
        integer(IKC) :: rand, lb, ub
        do itry = 1, ntry
        call disp%show("digits(0_IKC)")
        call disp%show( digits(0_IKC) )
        call disp%show("lb = -3_IKC; ub = 5_IKC")
                        lb = -3_IKC; ub = 5_IKC
        call disp%show("rand = getUnifRand(rng, lb, ub)")
                        rand = getUnifRand(rng, lb, ub)
        call disp%show("[lb, rand, ub]")
        call disp%show( [lb, rand, ub] )
        call disp%show("call setAsserted(lb <= rand .and. rand <= ub)")
                        call setAsserted(lb <= rand .and. rand <= ub)
        end do
        call disp%skip()
    end block

    block
        use pm_kind, only: IKC => IKL
        integer(IKC) :: rand, lb, ub
        do itry = 1, ntry
        call disp%show("digits(0_IKC)")
        call disp%show( digits(0_IKC) )
        call disp%show("lb = -huge(0_IKC); ub = huge(0_IKC) / 2_IKC")
                        lb = -huge(0_IKC); ub = huge(0_IKC) / 2_IKC
        call disp%show("rand = getUnifRand(rng, lb, ub)")
                        rand = getUnifRand(rng, lb, ub)
        call disp%show("[lb, rand, ub]")
        call disp%show( [lb, rand, ub] )
        call disp%show("call setAsserted(lb <= rand .and. rand <= ub)")
                        call setAsserted(lb <= rand .and. rand <= ub)
        end do
        call disp%skip()
    end block

    block
        character(2) :: rand, lb, ub
        do itry = 1, ntry
        call disp%show("lb = 'ai'; ub = 'by'")
                        lb = 'ai'; ub = 'by'
        call disp%show("rand = getUnifRand(rng, lb, ub)")
                        rand = getUnifRand(rng, lb, ub)
        call disp%show("[lb, rand, ub]")
        call disp%show( [lb, rand, ub] )
        call disp%show("call setAsserted(lb <= rand .and. rand <= ub)")
                        call setAsserted(lb <= rand .and. rand <= ub)
        end do
        call disp%skip()
    end block

    block
        use pm_logicalCompare, only: operator(<=)
        logical :: rand, lb, ub
        do itry = 1, 10
        call disp%show("lb = .false.; ub = .true.")
                        lb = .false.; ub = .true.
        call disp%show("rand = getUnifRand(rng, lb, ub)")
                        rand = getUnifRand(rng, lb, ub)
        call disp%show("[lb, rand, ub]")
        call disp%show( [lb, rand, ub] )
        call disp%show("call setAsserted(lb <= rand .and. rand <= ub)")
                        call setAsserted(lb <= rand .and. rand <= ub)
        end do
        call disp%skip()
    end block

    block
        use pm_complexCompareAll, only: operator(<=), operator(<)
        complex :: rand, lb, ub
        do itry = 1, ntry
        call disp%show("lb = (-1., +1.); ub = (1., +2.)")
                        lb = (-1., +1.); ub = (1., +2.)
        call disp%show("rand = getUnifRand(rng, lb, ub)")
                        rand = getUnifRand(rng, lb, ub)
        call disp%show("[lb, rand, ub]")
        call disp%show( [lb, rand, ub] )
        call disp%show("call setAsserted(lb <= rand .and. rand < ub)")
                        call setAsserted(lb <= rand .and. rand < ub)
        end do
        call disp%skip()
    end block

    block
        use pm_kind, only: RKC => RKH
        real(RKC) :: rand, lb, ub
        call disp%show("lb = 2._RKC; ub = lb + spacing(lb)")
                        lb = 2._RKC; ub = lb + spacing(lb)
        do itry = 1, ntry
        call disp%show("call setUnifRand(rng, rand, lb, ub)")
                        call setUnifRand(rng, rand, lb, ub)
        call disp%show("[lb, rand, ub], format = getFormat(width = 42_IK, ndigit = 35_IK)")
        call disp%show( [lb, rand, ub], format = getFormat(width = 42_IK, ndigit = 35_IK) )
        call disp%show("call setAsserted(lb <= rand .and. rand < ub)")
                        call setAsserted(lb <= rand .and. rand < ub)
        end do
        call disp%skip()
    end block

    block
        real :: rand, lb, ub
        do itry = 1, ntry
        call disp%show("lb = -3.; ub = 5.")
                        lb = -3.; ub = 5.
        call disp%show("rand = getUnifRand(rng, lb, ub)")
                        rand = getUnifRand(rng, lb, ub)
        call disp%show("[lb, rand, ub]")
        call disp%show( [lb, rand, ub] )
        call disp%show("call setAsserted(lb <= rand .and. rand < ub)")
                        call setAsserted(lb <= rand .and. rand < ub)
        end do
        call disp%skip()
    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example rand array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i
        integer(IK), parameter :: NP = 5000_IK
        integer, dimension(NP) :: rand
        call setUnifRand(rng, rand, -2, 3)
        open(newunit = fileUnit, file = "splitmix64_type.IK.txt")
        write(fileUnit,"(1(g0,:,' '))") rand
        close(fileUnit)
    end block

    block
        integer(IK) :: fileUnit, i
        integer(IK), parameter :: NP = 5000_IK
        complex, dimension(NP) :: rand
        call setUnifRand(rng, rand, (-2., +2.), (3., 5.))
        open(newunit = fileUnit, file = "splitmix64_type.CK.txt")
        write(fileUnit,"(2(g0,:,' '))") rand
        close(fileUnit)
    end block

    block
        integer(IK) :: fileUnit
        integer(IK), parameter :: NP = 5000_IK
        real, dimension(NP) :: rand
        call setUnifRand(rng, rand)
        open(newunit = fileUnit, file = "splitmix64_type.RK.txt")
        write(fileUnit,"(1(g0,:,' '))") rand
        close(fileUnit)
    end block

end program example