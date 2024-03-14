program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RKC => RK32 ! all real kinds are supported.
    use pm_distPois, only: getPoisRand
    use pm_arraySpace, only: setLinSpace
    use pm_arraySpace, only: setLogSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 1000_IK
    integer(IK) :: rand(NP)
    real(RKC) :: lambda(NP)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call setLinSpace(lambda, x1 = 0.1_RKC, x2 = 100._RKC)

    call disp%skip()
    call disp%show("lambda(1)")
    call disp%show( lambda(1) )
    call disp%show("rand(1) = getPoisRand(lambda(1))")
                    rand(1) = getPoisRand(lambda(1))
    call disp%show("rand(1)")
    call disp%show( rand(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("lambda(1)")
    call disp%show( lambda(1) )
    call disp%show("rand(1:2) = getPoisRand(lambda(1))")
                    rand(1:2) = getPoisRand(lambda(1))
    call disp%show("rand(1:2)")
    call disp%show( rand(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("lambda(1)")
    call disp%show( lambda(1) )
    call disp%show("rand(1:2) = getPoisRand(lambda(1))")
                    rand(1:2) = getPoisRand(lambda(1))
    call disp%show("rand(1:2)")
    call disp%show( rand(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("lambda(NP-2:NP)")
    call disp%show( lambda(NP-2:NP) )
    call disp%show("rand(NP-2:NP) = getPoisRand(lambda(NP-2:NP))")
                    rand(NP-2:NP) = getPoisRand(lambda(NP-2:NP))
    call disp%show("rand(NP-2:NP)")
    call disp%show( rand(NP-2:NP) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example rand array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i
        integer(IK) , parameter :: NP = 5000_IK
        real(RKC)   , parameter :: lambda(4) = [.1_RKC, 1._RKC, 4._RKC, 11._RKC]
        open(newunit = fileUnit, file = "getPoisRand.IK.txt")
        do i = 1, NP
            write(fileUnit,"(*(g0,:,' '))") getPoisRand(lambda)
        end do
        close(fileUnit)
    end block

end program example