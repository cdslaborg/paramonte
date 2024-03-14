program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RKC => RK32 ! all real kinds are supported.
    use pm_distUnif, only: xoshiro256ssw_type
    use pm_distPois, only: setPoisRand
    use pm_arraySpace, only: setLinSpace
    use pm_arraySpace, only: setLogSpace
    use pm_io, only: display_type

    implicit none

    type(xoshiro256ssw_type) :: rng
    integer(IK), parameter  :: NP = 1000_IK
    integer(IK) :: rand(NP)
    real(RKC) :: lambda(NP)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call setLinSpace(lambda, x1 = 0.1_RKC, x2 = 100._RKC)

    call disp%skip()
    call disp%show("lambda(1)")
    call disp%show( lambda(1) )
    call disp%show("call setPoisRand(rand(1), exp(-lambda(1)))")
                    call setPoisRand(rand(1), exp(-lambda(1)))
    call disp%show("rand(1)")
    call disp%show( rand(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("lambda(1)")
    call disp%show( lambda(1) )
    call disp%show("rng = xoshiro256ssw_type()")
                    rng = xoshiro256ssw_type()
    call disp%show("call setPoisRand(rng, rand(1:2), exp(-lambda(1)))")
                    call setPoisRand(rng, rand(1:2), exp(-lambda(1)))
    call disp%show("rand(1:2)")
    call disp%show( rand(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("lambda(1)")
    call disp%show( lambda(1) )
    call disp%show("call setPoisRand(rand(1:2), exp(-lambda(1)))")
                    call setPoisRand(rand(1:2), exp(-lambda(1)))
    call disp%show("rand(1:2)")
    call disp%show( rand(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("lambda(NP-2:NP)")
    call disp%show( lambda(NP-2:NP) )
    call disp%show("call setPoisRand(rand(NP-2:NP), lambda(NP-2:NP), log(lambda(NP-2:NP)), sqrt(lambda(NP-2:NP)))")
                    call setPoisRand(rand(NP-2:NP), lambda(NP-2:NP), log(lambda(NP-2:NP)), sqrt(lambda(NP-2:NP)))
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
        integer(IK) :: rand(NP, 4)
        call setPoisRand(rand(:, 1), exp(-lambda(1)))
        call setPoisRand(rand(:, 2), exp(-lambda(2)))
        call setPoisRand(rand(:, 3), exp(-lambda(3)))
        call setPoisRand(rand(:, 4), lambda(4), log(lambda(4)), sqrt(lambda(4)))
        open(newunit = fileUnit, file = "setPoisRand.IK.txt")
        write(fileUnit,"(4(g0,:,' '))") (rand(i,:), i = 1, NP)
        close(fileUnit)
    end block

end program example