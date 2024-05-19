program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RKG => RKS ! all real kinds are supported.
    use pm_distUnif, only: xoshiro256ssw_type
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: setLogSpace
    use pm_distGeom, only: setGeomRand
    use pm_io, only: display_type

    implicit none

    type(xoshiro256ssw_type) :: rng
    integer(IK), parameter  :: NP = 1000_IK
    integer(IK) :: rand(NP)
    real(RKG) :: logProbFailure(NP)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    logProbFailure = log(1 - getLinSpace(x1 = 0.001_RKG, x2 = .999_RKG, count = NP))

    call disp%skip()
    call disp%show("logProbFailure(1)")
    call disp%show( logProbFailure(1) )
    call disp%show("call setGeomRand(rand(1), logProbFailure(1))")
                    call setGeomRand(rand(1), logProbFailure(1))
    call disp%show("rand(1)")
    call disp%show( rand(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logProbFailure(1)")
    call disp%show( logProbFailure(1) )
    call disp%show("rng = xoshiro256ssw_type()")
                    rng = xoshiro256ssw_type()
    call disp%show("call setGeomRand(rng, rand(1:2), logProbFailure(1))")
                    call setGeomRand(rng, rand(1:2), logProbFailure(1))
    call disp%show("rand(1:2)")
    call disp%show( rand(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("logProbFailure(1)")
    call disp%show( logProbFailure(1) )
    call disp%show("call setGeomRand(rand(1:2), logProbFailure(1))")
                    call setGeomRand(rand(1:2), logProbFailure(1))
    call disp%show("rand(1:2)")
    call disp%show( rand(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("logProbFailure(NP-2:NP)")
    call disp%show( logProbFailure(NP-2:NP) )
    call disp%show("call setGeomRand(rand(NP-2:NP), logProbFailure(NP-2:NP))")
                    call setGeomRand(rand(NP-2:NP), logProbFailure(NP-2:NP))
    call disp%show("rand(NP-2:NP)")
    call disp%show( rand(NP-2:NP) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example rand array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i
        integer(IK) , parameter :: NP = 5000_IK
        real(RKG)   , parameter :: probSuccess(4) = [.1_RKG, .2_RKG, .5_RKG, .8_RKG]
        open(newunit = fileUnit, file = "setGeomRand.IK.txt")
        do i = 1, NP
            call setGeomRand(rand(1:4), log(1 - probSuccess))
            write(fileUnit,"(*(g0,:,','))") rand(1:4)
        end do
        close(fileUnit)
    end block

end program example