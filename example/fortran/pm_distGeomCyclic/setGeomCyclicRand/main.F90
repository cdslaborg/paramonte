program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RKG => RKS ! all real kinds are supported.
    use pm_distUnif, only: xoshiro256ssw_type
    use pm_distGeomCyclic, only: setGeomCyclicRand
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: setLogSpace
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
    call disp%show("call setGeomCyclicRand(rand(1), logProbFailure(1), period = 10_IK)")
                    call setGeomCyclicRand(rand(1), logProbFailure(1), period = 10_IK)
    call disp%show("rand(1)")
    call disp%show( rand(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logProbFailure(1)")
    call disp%show( logProbFailure(1) )
    call disp%show("rng = xoshiro256ssw_type()")
                    rng = xoshiro256ssw_type()
    call disp%show("call setGeomCyclicRand(rng, rand(1:2), logProbFailure(1), period = 10_IK)")
                    call setGeomCyclicRand(rng, rand(1:2), logProbFailure(1), period = 10_IK)
    call disp%show("rand(1:2)")
    call disp%show( rand(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("logProbFailure(1)")
    call disp%show( logProbFailure(1) )
    call disp%show("call setGeomCyclicRand(rand(1:2), logProbFailure(1), period = 100_IK)")
                    call setGeomCyclicRand(rand(1:2), logProbFailure(1), period = 100_IK)
    call disp%show("rand(1:2)")
    call disp%show( rand(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("logProbFailure(NP-2:NP)")
    call disp%show( logProbFailure(NP-2:NP) )
    call disp%show("call setGeomCyclicRand(rand(NP-2:NP), logProbFailure(NP-2:NP), period = 100_IK)")
                    call setGeomCyclicRand(rand(NP-2:NP), logProbFailure(NP-2:NP), period = 100_IK)
    call disp%show("rand(NP-2:NP)")
    call disp%show( rand(NP-2:NP) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example rand array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        integer :: fileUnit, i
        integer(IK) :: rand(4)
        open(newunit = fileUnit, file = "setGeomCyclicRand.IK.txt")
        do i = 1, 5000
            call setGeomCyclicRand(rand, log(1 - [.05_RKG, .25_RKG, .05_RKG, .25_RKG]), period = [10_IK, 10_IK, 10000_IK, 10000_IK])
            write(fileUnit, "(*(g0,:,','))") rand
        end do
        close(fileUnit)

    end block

end program example