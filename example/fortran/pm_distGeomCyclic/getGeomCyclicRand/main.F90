program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RKG => RKS ! all real kinds are supported.
    use pm_distGeomCyclic, only: getGeomCyclicRand
    use pm_arraySpace, only: setLinSpace
    use pm_arraySpace, only: setLogSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 1000_IK
    integer(IK) :: rand(NP)
    real(RKG) :: probSuccess(NP)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call setLinSpace(probSuccess, x1 = 0.001_RKG, x2 = 1._RKG)

    call disp%skip()
    call disp%show("probSuccess(1)")
    call disp%show( probSuccess(1) )
    call disp%show("rand(1) = getGeomCyclicRand(probSuccess(1), period = 10_IK)")
                    rand(1) = getGeomCyclicRand(probSuccess(1), period = 10_IK)
    call disp%show("rand(1)")
    call disp%show( rand(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("probSuccess(1)")
    call disp%show( probSuccess(1) )
    call disp%show("rand(1:2) = getGeomCyclicRand(probSuccess(1), period = 10_IK)")
                    rand(1:2) = getGeomCyclicRand(probSuccess(1), period = 10_IK)
    call disp%show("rand(1:2)")
    call disp%show( rand(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("probSuccess(1)")
    call disp%show( probSuccess(1) )
    call disp%show("rand(1:2) = getGeomCyclicRand(probSuccess(1), period = 100_IK)")
                    rand(1:2) = getGeomCyclicRand(probSuccess(1), period = 100_IK)
    call disp%show("rand(1:2)")
    call disp%show( rand(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("probSuccess(NP-2:NP)")
    call disp%show( probSuccess(NP-2:NP) )
    call disp%show("rand(NP-2:NP) = getGeomCyclicRand(probSuccess(NP-2:NP), period = 100_IK)")
                    rand(NP-2:NP) = getGeomCyclicRand(probSuccess(NP-2:NP), period = 100_IK)
    call disp%show("rand(NP-2:NP)")
    call disp%show( rand(NP-2:NP) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example rand array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: rand(4)
        integer :: fileUnit, i
        open(newunit = fileUnit, file = "getGeomCyclicRand.IK.txt")
        do i = 1, 5000
            rand = getGeomCyclicRand([.05_RKG, .25_RKG, .05_RKG, .25_RKG], period = [10_IK, 10_IK, 10000_IK, 10000_IK])
            write(fileUnit, "(*(g0,:,','))") rand
        end do
        close(fileUnit)
    end block

end program example