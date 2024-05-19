program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RKG => RKS ! all real kinds are supported.
    use pm_distGeom, only: getGeomRand
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
    call disp%show("rand(1) = getGeomRand(probSuccess(1))")
                    rand(1) = getGeomRand(probSuccess(1))
    call disp%show("rand(1)")
    call disp%show( rand(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("probSuccess(1)")
    call disp%show( probSuccess(1) )
    call disp%show("rand(1:2) = getGeomRand(probSuccess(1))")
                    rand(1:2) = getGeomRand(probSuccess(1))
    call disp%show("rand(1:2)")
    call disp%show( rand(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("probSuccess(1)")
    call disp%show( probSuccess(1) )
    call disp%show("rand(1:2) = getGeomRand(probSuccess(1))")
                    rand(1:2) = getGeomRand(probSuccess(1))
    call disp%show("rand(1:2)")
    call disp%show( rand(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("probSuccess(NP-2:NP)")
    call disp%show( probSuccess(NP-2:NP) )
    call disp%show("rand(NP-2:NP) = getGeomRand(probSuccess(NP-2:NP))")
                    rand(NP-2:NP) = getGeomRand(probSuccess(NP-2:NP))
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
        open(newunit = fileUnit, file = "getGeomRand.IK.txt")
        do i = 1, NP
            write(fileUnit,"(*(g0,:,','))") getGeomRand(probSuccess)
        end do
        close(fileUnit)
    end block

end program example