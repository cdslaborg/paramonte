program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKG => RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distGeom, only: getGeomLogPMF

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getGeomLogPMF(1_IK, probSuccess = .5_RKG)")
    call disp%show( getGeomLogPMF(1_IK, probSuccess = .5_RKG) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGeomLogPMF(2_IK, probSuccess = .5_RKG)")
    call disp%show( getGeomLogPMF(2_IK, probSuccess = .5_RKG) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGeomLogPMF(2_IK, probSuccess = .2_RKG)")
    call disp%show( getGeomLogPMF(2_IK, probSuccess = .2_RKG) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGeomLogPMF([integer(IK) :: 1, 2, 3], probSuccess = .9_RKG)")
    call disp%show( getGeomLogPMF([integer(IK) :: 1, 2, 3], probSuccess = .9_RKG) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGeomLogPMF([integer(IK) :: 1, 2, 3, 4], probSuccess = [0.01_RKG, 0.1_RKG, .1_RKG, 1._RKG])")
    call disp%show( getGeomLogPMF([integer(IK) :: 1, 2, 3, 4], probSuccess = [0.01_RKG, 0.1_RKG, .1_RKG, 1._RKG]) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arrayRange, only: getRange
        integer(IK), allocatable :: count(:)
        integer :: fileUnit, i

        count = getRange(1_IK, 10_IK)
        open(newunit = fileUnit, file = "getGeomLogPMF.IK.txt")
        do i = 1, size(count)
            write(fileUnit, "(*(g0,:,' '))") count(i), exp(getGeomLogPMF(count(i), [.1_RKG, .2_RKG, .5_RKG, .8_RKG]))
        end do
        close(fileUnit)

    end block

end program example