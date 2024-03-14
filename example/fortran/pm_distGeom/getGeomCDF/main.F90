program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKC => RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distGeom, only: getGeomCDF

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getGeomCDF(1_IK, probSuccess = .5_RKC)")
    call disp%show( getGeomCDF(1_IK, probSuccess = .5_RKC) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGeomCDF(2_IK, probSuccess = .5_RKC)")
    call disp%show( getGeomCDF(2_IK, probSuccess = .5_RKC) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGeomCDF(2_IK, probSuccess = .2_RKC)")
    call disp%show( getGeomCDF(2_IK, probSuccess = .2_RKC) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGeomCDF([integer(IK) :: 1, 2, 3], probSuccess = .9_RKC)")
    call disp%show( getGeomCDF([integer(IK) :: 1, 2, 3], probSuccess = .9_RKC) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGeomCDF([integer(IK) :: 1, 2, 3, 4], probSuccess = [0.01_RKC, 0.1_RKC, .1_RKC, 1._RKC])")
    call disp%show( getGeomCDF([integer(IK) :: 1, 2, 3, 4], probSuccess = [0.01_RKC, 0.1_RKC, .1_RKC, 1._RKC]) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arrayRange, only: getRange
        integer(IK), allocatable :: probSuccess(:)
        integer :: fileUnit, i

        probSuccess = getRange(1_IK, 10_IK)
        open(newunit = fileUnit, file = "getGeomCDF.IK.txt")
        do i = 1, size(probSuccess)
            write(fileUnit, "(*(g0,:,' '))" ) probSuccess(i), getGeomCDF(probSuccess(i), [.1_RKC, .2_RKC, .5_RKC, .8_RKC])
        end do
        close(fileUnit)

    end block

end program example