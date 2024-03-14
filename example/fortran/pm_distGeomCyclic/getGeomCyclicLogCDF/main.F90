program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKC => RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distGeomCyclic, only: getGeomCyclicLogCDF

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getGeomCyclicLogCDF(1_IK, probSuccess = .5_RKC, period = 2_IK)")
    call disp%show( getGeomCyclicLogCDF(1_IK, probSuccess = .5_RKC, period = 2_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGeomCyclicLogCDF(2_IK, probSuccess = .5_RKC, period = 2_IK)")
    call disp%show( getGeomCyclicLogCDF(2_IK, probSuccess = .5_RKC, period = 2_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGeomCyclicLogCDF(2_IK, probSuccess = .2_RKC, period = 2_IK)")
    call disp%show( getGeomCyclicLogCDF(2_IK, probSuccess = .2_RKC, period = 2_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGeomCyclicLogCDF([integer(IK) :: 1, 2, 3], probSuccess = .9_RKC, period = 3_IK)")
    call disp%show( getGeomCyclicLogCDF([integer(IK) :: 1, 2, 3], probSuccess = .9_RKC, period = 3_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGeomCyclicLogCDF([integer(IK) :: 1, 2, 3, 4], probSuccess = [0.01_RKC, 0.1_RKC, .1_RKC, 1._RKC], period = 10_IK)")
    call disp%show( getGeomCyclicLogCDF([integer(IK) :: 1, 2, 3, 4], probSuccess = [0.01_RKC, 0.1_RKC, .1_RKC, 1._RKC], period = 10_IK) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arrayRange, only: getRange
        integer(IK), allocatable :: count(:)
        integer :: fileUnit, i

        count = getRange(1_IK, 10_IK)
        open(newunit = fileUnit, file = "getGeomCyclicLogCDF.IK.txt")
        do i = 1, size(count)
            write(fileUnit, "(*(g0,:,','))" ) count(i), exp(getGeomCyclicLogCDF(count(i), [.05_RKC, .25_RKC, .05_RKC, .25_RKC], period = [10_IK, 10_IK, 10000_IK, 10000_IK]))
        end do
        close(fileUnit)

    end block

end program example