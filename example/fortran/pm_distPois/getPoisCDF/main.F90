program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKG => RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distPois, only: getPoisCDF

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getPoisCDF(0_IK, lambda = 2._RKG)")
    call disp%show( getPoisCDF(0_IK, lambda = 2._RKG) )
    call disp%skip()

    call disp%skip()
    call disp%show("getPoisCDF(1_IK, lambda = 2._RKG)")
    call disp%show( getPoisCDF(1_IK, lambda = 2._RKG) )
    call disp%skip()

    call disp%skip()
    call disp%show("getPoisCDF(2_IK, lambda = 2._RKG)")
    call disp%show( getPoisCDF(2_IK, lambda = 2._RKG) )
    call disp%skip()

    call disp%skip()
    call disp%show("getPoisCDF([integer(IK) :: 0, 1, 2], lambda = 2._RKG)")
    call disp%show( getPoisCDF([integer(IK) :: 0, 1, 2], lambda = 2._RKG) )
    call disp%skip()

    call disp%skip()
    call disp%show("getPoisCDF([integer(IK) :: 0, 1, 2], lambda = [0.1_RKG, 1._RKG, 10._RKG])")
    call disp%show( getPoisCDF([integer(IK) :: 0, 1, 2], lambda = [0.1_RKG, 1._RKG, 10._RKG]) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arrayRange, only: getRange
        integer(IK), allocatable :: count(:)
        integer :: fileUnit, i

        count = getRange(0_IK, 20_IK)
        open(newunit = fileUnit, file = "getPoisCDF.IK.txt")
        do i = 1, size(count)
            write(fileUnit, "(*(g0,:,' '))") count(i), exp(getPoisCDF(count(i), [.1_RKG, 1._RKG, 4._RKG, 10._RKG]))
        end do
        close(fileUnit)

    end block

end program example