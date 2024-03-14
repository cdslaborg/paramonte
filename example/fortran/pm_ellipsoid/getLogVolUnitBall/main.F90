program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_ellipsoid, only: getLogVolUnitBall

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getLogVolUnitBall(ndim = 2.)")
    call disp%show( getLogVolUnitBall(ndim = 2.) )
    call disp%skip()

    call disp%skip()
    call disp%show("getLogVolUnitBall(ndim = [real :: 0, 1, 2, 3, 4])")
    call disp%show( getLogVolUnitBall(ndim = [real :: 0, 1, 2, 3, 4]) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLinSpace
        integer(IK), parameter :: nsim = 1000_IK
        real :: volUnitBall(nsim), ndim(nsim)
        integer(IK) :: fileUnit, i
        open(newunit = fileUnit, file = "getLogVolUnitBall.RK.txt")
        call setLinSpace(ndim, 0., 20.)
        volUnitBall = exp(getLogVolUnitBall(ndim))
        do i = 1, nsim
            write(fileUnit, "(*(g0,:,','))") ndim(i), volUnitBall(i)
        end do
        close(fileUnit)
    end block

end program example