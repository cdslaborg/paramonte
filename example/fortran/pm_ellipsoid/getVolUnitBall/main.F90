program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_ellipsoid, only: getVolUnitBall

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getVolUnitBall(ndim = 2.)")
    call disp%show( getVolUnitBall(ndim = 2.) )
    call disp%skip()

    call disp%skip()
    call disp%show("getVolUnitBall(ndim = [real :: 0, 1, 2, 3, 4])")
    call disp%show( getVolUnitBall(ndim = [real :: 0, 1, 2, 3, 4]) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_kind, only: RKD
        use pm_arrayRange, only: getRange
        integer(IK), parameter :: nsim = 30_IK
        real(RKD) :: volUnitBall(0 : nsim), ndim(0 : nsim)
        integer(IK) :: fileUnit, i
        open(newunit = fileUnit, file = "getVolUnitBall.RK.txt")
        ndim(0 : nsim) = getRange(0_IK, nsim)
        volUnitBall = getVolUnitBall(ndim)
        do i = 0, nsim
            write(fileUnit, "(*(g0,:,','))") ndim(i), volUnitBall(i)
        end do
        close(fileUnit)
    end block

end program example