program example

    use pm_kind, only: SK, IK, LK, RKG => RKH
    use pm_ellipsoid, only: setVolUnitBall
    use pm_io, only: display_type

    implicit none

    real(RKG) :: volUnitBall(5)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("call setVolUnitBall(volUnitBall(1), ndim = 2_IK)")
                    call setVolUnitBall(volUnitBall(1), ndim = 2_IK)
    call disp%show("volUnitBall(1)")
    call disp%show( volUnitBall(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setVolUnitBall(volUnitBall(1:5), ndim = [integer(IK) :: 0, 1, 2, 3, 4])")
                    call setVolUnitBall(volUnitBall(1:5), ndim = [integer(IK) :: 0, 1, 2, 3, 4])
    call disp%show("volUnitBall(1:5)")
    call disp%show( volUnitBall(1:5) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_kind, only: RKD
        use pm_arrayRange, only: getRange
        integer(IK), parameter :: nsim = 30_IK
        integer(IK) :: ndim(0 : nsim)
        real(RKD) :: volUnitBall(0 : nsim)
        integer(IK) :: fileUnit, i
        open(newunit = fileUnit, file = "setVolUnitBall.RK.txt")
        ndim(0 : nsim) = getRange(0_IK, nsim)
        call setVolUnitBall(volUnitBall, ndim)
        do i = 0, nsim
            write(fileUnit, "(*(g0,:,','))") ndim(i), volUnitBall(i)
        end do
        close(fileUnit)
    end block

end program example