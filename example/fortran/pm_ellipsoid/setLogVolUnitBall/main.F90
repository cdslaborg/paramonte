program example

    use pm_kind, only: SK, IK, LK, RKG => RKH
    use pm_io, only: display_type
    use pm_ellipsoid, only: setLogVolUnitBall

    implicit none

    real(RKG) :: logVolUnitBall(5)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("call setLogVolUnitBall(logVolUnitBall(1), ndim = 2_IK)")
                    call setLogVolUnitBall(logVolUnitBall(1), ndim = 2_IK)
    call disp%show("logVolUnitBall(1)")
    call disp%show( logVolUnitBall(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogVolUnitBall(logVolUnitBall(1), ndim = 2._RKG)")
                    call setLogVolUnitBall(logVolUnitBall(1), ndim = 2._RKG)
    call disp%show("logVolUnitBall(1)")
    call disp%show( logVolUnitBall(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogVolUnitBall(logVolUnitBall(1:5), ndim = [integer(IK) :: 0, 1, 2, 3, 4])")
                    call setLogVolUnitBall(logVolUnitBall(1:5), ndim = [integer(IK) :: 0, 1, 2, 3, 4])
    call disp%show("logVolUnitBall(1:5)")
    call disp%show( logVolUnitBall(1:5) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogVolUnitBall(logVolUnitBall(1:5), ndim = [real(RKG) :: 0, 1, 2, 3, 4.5])")
                    call setLogVolUnitBall(logVolUnitBall(1:5), ndim = [real(RKG) :: 0, 1, 2, 3, 4.5])
    call disp%show("logVolUnitBall(1:5)")
    call disp%show( logVolUnitBall(1:5) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLinSpace
        integer(IK), parameter :: nsim = 1000_IK
        integer(IK) :: fileUnit, i
        real :: volUnitBall(nsim), ndim(nsim)
        open(newunit = fileUnit, file = "setLogVolUnitBall.RK.txt")
        call setLinSpace(ndim, 0., 20.)
        call setLogVolUnitBall(volUnitBall, ndim)
        do i = 1, nsim
            write(fileUnit, "(*(g0,:,','))") ndim(i), exp(volUnitBall(i))
        end do
        close(fileUnit)
    end block

end program example