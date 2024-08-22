program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_mathBeta, only: getBetaInc

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized (lower) Incomplete Beta Function using its series representation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getBetaInc(x = 0._RK, alpha = 2._RK, beta = 3._RK)")
    call disp%show( getBetaInc(x = 0._RK, alpha = 2._RK, beta = 3._RK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getBetaInc(x = .5_RK, alpha = 2._RK, beta = 3._RK) ! 0.6875")
    call disp%show( getBetaInc(x = .5_RK, alpha = 2._RK, beta = 3._RK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getBetaInc(x = 1._RK, alpha = 2._RK, beta = 3._RK)")
    call disp%show( getBetaInc(x = 1._RK, alpha = 2._RK, beta = 3._RK) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized (lower) Incomplete Beta Function for a vector of points.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getBetaInc(x = [0._RK, 0.5_RK, 1._RK], alpha = 2._RK, beta = 3._RK)")
    call disp%show( getBetaInc(x = [0._RK, 0.5_RK, 1._RK], alpha = 2._RK, beta = 3._RK) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized Incomplete Beta Function for a vector of input parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getBetaInc(x = [0._RK, 0.5_RK, 1._RK], alpha = [0.1_RK, 1._RK, 10._RK], beta = 3._RK) ! 0, 0.875, 1")
    call disp%show( getBetaInc(x = [0._RK, 0.5_RK, 1._RK], alpha = [0.1_RK, 1._RK, 10._RK], beta = 3._RK) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array of the regularized Incomplete Beta function for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !block
    !    use pm_arraySpace, only: setLinSpace
    !    integer(IK) , parameter :: nx = 1000_IK
    !    real(RK) :: x(nx)
    !    integer :: fileUnit, i
    !    call setLinSpace(x, 0._RK, 1._RK)
    !    open(newunit = fileUnit, file = "getBetaInc.RK.txt")
    !    do i = 1, nx
    !        write(fileUnit, "(*(g0,:,' '))") x_RK(i), getBetaInc(x(i), 5._RK, [1.0_RK, 5.0_RK, 10._RK])
    !    end do
    !    close(fileUnit)
    !end block

    block
        use pm_kind, only: RKG => RKH, RKH
        use pm_arraySpace, only: setLinSpace
        integer(IK) , parameter :: nx = 1000
        real(RKG)   , parameter :: alpha(*) = [0.1_RKG, 10._RKG, 1._RKG, 0.1_RKG, 10._RKG], beta(*) = [0.1_RKG, 0.1_RKG, 1._RKG, 10._RKG, 10._RKG]
        real(RKG)   :: betaInc(max(size(alpha), size(beta)))
        real(RKG)   :: x(nx)
        integer     :: fileUnit, i
        call setLinSpace(x, 0._RKG, 1._RKG, fopen = .true._LK, lopen = .true._LK)
        open(newunit = fileUnit, file = "getBetaInc.RK.txt")
        do i = 1, nx
            betaInc = getBetaInc(x(i), alpha, beta)
            write(fileUnit, "(*(g0,:,','))") x(i), merge(1 + betaInc, betaInc, betaInc < 0)
        end do
        close(fileUnit)
    end block

end program example