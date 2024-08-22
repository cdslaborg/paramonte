program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKG => RKH ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_mathBeta, only: setBetaInc
    use pm_mathBeta, only: getLogBeta
    use pm_err, only: setAsserted

    implicit none

    integer(IK) :: info(3)
    real(RKG)   :: betaInc(3), alpha, beta, reltol
    real(RKG)   , allocatable :: alphas(:), betas(:)
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized (lower) Incomplete Beta Function using its series representation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("alpha = 2.; beta = 3.")
                    alpha = 2.; beta = 3.
    call disp%show("call setBetaInc(betaInc(1), 0._RKG, alpha, beta, getLogBeta(alpha, beta), signed = .false._LK, info = info(1))")
                    call setBetaInc(betaInc(1), 0._RKG, alpha, beta, getLogBeta(alpha, beta), signed = .false._LK, info = info(1))
    call disp%show("call setAsserted(.not. info(1) < 0)")
                    call setAsserted(.not. info(1) < 0)
    call disp%show("betaInc(1)")
    call disp%show( betaInc(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("alpha = 2.; beta = 3.")
                    alpha = 2.; beta = 3.
    call disp%show("call setBetaInc(betaInc(1), .5_RKG, alpha, beta, getLogBeta(alpha, beta), signed = .false._LK, info = info(1)) ! 0.6875")
                    call setBetaInc(betaInc(1), .5_RKG, alpha, beta, getLogBeta(alpha, beta), signed = .false._LK, info = info(1))
    call disp%show("call setAsserted(.not. info(1) < 0)")
                    call setAsserted(.not. info(1) < 0)
    call disp%show("betaInc(1)")
    call disp%show( betaInc(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("alpha = 2.; beta = 3.")
                    alpha = 2.; beta = 3.
    call disp%show("call setBetaInc(betaInc(1), 1._RKG, alpha, beta, getLogBeta(alpha, beta), signed = .false._LK, info = info(1))")
                    call setBetaInc(betaInc(1), 1._RKG, alpha, beta, getLogBeta(alpha, beta), signed = .false._LK, info = info(1))
    call disp%show("call setAsserted(.not. info(1) < 0)")
                    call setAsserted(.not. info(1) < 0)
    call disp%show("betaInc(1)")
    call disp%show( betaInc(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized (lower) Incomplete Beta Function for a vector of points.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("alpha = 2.; beta = 3.")
                    alpha = 2.; beta = 3.
    call disp%show("call setBetaInc(betaInc, [0._RKG, 0.5_RKG, 1._RKG], alpha, beta, getLogBeta(alpha, beta), signed = .false._LK, info = info)")
                    call setBetaInc(betaInc, [0._RKG, 0.5_RKG, 1._RKG], alpha, beta, getLogBeta(alpha, beta), signed = .false._LK, info = info)
    call disp%show("info")
    call disp%show( info )
    call disp%show("betaInc")
    call disp%show( betaInc )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized Incomplete Beta Function for a vector of input parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("alphas = [real(RKG) :: 0.1, 1., 10.]; beta = 3.")
                    alphas = [real(RKG) :: 0.1, 1., 10.]; beta = 3.
    call disp%show("call setBetaInc(betaInc, [0._RKG, 0.5_RKG, 1._RKG], alphas, beta, getLogBeta(alphas, beta), signed = .false._LK, info = info)")
                    call setBetaInc(betaInc, [0._RKG, 0.5_RKG, 1._RKG], alphas, beta, getLogBeta(alphas, beta), signed = .false._LK, info = info)
    call disp%show("info")
    call disp%show( info )
    call disp%show("betaInc")
    call disp%show( betaInc )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized (lower) Incomplete Beta Function using brute-force integration.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("reltol = epsilon(reltol)**.7")
                    reltol = epsilon(reltol)**.7
    call disp%show("alphas = 5000.**[1, -1]; betas = 5000.**[1, -1]")
                    alphas = 5000.**[1, -1]; betas = 5000.**[1, -1]
    call disp%show("call setBetaInc(betaInc(1:2), .5_RKG, alphas, betas, getLogBeta(alphas, betas), reltol, signed = .false._LK, info = info(1:2))")
                    call setBetaInc(betaInc(1:2), .5_RKG, alphas, betas, getLogBeta(alphas, betas), reltol, signed = .false._LK, info = info(1:2))
    call disp%show("info(1:2)")
    call disp%show( info(1:2) )
    call disp%show("betaInc(1:2)")
    call disp%show( betaInc(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("reltol = epsilon(reltol)**.7")
                    reltol = epsilon(reltol)**.7
    call disp%show("reltol")
    call disp%show( reltol )
    call disp%show("alphas = 5000.**[1, -1]; betas = 5000.**[1, -1]")
                    alphas = 5000.**[1, -1]; betas = 5000.**[1, -1]
    call disp%show("call setBetaInc(betaInc(1:2), .5_RKG, alphas, betas, getLogBeta(alphas, betas), signed = .false._LK, info = info(1:2))")
                    call setBetaInc(betaInc(1:2), .5_RKG, alphas, betas, getLogBeta(alphas, betas), signed = .false._LK, info = info(1:2))
    call disp%show("info(1:2)")
    call disp%show( info(1:2) )
    call disp%show("betaInc(1:2)")
    call disp%show( betaInc(1:2) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array of the regularized Incomplete Beta function for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !block
    !    use pm_arraySpace, only: setLinSpace
    !    integer(IK) , parameter :: n1000_IK
    !    real(RK) :: betaInc(3), x(nx)
    !    integer :: fileUnit, i
    !    call setLinSpace(x, 0._RKG, 1._RKG)
    !    open(newunit = fileUnit, file = "setBetaInc.RK.txt")
    !    do i = 1, nx
    !        call setBetaInc(betaInc, x(i), 5._RKG, [1.0_RKG, 5.0_RKG, 10._RKG], getLogBeta(5._RKG, [1.0_RKG, 5.0_RKG, 10._RKG]), err)
    !        write(fileUnit, "(*(g0,:,' '))") x_RK(i), betaInc
    !    end do
    !    close(fileUnit)
    !end block

    block
        use pm_kind, only: RKG => RKH, RKH
        use pm_arraySpace, only: setLinSpace
        integer(IK) , parameter :: nx = 1000
        real(RKG)   , parameter :: alphas(*) = [0.1_RKG, 10._RKG, 1._RKG, 0.1_RKG, 10._RKG], betas(*) = [0.1_RKG, 0.1_RKG, 1._RKG, 10._RKG, 10._RKG]
        real(RKG)   :: betaInc(max(size(alphas), size(betas)))
        integer(IK) :: fileUnit, i, info(size(betaInc))
        real(RKG)   :: x(nx)
        call setLinSpace(x, 0._RKG, 1._RKG, fopen = .true._LK, lopen = .true._LK)
        open(newunit = fileUnit, file = "setBetaInc.RK.txt")
        do i = 1, nx
            call setBetaInc(betaInc, x(i), alphas, betas, getLogBeta(alphas, betas), signed = .false._LK, info = info)
            if (any(info /= 0)) error stop
            write(fileUnit, "(*(g0,:,','))") x(i), betaInc
        end do
        close(fileUnit)
    end block

end program example