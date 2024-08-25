program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RKG => RKS ! all real kinds are supported.
    use pm_distUnif, only: xoshiro256ssw_type
    use pm_distBeta, only: setBetaRand
    use pm_arraySpace, only: setLinSpace
    use pm_arraySpace, only: setLogSpace
    use pm_io, only: display_type

    implicit none

    type(xoshiro256ssw_type) :: rng
    integer(IK), parameter  :: NP = 1000_IK
    real(RKG), dimension(NP) :: alpha, beta, rand

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call setLogSpace(alpha, logx1 = log(0.1_RKG), logx2 = log(10._RKG))
    call setLogSpace(beta, logx1 = log(0.1_RKG), logx2 = log(10._RKG))

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate random numbers from the Beta distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("alpha(1)")
    call disp%show( alpha(1) )
    call disp%show("beta(1)")
    call disp%show( beta(1) )
    call disp%show("call setBetaRand(rand(1), 1._RKG, beta = beta(1))")
                    call setBetaRand(rand(1), 1._RKG, beta = beta(1))
    call disp%show("rand(1)")
    call disp%show( rand(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("alpha(1)")
    call disp%show( alpha(1) )
    call disp%show("beta(1)")
    call disp%show( beta(1) )
    call disp%show("rng = xoshiro256ssw_type()")
                    rng = xoshiro256ssw_type()
    call disp%show("call setBetaRand(rng, rand(1:2), 1._RKG, beta = beta(1))")
                    call setBetaRand(rng, rand(1:2), 1._RKG, beta = beta(1))
    call disp%show("rand(1:2)")
    call disp%show( rand(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("alpha(1)")
    call disp%show( alpha(1) )
    call disp%show("beta(1)")
    call disp%show( beta(1) )
    call disp%show("call setBetaRand(rand(1:2), alpha(1), beta = beta(1))")
                    call setBetaRand(rand(1:2), alpha(1), beta = beta(1))
    call disp%show("rand(1:2)")
    call disp%show( rand(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("alpha(1:NP:NP/3)")
    call disp%show( alpha(1:NP:NP/3) )
    call disp%show("beta(1:NP:NP/3)")
    call disp%show( beta(1:NP:NP/3) )
    call disp%show("call setBetaRand(rand(1:NP:NP/3), alpha(1:NP:NP/3), beta = beta(1:NP:NP/3))")
                    call setBetaRand(rand(1:NP:NP/3), alpha(1:NP:NP/3), beta = beta(1:NP:NP/3))
    call disp%show("rand(1:NP:NP/3)")
    call disp%show( rand(1:NP:NP/3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example rand array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_io, only: getErrTableWrite
        real(RKG) :: rand(5000, 4)
        call setBetaRand(rand(:, 1), alpha = 0.5_RKG, beta = 0.5_RKG)
        call setBetaRand(rand(:, 2), alpha = 2.0_RKG, beta = 2.0_RKG)
        call setBetaRand(rand(:, 3), alpha = 2.0_RKG, beta = 5.0_RKG)
        call setBetaRand(rand(:, 4), alpha = 5.0_RKG, beta = 2.0_RKG)
        if (0 /= getErrTableWrite(SK_"setBetaRand.RK.txt", rand)) error stop "Failed to write table."
    end block

end program example