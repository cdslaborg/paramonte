program example

    use pm_kind, only: SK
    use pm_kind, only: IK
    use pm_kind, only: LK
    use pm_kind, only: RKS, RKD, RKH
    use pm_io, only: display_type
    use pm_mathGammaNR, only: getGammaIncUppNR

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized Upper Incomplete Gamma Function using its series representation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUppNR(x = 1.5_RKS, kappa = 2._RKS)")
    call disp%show( getGammaIncUppNR(x = 1.5_RKS, kappa = 2._RKS) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUppNR(x = 1.5_RKD, kappa = 2._RKD)")
    call disp%show( getGammaIncUppNR(x = 1.5_RKD, kappa = 2._RKD) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUppNR(x = 1.5_RKH, kappa = 2._RKH)")
    call disp%show( getGammaIncUppNR(x = 1.5_RKH, kappa = 2._RKH) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized Upper Incomplete Gamma Function for a vector of points.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUppNR(x = 1.5_RKH, kappa = [real(RKH) :: 2, 20, 200])")
    call disp%show( getGammaIncUppNR(x = 1.5_RKH, kappa = [real(RKH) :: 2, 20, 200]) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized Upper Incomplete Gamma Function for a vector of kappa parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUppNR(x = 1._RKS, kappa = [0.1_RKS, 1._RKS, 10._RKS])")
    call disp%show( getGammaIncUppNR(x = 1._RKS, kappa = [0.1_RKS, 1._RKS, 10._RKS]) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUppNR(x = 1._RKD, kappa = [0.1_RKD, 1._RKD, 10._RKD])")
    call disp%show( getGammaIncUppNR(x = 1._RKD, kappa = [0.1_RKD, 1._RKD, 10._RKD]) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUppNR(x = 1._RKH, kappa = [0.1_RKH, 1._RKH, 10._RKH])")
    call disp%show( getGammaIncUppNR(x = 1._RKH, kappa = [0.1_RKH, 1._RKH, 10._RKH]) )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized Upper Incomplete Gamma Function for a vector of points and kappa parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()


    call disp%skip()
    call disp%show("getGammaIncUppNR(x = [1._RKS, 1._RKS, 1._RKS], kappa = [0.1_RKS, 1._RKS, 10._RKS])")
    call disp%show( getGammaIncUppNR(x = [1._RKS, 1._RKS, 1._RKS], kappa = [0.1_RKS, 1._RKS, 10._RKS]) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUppNR(x = [1._RKD, 1._RKD, 1._RKD], kappa = [0.1_RKD, 1._RKD, 10._RKD])")
    call disp%show( getGammaIncUppNR(x = [1._RKD, 1._RKD, 1._RKD], kappa = [0.1_RKD, 1._RKD, 10._RKD]) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUppNR(x = [1._RKH, 1._RKH, 1._RKH], kappa = [0.1_RKH, 1._RKH, 10._RKH])")
    call disp%show( getGammaIncUppNR(x = [1._RKH, 1._RKH, 1._RKH], kappa = [0.1_RKH, 1._RKH, 10._RKH]) )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized Upper Incomplete Gamma Function for a fixed point and kappa parameter, but with different accuracies.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUppNR(x = 1._RKH, kappa = 2._RKH, tol = [1.e-8_RKH, 1.e-16_RKH, 1.e-32_RKH])")
    call disp%show( getGammaIncUppNR(x = 1._RKH, kappa = 2._RKH, tol = [1.e-8_RKH, 1.e-16_RKH, 1.e-32_RKH]) )
    call disp%skip()


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array of the regularized Upper Incomplete Gamma function for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arraySpace, only: setLinSpace
        integer(IK) , parameter :: NP = 1000_IK
        real(RKS) :: x_RKS(NP)
        integer :: fileUnit, i

        call setLinSpace(x_RKS, 0._RKS, 10._RKS)
        open(newunit = fileUnit, file = "getGammaIncUppNR.RK.txt")
        do i = 1, NP
            write(fileUnit, "(*(g0,:,' '))") x_RKS(i), getGammaIncUppNR(x_RKS(i), kappa = [1.0_RKS, 2.5_RKS, 5.0_RKS])
        end do
        close(fileUnit)

    end block

end program example