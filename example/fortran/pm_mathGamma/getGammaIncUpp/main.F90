program example

    use pm_kind, only: SK
    use pm_kind, only: IK
    use pm_kind, only: LK
    use pm_kind, only: RK32, RK64, RK128
    use pm_io, only: display_type
    use pm_mathGamma, only: getGammaIncUpp

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
    call disp%show("getGammaIncUpp(x = 1.5_RK32, kappa = 2._RK32)")
    call disp%show( getGammaIncUpp(x = 1.5_RK32, kappa = 2._RK32) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUpp(x = 1.5_RK64, kappa = 2._RK64)")
    call disp%show( getGammaIncUpp(x = 1.5_RK64, kappa = 2._RK64) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUpp(x = 1.5_RK128, kappa = 2._RK128)")
    call disp%show( getGammaIncUpp(x = 1.5_RK128, kappa = 2._RK128) )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized Upper Incomplete Gamma Function for a vector of points.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()


    call disp%skip()
    call disp%show("getGammaIncUpp(x = [0._RK32, 1._RK32, 10._RK32], kappa = 2._RK32)")
    call disp%show( getGammaIncUpp(x = [0._RK32, 1._RK32, 10._RK32], kappa = 2._RK32) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUpp(x = [0._RK64, 1._RK64, 10._RK64], kappa = 2._RK64)")
    call disp%show( getGammaIncUpp(x = [0._RK64, 1._RK64, 10._RK64], kappa = 2._RK64) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUpp(x = [0._RK128, 1._RK128, 10._RK128], kappa = 2._RK128)")
    call disp%show( getGammaIncUpp(x = [0._RK128, 1._RK128, 10._RK128], kappa = 2._RK128) )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized Upper Incomplete Gamma Function for a vector of kappa parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()


    call disp%skip()
    call disp%show("getGammaIncUpp(x = 1._RK32, kappa = [0.1_RK32, 1._RK32, 10._RK32])")
    call disp%show( getGammaIncUpp(x = 1._RK32, kappa = [0.1_RK32, 1._RK32, 10._RK32]) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUpp(x = 1._RK64, kappa = [0.1_RK64, 1._RK64, 10._RK64])")
    call disp%show( getGammaIncUpp(x = 1._RK64, kappa = [0.1_RK64, 1._RK64, 10._RK64]) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUpp(x = 1._RK128, kappa = [0.1_RK128, 1._RK128, 10._RK128])")
    call disp%show( getGammaIncUpp(x = 1._RK128, kappa = [0.1_RK128, 1._RK128, 10._RK128]) )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized Upper Incomplete Gamma Function for a vector of points and kappa parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()


    call disp%skip()
    call disp%show("getGammaIncUpp(x = [1._RK32, 1._RK32, 1._RK32], kappa = [0.1_RK32, 1._RK32, 10._RK32])")
    call disp%show( getGammaIncUpp(x = [1._RK32, 1._RK32, 1._RK32], kappa = [0.1_RK32, 1._RK32, 10._RK32]) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUpp(x = [1._RK64, 1._RK64, 1._RK64], kappa = [0.1_RK64, 1._RK64, 10._RK64])")
    call disp%show( getGammaIncUpp(x = [1._RK64, 1._RK64, 1._RK64], kappa = [0.1_RK64, 1._RK64, 10._RK64]) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUpp(x = [1._RK128, 1._RK128, 1._RK128], kappa = [0.1_RK128, 1._RK128, 10._RK128])")
    call disp%show( getGammaIncUpp(x = [1._RK128, 1._RK128, 1._RK128], kappa = [0.1_RK128, 1._RK128, 10._RK128]) )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized Upper Incomplete Gamma Function for a fixed point and kappa parameter, but with different accuracies.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUpp(x = 1._RK128, kappa = 2._RK128, tol = [1.e-8_RK128, 1.e-16_RK128, 1.e-32_RK128])")
    call disp%show( getGammaIncUpp(x = 1._RK128, kappa = 2._RK128, tol = [1.e-8_RK128, 1.e-16_RK128, 1.e-32_RK128]) )
    call disp%skip()


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array of the regularized Upper Incomplete Gamma function for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arraySpace, only: setLinSpace
        integer(IK) , parameter :: NP = 1000_IK
        real(RK32) :: x_RK32(NP)
        integer :: fileUnit, i

        call setLinSpace(x_RK32, 0._RK32, 10._RK32)
        open(newunit = fileUnit, file = "getGammaIncUpp.RK.txt")
        do i = 1, NP
            write(fileUnit, "(*(g0,:,' '))" ) x_RK32(i), getGammaIncUpp(x_RK32(i), kappa = [1.0_RK32, 2.5_RK32, 5.0_RK32])
        end do
        close(fileUnit)

    end block

end program example