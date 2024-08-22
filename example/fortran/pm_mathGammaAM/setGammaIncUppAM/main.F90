program example

    use pm_kind, only: SK
    use pm_kind, only: IK
    use pm_kind, only: LK
    use pm_kind, only: RKS, RKD, RKH
    use pm_io, only: display_type
    use pm_mathGammaAM, only: setGammaIncUppAM

    implicit none

    integer(IK) , parameter :: NP = 1000_IK
    real(RKH)               :: gamIncUpp_RKH, x_RKH, kappa_RKH
    real(RKD)               :: gamIncUpp_RKD, x_RKD, kappa_RKD
    real(RKS)               :: gamIncUpp_RKS, x_RKS, kappa_RKS
    integer(IK)             :: info

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    kappa_RKH = 1.5_RKH
    kappa_RKD = 1.5_RKD
    kappa_RKS = 1.5_RKS

    x_RKH = 2._RKH
    x_RKD = 2._RKD
    x_RKS = 2._RKS

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized Upper Incomplete Gamma Function using its series representation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("x_RKS")
    call disp%show( x_RKS )
    call disp%show("kappa_RKS")
    call disp%show( kappa_RKS )
    call disp%show("call setGammaIncUppAM(gamIncUpp_RKS, x_RKS, log_gamma(kappa_RKS), kappa = kappa_RKS, info = info)")
                    call setGammaIncUppAM(gamIncUpp_RKS, x_RKS, log_gamma(kappa_RKS), kappa = kappa_RKS, info = info)
    call disp%show("gamIncUpp_RKS")
    call disp%show( gamIncUpp_RKS )
    call disp%show("info")
    call disp%show( info )
    call disp%skip()

    call disp%skip()
    call disp%show("x_RKD")
    call disp%show( x_RKD)
    call disp%show("kappa_RKD")
    call disp%show( kappa_RKD)
    call disp%show("call setGammaIncUppAM(gamIncUpp_RKD, x_RKD, log_gamma(kappa_RKD), kappa = kappa_RKD, info = info)")
                    call setGammaIncUppAM(gamIncUpp_RKD, x_RKD, log_gamma(kappa_RKD), kappa = kappa_RKD, info = info)
    call disp%show("gamIncUpp_RKD")
    call disp%show( gamIncUpp_RKD)
    call disp%show("info")
    call disp%show( info )
    call disp%skip()

    call disp%skip()
    call disp%show("x_RKH")
    call disp%show( x_RKH )
    call disp%show("kappa_RKH")
    call disp%show( kappa_RKH )
    call disp%show("call setGammaIncUppAM(gamIncUpp_RKH, x_RKH, log_gamma(kappa_RKH), kappa = kappa_RKH, info = info)")
                    call setGammaIncUppAM(gamIncUpp_RKH, x_RKH, log_gamma(kappa_RKH), kappa = kappa_RKH, info = info)
    call disp%show("gamIncUpp_RKH")
    call disp%show( gamIncUpp_RKH )
    call disp%show("info")
    call disp%show( info )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array of the regularized Upper Incomplete Gamma function for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arraySpace, only: setLinSpace
        real(RKS) :: x_RKS(NP)
        integer :: fileUnit, i

        call setLinSpace(x_RKS, 0._RKS, 8._RKS)
        open(newunit = fileUnit, file = "setGammaIncUppAM.RK.txt")
        do i = 1, NP
            call setGammaIncUppAM(gamIncUpp_RKS, x_RKS(i), log_gamma(kappa_RKS), kappa = kappa_RKS, info = info)
            write(fileUnit,"(2(g0,:,' '))") x_RKS(i), gamIncUpp_RKS
        end do
        close(fileUnit)

    end block

end program example