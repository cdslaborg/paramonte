program example

    use pm_kind, only: SK
    use pm_kind, only: IK
    use pm_kind, only: LK
    use pm_kind, only: RK32, RK64, RK128
    use pm_io, only: display_type
    use pm_mathGamma, only: setGammaIncLow

    implicit none

    integer(IK) , parameter :: NP = 1000_IK
    real(RK128)             :: gamIncLow_RK128, x_RK128, kappa_RK128
    real(RK64 )             :: gamIncLow_RK64 , x_RK64 , kappa_RK64
    real(RK32 )             :: gamIncLow_RK32 , x_RK32 , kappa_RK32
    integer(IK)             :: info

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    kappa_RK128 = 1.5_RK128
    kappa_RK64  = 1.5_RK64
    kappa_RK32  = 1.5_RK32

    x_RK128 = 2._RK128
    x_RK64  = 2._RK64
    x_RK32  = 2._RK32

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized Lower Incomplete Gamma Function using its series representation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("x_RK32")
    call disp%show( x_RK32 )
    call disp%show("kappa_RK32")
    call disp%show( kappa_RK32 )
    call disp%show("call setGammaIncLow(gamIncLow_RK32, x_RK32, logGammaKappa = log_gamma(kappa_RK32), kappa = kappa_RK32, info = info)")
                    call setGammaIncLow(gamIncLow_RK32, x_RK32, logGammaKappa = log_gamma(kappa_RK32), kappa = kappa_RK32, info = info)
    call disp%show("gamIncLow_RK32")
    call disp%show( gamIncLow_RK32 )
    call disp%show("info")
    call disp%show( info )
    call disp%skip()

    call disp%skip()
    call disp%show("x_RK64")
    call disp%show( x_RK64 )
    call disp%show("kappa_RK64")
    call disp%show( kappa_RK64 )
    call disp%show("call setGammaIncLow(gamIncLow_RK64, x_RK64, logGammaKappa = log_gamma(kappa_RK64), kappa = kappa_RK64, info = info)")
                    call setGammaIncLow(gamIncLow_RK64, x_RK64, logGammaKappa = log_gamma(kappa_RK64), kappa = kappa_RK64, info = info)
    call disp%show("gamIncLow_RK64")
    call disp%show( gamIncLow_RK64 )
    call disp%show("info")
    call disp%show( info )
    call disp%skip()

    call disp%skip()
    call disp%show("x_RK128")
    call disp%show( x_RK128 )
    call disp%show("kappa_RK128")
    call disp%show( kappa_RK128 )
    call disp%show("call setGammaIncLow(gamIncLow_RK128, x_RK128, logGammaKappa = log_gamma(kappa_RK128), kappa = kappa_RK128, info = info)")
                    call setGammaIncLow(gamIncLow_RK128, x_RK128, logGammaKappa = log_gamma(kappa_RK128), kappa = kappa_RK128, info = info)
    call disp%show("gamIncLow_RK128")
    call disp%show( gamIncLow_RK128 )
    call disp%show("info")
    call disp%show( info )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array of the regularized Lower Incomplete Gamma function for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arraySpace, only: setLinSpace
        real(RK32) :: x_RK32(NP)
        integer :: fileUnit, i

        call setLinSpace(x_RK32, 0._RK32, 8._RK32)
        open(newunit = fileUnit, file = "setGammaIncLow.RK.txt")
        do i = 1, NP
            call setGammaIncLow(gamIncLow_RK32, x_RK32(i), logGammaKappa = log_gamma(kappa_RK32), kappa = kappa_RK32, info = info)
            write(fileUnit,"(2(g0,:,' '))") x_RK32(i), gamIncLow_RK32
        end do
        close(fileUnit)

    end block

end program example