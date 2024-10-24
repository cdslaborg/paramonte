program example

    use pm_kind, only: SK
    use pm_kind, only: IK
    use pm_kind, only: LK
    use pm_kind, only: RKS, RKD, RKH
    use pm_io, only: display_type
    use pm_mathGammaNR, only: setGammaIncLowSeriesNR

    implicit none

    integer(IK) , parameter :: NP = 1000_IK
    real(RKH)               :: gamIncLow_RKH, x_RKH, kappa_RKH
    real(RKD)               :: gamIncLow_RKD, x_RKD, kappa_RKD
    real(RKS)               :: gamIncLow_RKS, x_RKS, kappa_RKS
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
    call disp%show("x_RKS")
    call disp%show( x_RKS )
    call disp%show("kappa_RKS")
    call disp%show( kappa_RKS )
    call disp%show("call setGammaIncLowSeriesNR(gamIncLow_RKS, x_RKS, logGammaKappa = log_gamma(kappa_RKS), kappa = kappa_RKS, info = info)")
                    call setGammaIncLowSeriesNR(gamIncLow_RKS, x_RKS, logGammaKappa = log_gamma(kappa_RKS), kappa = kappa_RKS, info = info)
    call disp%show("gamIncLow_RKS")
    call disp%show( gamIncLow_RKS )
    call disp%show("info")
    call disp%show( info )
    call disp%skip()

    call disp%skip()
    call disp%show("x_RKD")
    call disp%show( x_RKD)
    call disp%show("kappa_RKD")
    call disp%show( kappa_RKD)
    call disp%show("call setGammaIncLowSeriesNR(gamIncLow_RKD, x_RKD, logGammaKappa = log_gamma(kappa_RKD), kappa = kappa_RKD, info = info)")
                    call setGammaIncLowSeriesNR(gamIncLow_RKD, x_RKD, logGammaKappa = log_gamma(kappa_RKD), kappa = kappa_RKD, info = info)
    call disp%show("gamIncLow_RKD")
    call disp%show( gamIncLow_RKD)
    call disp%show("info")
    call disp%show( info )
    call disp%skip()

    call disp%skip()
    call disp%show("x_RKH")
    call disp%show( x_RKH )
    call disp%show("kappa_RKH")
    call disp%show( kappa_RKH )
    call disp%show("call setGammaIncLowSeriesNR(gamIncLow_RKH, x_RKH, logGammaKappa = log_gamma(kappa_RKH), kappa = kappa_RKH, info = info)")
                    call setGammaIncLowSeriesNR(gamIncLow_RKH, x_RKH, logGammaKappa = log_gamma(kappa_RKH), kappa = kappa_RKH, info = info)
    call disp%show("gamIncLow_RKH")
    call disp%show( gamIncLow_RKH )
    call disp%show("info")
    call disp%show( info )
    call disp%skip()

    block
        use pm_kind, only: RKG => RKS
        integer(IK) :: i, rprecision
        integer(IK), allocatable :: exprange(:), info(:)
        real(RKG), allocatable :: gamIncLow(:)
        call disp%skip()
        call disp%show("rprecision = precision(0._RKG) / 2")
                        rprecision = precision(0._RKG) / 2
        call disp%show("rprecision")
        call disp%show( rprecision )
        call disp%show("exprange = [(i, i = -rprecision, rprecision)]")
                        exprange = [(i, i = -rprecision, rprecision)]
        call disp%show("exprange")
        call disp%show( exprange )
        call disp%show("allocate(gamIncLow(size(exprange)), info(size(exprange)))")
                        allocate(gamIncLow(size(exprange)), info(size(exprange)))
        call disp%show("call setGammaIncLowSeriesNR(gamIncLow, 10._RKG**exprange, logGammaKappa = log_gamma(10._RKG**exprange), kappa = 10._RKG**exprange, info = info)")
                        call setGammaIncLowSeriesNR(gamIncLow, 10._RKG**exprange, logGammaKappa = log_gamma(10._RKG**exprange), kappa = 10._RKG**exprange, info = info)
        call disp%show("gamIncLow")
        call disp%show( gamIncLow )
        call disp%show("info")
        call disp%show( info )
        call disp%skip()
    end block

    block
        use pm_kind, only: RKG => RKD
        integer(IK) :: i, rprecision
        integer(IK), allocatable :: exprange(:), info(:)
        real(RKG), allocatable :: gamIncLow(:)
        call disp%skip()
        call disp%show("rprecision = precision(0._RKG) / 2")
                        rprecision = precision(0._RKG) / 2
        call disp%show("rprecision")
        call disp%show( rprecision )
        call disp%show("exprange = [(i, i = -rprecision, rprecision)]")
                        exprange = [(i, i = -rprecision, rprecision)]
        call disp%show("exprange")
        call disp%show( exprange )
        call disp%show("allocate(gamIncLow(size(exprange)), info(size(exprange)))")
                        allocate(gamIncLow(size(exprange)), info(size(exprange)))
        call disp%show("call setGammaIncLowSeriesNR(gamIncLow, 10._RKG**exprange, logGammaKappa = log_gamma(10._RKG**exprange), kappa = 10._RKG**exprange, info = info)")
                        call setGammaIncLowSeriesNR(gamIncLow, 10._RKG**exprange, logGammaKappa = log_gamma(10._RKG**exprange), kappa = 10._RKG**exprange, info = info)
        call disp%show("gamIncLow")
        call disp%show( gamIncLow )
        call disp%show("info")
        call disp%show( info )
        call disp%skip()
    end block

    block
        use pm_kind, only: RKG => RKH
        integer(IK) :: i, rprecision
        integer(IK), allocatable :: exprange(:), info(:)
        real(RKG), allocatable :: gamIncLow(:)
        call disp%skip()
        call disp%show("rprecision = precision(0._RKG) / 2")
                        rprecision = precision(0._RKG) / 2
        call disp%show("rprecision")
        call disp%show( rprecision )
        call disp%show("exprange = [(i, i = -rprecision, rprecision)]")
                        exprange = [(i, i = -rprecision, rprecision)]
        call disp%show("exprange")
        call disp%show( exprange )
        call disp%show("allocate(gamIncLow(size(exprange)), info(size(exprange)))")
                        allocate(gamIncLow(size(exprange)), info(size(exprange)))
        call disp%show("call setGammaIncLowSeriesNR(gamIncLow, 10._RKG**exprange, logGammaKappa = log_gamma(10._RKG**exprange), kappa = 10._RKG**exprange, info = info)")
                        call setGammaIncLowSeriesNR(gamIncLow, 10._RKG**exprange, logGammaKappa = log_gamma(10._RKG**exprange), kappa = 10._RKG**exprange, info = info)
        call disp%show("gamIncLow")
        call disp%show( gamIncLow )
        call disp%show("info")
        call disp%show( info )
        call disp%skip()
    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array of the regularized Lower Incomplete Gamma function for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arraySpace, only: setLinSpace
        real(RKS) :: x_RKS(NP)
        integer :: fileUnit, i

        call setLinSpace(x_RKS, 0._RKS, 8._RKS)
        open(newunit = fileUnit, file = "setGammaIncLowSeriesNR.RK.txt")
        do i = 1, NP
            call setGammaIncLowSeriesNR(gamIncLow_RKS, x_RKS(i), logGammaKappa = log_gamma(kappa_RKS), kappa = kappa_RKS, info = info)
            write(fileUnit,"(2(g0,:,' '))") x_RKS(i), gamIncLow_RKS
        end do
        close(fileUnit)

    end block

end program example