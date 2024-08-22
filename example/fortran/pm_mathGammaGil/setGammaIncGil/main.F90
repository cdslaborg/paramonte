program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKS, RKD, RKH
    use pm_mathGammaGil, only: setGammaIncGil
    use pm_io, only: display_type

    implicit none

    integer(IK) , parameter :: NP = 1000_IK
    real(RKH)               :: gamIncLow_RKH, gamIncUpp_RKH, x_RKH, kappa_RKH
    real(RKD)               :: gamIncLow_RKD, gamIncUpp_RKD, x_RKD, kappa_RKD
    real(RKS)               :: gamIncLow_RKS, gamIncUpp_RKS, x_RKS, kappa_RKS
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
    call disp%show("! Compute the regularized Lower Incomplete Gamma Function using its series representation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("x_RKS")
    call disp%show( x_RKS )
    call disp%show("kappa_RKS")
    call disp%show( kappa_RKS )
    call disp%show("call setGammaIncGil(gamIncLow_RKS, gamIncUpp_RKS, x_RKS, kappa = kappa_RKS, info = info)")
                    call setGammaIncGil(gamIncLow_RKS, gamIncUpp_RKS, x_RKS, kappa = kappa_RKS, info = info)
    call disp%show("[gamIncLow_RKS, gamIncUpp_RKS]")
    call disp%show( [gamIncLow_RKS, gamIncUpp_RKS] )
    call disp%show("info")
    call disp%show( info )
    call disp%skip()

    call disp%skip()
    call disp%show("x_RKD")
    call disp%show( x_RKD )
    call disp%show("kappa_RKD")
    call disp%show( kappa_RKD )
    call disp%show("call setGammaIncGil(gamIncLow_RKD, gamIncUpp_RKD, x_RKD, kappa = kappa_RKD, info = info)")
                    call setGammaIncGil(gamIncLow_RKD, gamIncUpp_RKD, x_RKD, kappa = kappa_RKD, info = info)
    call disp%show("[gamIncLow_RKD, gamIncUpp_RKD]")
    call disp%show( [gamIncLow_RKD, gamIncUpp_RKD] )
    call disp%show("info")
    call disp%show( info )
    call disp%skip()

    call disp%skip()
    call disp%show("x_RKH")
    call disp%show( x_RKH )
    call disp%show("kappa_RKH")
    call disp%show( kappa_RKH )
    call disp%show("call setGammaIncGil(gamIncLow_RKH, gamIncUpp_RKH, x_RKH, kappa = kappa_RKH, info = info)")
                    call setGammaIncGil(gamIncLow_RKH, gamIncUpp_RKH, x_RKH, kappa = kappa_RKH, info = info)
    call disp%show("[gamIncLow_RKH, gamIncUpp_RKH]")
    call disp%show( [gamIncLow_RKH, gamIncUpp_RKH] )
    call disp%show("info")
    call disp%show( info )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized Lower Incomplete Gamma Function for a vector of points and shape parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: RKG => RKD
        integer(IK) :: i, rprecision
        integer(IK), allocatable :: exprange(:), info(:)
        real(RKG), allocatable :: gamIncLow(:), gamIncUpp(:)
        call disp%skip()
        call disp%show("rprecision = precision(0._RKG) * 2")
                        rprecision = precision(0._RKG) * 2
        call disp%show("rprecision")
        call disp%show( rprecision )
        call disp%show("exprange = [(i, i = -rprecision, rprecision)]")
                        exprange = [(i, i = -rprecision, rprecision)]
        call disp%show("exprange")
        call disp%show( exprange )
        call disp%show("allocate(gamIncLow(size(exprange)), gamIncUpp(size(exprange)), info(size(exprange)))")
                        allocate(gamIncLow(size(exprange)), gamIncUpp(size(exprange)), info(size(exprange)))
        call disp%show("call setGammaIncGil(gamIncLow, gamIncUpp, x = 10._RKG**exprange, kappa = 10._RKG**exprange, info = info)")
                        call setGammaIncGil(gamIncLow, gamIncUpp, x = 10._RKG**exprange, kappa = 10._RKG**exprange, info = info)
        call disp%show("reshape([10._RKG**exprange, gamIncLow, real(info, RKG)], shape = [size(info), 3])")
        call disp%show( reshape([10._RKG**exprange, gamIncLow, real(info, RKG)], shape = [size(info), 3]) )
        call disp%skip()
    end block

    block
        use pm_kind, only: RKG => RKH
        integer(IK) :: i, rprecision
        integer(IK), allocatable :: exprange(:), info(:)
        real(RKG), allocatable :: gamIncLow(:), gamIncUpp(:)
        call disp%skip()
        call disp%show("rprecision = precision(0._RKG) * 2")
                        rprecision = precision(0._RKG) * 2
        call disp%show("rprecision")
        call disp%show( rprecision )
        call disp%show("exprange = [(i, i = -rprecision, rprecision)]")
                        exprange = [(i, i = -rprecision, rprecision)]
        call disp%show("exprange")
        call disp%show( exprange )
        call disp%show("allocate(gamIncLow(size(exprange)), gamIncUpp(size(exprange)), info(size(exprange)))")
                        allocate(gamIncLow(size(exprange)), gamIncUpp(size(exprange)), info(size(exprange)))
        call disp%show("call setGammaIncGil(gamIncLow, gamIncUpp, x = 10._RKG**exprange, kappa = 10._RKG**exprange, info = info)")
                        call setGammaIncGil(gamIncLow, gamIncUpp, x = 10._RKG**exprange, kappa = 10._RKG**exprange, info = info)
        call disp%show("reshape([10._RKG**exprange, gamIncLow, real(info, RKG)], shape = [size(info), 3])")
        call disp%show( reshape([10._RKG**exprange, gamIncLow, real(info, RKG)], shape = [size(info), 3]) )
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
        open(newunit = fileUnit, file = "setGammaIncGil.RK.txt")
        do i = 1, NP
            call setGammaIncGil(gamIncLow_RKS, gamIncUpp_RKS, x_RKS(i), kappa = kappa_RKS, info = info)
            write(fileUnit,"(2(g0,:,' '))") x_RKS(i), gamIncLow_RKS
        end do
        close(fileUnit)

    end block

end program example