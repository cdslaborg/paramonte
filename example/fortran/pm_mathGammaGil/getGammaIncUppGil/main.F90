program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKS, RKD, RKH
    use pm_mathGammaGil, only: getGammaIncUppGil
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized Lower Incomplete Gamma Function using its series representation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUppGil(x = 1.5_RKS, kappa = 2._RKS)")
    call disp%show( getGammaIncUppGil(x = 1.5_RKS, kappa = 2._RKS) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUppGil(x = 1.5_RKD, kappa = 2._RKD)")
    call disp%show( getGammaIncUppGil(x = 1.5_RKD, kappa = 2._RKD) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUppGil(x = 1.5_RKH, kappa = 2._RKH)")
    call disp%show( getGammaIncUppGil(x = 1.5_RKH, kappa = 2._RKH) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized Lower Incomplete Gamma Function for a vector of points.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUppGil(x = [0._RKS, 1._RKS, 10._RKS], kappa = 2._RKS)")
    call disp%show( getGammaIncUppGil(x = [0._RKS, 1._RKS, 10._RKS], kappa = 2._RKS) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUppGil(x = [0._RKD, 1._RKD, 10._RKD], kappa = 2._RKD)")
    call disp%show( getGammaIncUppGil(x = [0._RKD, 1._RKD, 10._RKD], kappa = 2._RKD) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUppGil(x = [0._RKH, 1._RKH, 10._RKH], kappa = 2._RKH)")
    call disp%show( getGammaIncUppGil(x = [0._RKH, 1._RKH, 10._RKH], kappa = 2._RKH) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized Lower Incomplete Gamma Function for a vector of shape parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUppGil(x = 1._RKS, kappa = [0.1_RKS, 1._RKS, 10._RKS])")
    call disp%show( getGammaIncUppGil(x = 1._RKS, kappa = [0.1_RKS, 1._RKS, 10._RKS]) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUppGil(x = 1._RKD, kappa = [0.1_RKD, 1._RKD, 10._RKD])")
    call disp%show( getGammaIncUppGil(x = 1._RKD, kappa = [0.1_RKD, 1._RKD, 10._RKD]) )
    call disp%skip()

    call disp%skip()
    call disp%show("getGammaIncUppGil(x = 1._RKH, kappa = [0.1_RKH, 1._RKH, 10._RKH])")
    call disp%show( getGammaIncUppGil(x = 1._RKH, kappa = [0.1_RKH, 1._RKH, 10._RKH]) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized Lower Incomplete Gamma Function for a vector of points and shape parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: RKG => RKS
        integer(IK) :: i, rprecision
        integer(IK), allocatable :: exprange(:), info(:)
        real(RKG), allocatable :: gamIncLow(:)
        call disp%skip()
        call disp%show("rprecision = precision(0._RKG) * 2")
                        rprecision = precision(0._RKG) * 2
        call disp%show("rprecision")
        call disp%show( rprecision )
        call disp%show("exprange = [(i, i = -rprecision, rprecision)]")
                        exprange = [(i, i = -rprecision, rprecision)]
        call disp%show("exprange")
        call disp%show( exprange )
        call disp%show("allocate(gamIncLow(size(exprange)), info(size(exprange)))")
                        allocate(gamIncLow(size(exprange)), info(size(exprange)))
        call disp%show("gamIncLow = getGammaIncUppGil(x = 10._RKG**exprange, kappa = 10._RKG**exprange)")
                        gamIncLow = getGammaIncUppGil(x = 10._RKG**exprange, kappa = 10._RKG**exprange)
        call disp%show("reshape([10._RKG**exprange, gamIncLow], shape = [size(info), 2])")
        call disp%show( reshape([10._RKG**exprange, gamIncLow], shape = [size(info), 2]) )
        call disp%skip()
    end block

    block
        use pm_kind, only: RKG => RKD
        integer(IK) :: i, rprecision
        integer(IK), allocatable :: exprange(:), info(:)
        real(RKG), allocatable :: gamIncLow(:)
        call disp%skip()
        call disp%show("rprecision = precision(0._RKG) * 2")
                        rprecision = precision(0._RKG) * 2
        call disp%show("rprecision")
        call disp%show( rprecision )
        call disp%show("exprange = [(i, i = -rprecision, rprecision)]")
                        exprange = [(i, i = -rprecision, rprecision)]
        call disp%show("exprange")
        call disp%show( exprange )
        call disp%show("allocate(gamIncLow(size(exprange)), info(size(exprange)))")
                        allocate(gamIncLow(size(exprange)), info(size(exprange)))
        call disp%show("gamIncLow = getGammaIncUppGil(x = 10._RKG**exprange, kappa = 10._RKG**exprange)")
                        gamIncLow = getGammaIncUppGil(x = 10._RKG**exprange, kappa = 10._RKG**exprange)
        call disp%show("reshape([10._RKG**exprange, gamIncLow], shape = [size(info), 2])")
        call disp%show( reshape([10._RKG**exprange, gamIncLow], shape = [size(info), 2]) )
        call disp%skip()
    end block

    block
        use pm_kind, only: RKG => RKH
        integer(IK) :: i, rprecision
        integer(IK), allocatable :: exprange(:), info(:)
        real(RKG), allocatable :: gamIncLow(:)
        call disp%skip()
        call disp%show("rprecision = precision(0._RKG) * 2")
                        rprecision = precision(0._RKG) * 2
        call disp%show("rprecision")
        call disp%show( rprecision )
        call disp%show("exprange = [(i, i = -rprecision, rprecision)]")
                        exprange = [(i, i = -rprecision, rprecision)]
        call disp%show("exprange")
        call disp%show( exprange )
        call disp%show("allocate(gamIncLow(size(exprange)), info(size(exprange)))")
                        allocate(gamIncLow(size(exprange)), info(size(exprange)))
        call disp%show("gamIncLow = getGammaIncUppGil(x = 10._RKG**exprange, kappa = 10._RKG**exprange)")
                        gamIncLow = getGammaIncUppGil(x = 10._RKG**exprange, kappa = 10._RKG**exprange)
        call disp%show("reshape([10._RKG**exprange, gamIncLow], shape = [size(info), 2])")
        call disp%show( reshape([10._RKG**exprange, gamIncLow], shape = [size(info), 2]) )
        call disp%skip()
    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array of the regularized Lower Incomplete Gamma function for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arraySpace, only: setLinSpace
        integer(IK) , parameter :: NP = 1000_IK
        real(RKS) :: x_RKS(NP)
        integer :: fileUnit, i

        call setLinSpace(x_RKS, 0._RKS, 10._RKS)
        open(newunit = fileUnit, file = "getGammaIncUppGil.RK.txt")
        do i = 1, NP
            write(fileUnit, "(*(g0,:,' '))") x_RKS(i), getGammaIncUppGil(x_RKS(i), kappa = [1.0_RKS, 2.5_RKS, 5.0_RKS])
        end do
        close(fileUnit)

    end block

end program example