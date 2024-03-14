program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKC => RK128 ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distBand, only: getBandEbreak
    use pm_distBand, only: setBandEnergy
    use pm_physUnit, only: ERGS2KEV, KEV2ERGS

    implicit none

    integer(IK) :: info
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block

        use pm_kind, only: RKC => RK64
        real(RKC) :: energy, photon, lbnew, ubnew, lb, ub, alpha, beta, ebreak

        call disp%skip()
        call disp%show("photon = 1.084876_RKC * ERGS2KEV; lbnew = 50._RKC; ubnew = 300._RKC; lb = 50._RKC; ub = 300._RKC; alpha = -9.469590e-01_RKC; beta = -3.722981_RKC; ebreak = getBandEbreak(alpha, beta, 1.928073e+02_RKC);")
                        photon = 1.084876_RKC * ERGS2KEV; lbnew = 50._RKC; ubnew = 300._RKC; lb = 50._RKC; ub = 300._RKC; alpha = -9.469590e-01_RKC; beta = -3.722981_RKC; ebreak = getBandEbreak(alpha, beta, 1.928073e+02_RKC);
        call disp%show("call setBandEnergy(energy, photon, lb, ub, alpha, beta, ebreak, info)")
                        call setBandEnergy(energy, photon, lb, ub, alpha, beta, ebreak, info)
        call disp%show("if (info < 0) error stop")
                        if (info < 0) error stop
        call disp%show("energy ! energy fluence 2.044544e-07")
        call disp%show( energy )
        call disp%skip()

    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arraySpace, only: setLogSpace
        integer(IK), parameter :: NP = 1000_IK
        real(RKC) :: energy(4), epeak(NP), alpha(4), beta(4)
        integer :: fileUnit, i
        integer(IK) :: info(4)

        alpha = [.5_RKC, 1.5_RKC, -.5_RKC, -1.1_RKC]
        beta = -[.5_RKC, 1._RKC, 2._RKC, 3._RKC]
        call setLogSpace(epeak, log(1._RKC), log(10000._RKC))
        open(newunit = fileUnit, file = "setBandEnergy.RK.txt")
        do i = 1, NP
            call setBandEnergy(energy, 0.01_RKC, 20000._RKC, 1._RKC, 50._RKC, 300._RKC, alpha, beta, getBandEbreak(alpha, beta, epeak(i)), info)
            if (any(info < 0)) error stop
            write(fileUnit, "(*(g0,:,' '))" ) epeak(i), energy
        end do
        close(fileUnit)

    end block

end program example