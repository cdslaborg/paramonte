program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKG => RKH ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distBand, only: getBandEbreak
    use pm_distBand, only: setBandEnergy
    use pm_physUnit, only: ERGS2KEV, KEV2ERGS

    implicit none

    integer(IK) :: info
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block

        use pm_kind, only: RKG => RKD
        real(RKG) :: energy, photon, lbnew, ubnew, lb, ub, alpha, beta, ebreak

        call disp%skip()
        call disp%show("photon = 1.084876_RKG * ERGS2KEV; lbnew = 50._RKG; ubnew = 300._RKG; lb = 50._RKG; ub = 300._RKG; alpha = -9.469590e-01_RKG; beta = -3.722981_RKG; ebreak = getBandEbreak(alpha, beta, 1.928073e+02_RKG);")
                        photon = 1.084876_RKG * ERGS2KEV; lbnew = 50._RKG; ubnew = 300._RKG; lb = 50._RKG; ub = 300._RKG; alpha = -9.469590e-01_RKG; beta = -3.722981_RKG; ebreak = getBandEbreak(alpha, beta, 1.928073e+02_RKG);
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
        real(RKG) :: energy(4), epeak(NP), alpha(4), beta(4)
        integer :: fileUnit, i
        integer(IK) :: info(4)

        alpha = [.5_RKG, 1.5_RKG, -.5_RKG, -1.1_RKG]
        beta = -[.5_RKG, 1._RKG, 2._RKG, 3._RKG]
        call setLogSpace(epeak, log(1._RKG), log(10000._RKG))
        open(newunit = fileUnit, file = "setBandEnergy.RK.txt")
        do i = 1, NP
            call setBandEnergy(energy, 0.01_RKG, 20000._RKG, 1._RKG, 50._RKG, 300._RKG, alpha, beta, getBandEbreak(alpha, beta, epeak(i)), info)
            if (any(info < 0)) error stop
            write(fileUnit, "(*(g0,:,' '))") epeak(i), energy
        end do
        close(fileUnit)

    end block

end program example