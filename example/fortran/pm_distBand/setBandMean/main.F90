program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKC => RK128 ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distBand, only: setBandMean
    use pm_distBand, only: getBandEbreak
    use pm_physUnit, only: ERGS2KEV, KEV2ERGS

    implicit none

    integer(IK) :: info
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block

        use pm_kind, only: RKC => RK64
        real(RKC) :: mean, fluence, lbnew, ubnew, lb, ub, alpha, beta, ebreak

        call disp%skip()
        call disp%show("fluence = 2.044544e-07_RKC; lbnew = 50._RKC; ubnew = 300._RKC; lb = 50._RKC; ub = 300._RKC; alpha = -9.469590e-01_RKC; beta = -3.722981_RKC; ebreak = getBandEbreak(alpha, beta, 1.928073e+02_RKC);")
                        fluence = 2.044544e-07_RKC; lbnew = 50._RKC; ubnew = 300._RKC; lb = 50._RKC; ub = 300._RKC; alpha = -9.469590e-01_RKC; beta = -3.722981_RKC; ebreak = getBandEbreak(alpha, beta, 1.928073e+02_RKC);
        call disp%show("call setBandMean(mean, lb, ub, alpha, beta, ebreak, info)")
                        call setBandMean(mean, lb, ub, alpha, beta, ebreak, info)
        call disp%show("if (info < 0) error stop")
                        if (info < 0) error stop
        call disp%show("ERGS2KEV * fluence / mean ! photon fluence 1.084876")
        call disp%show( ERGS2KEV * fluence / mean )
        call disp%skip()

    end block

end program example