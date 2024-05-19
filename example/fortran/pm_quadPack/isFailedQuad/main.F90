program example

    use pm_kind, only: SK, IK, RKH, RKG => RKH ! All other real kinds are also supported.
    use pm_quadTest, only: int1_type
    use pm_quadTest, only: int2_type
    use pm_quadTest, only: int3_type
    use pm_quadTest, only: int4_type
    use pm_quadTest, only: int5_type
    use pm_quadTest, only: int6_type
    use pm_quadTest, only: int7_type
    use pm_quadTest, only: int8_type
    use pm_quadTest, only: int9_type
    use pm_quadTest, only: intGamUpp_type
    use pm_quadTest, only: intSinCos_type
    use pm_quadTest, only: intNormPDF_type
    use pm_quadTest, only: intDoncker1_type
    use pm_quadTest, only: intDoncker2_type
    use pm_quadTest, only: intLogNormPDF_type
    use pm_quadTest, only: intPentaGammaInf_type
    use pm_quadTest, only: intGenExpGammaPDF_type
    use pm_quadTest, only: test_isFailedQuad
    use pm_quadTest, only: intCauchy1_type
    use pm_quadTest, only: intCauchy2_type
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compare integration via the Adaptive Global Gauss-Kronrod Quadrature method with (QAG) and without (QAGS) extrapolation acceleration.")
    call disp%show("! Note the significant efficiency improvement in the integration via QAGS when the integrand has integrable singularity at some points.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call test_isFailedQuad(disp, int1_type(-4._RKH, +4._RKH))
    call test_isFailedQuad(disp, int2_type(+2._RKH, +3._RKH))
    call test_isFailedQuad(disp, int3_type(ub = +3._RKH))
    call test_isFailedQuad(disp, int4_type())
    call test_isFailedQuad(disp, intSinCos_type(-4_IK, +4_IK, a = 10._RKH, b = -10._RKH))
    call test_isFailedQuad(disp, intNormPDF_type(-3._RKH, +3._RKH))
    call test_isFailedQuad(disp, intLogNormPDF_type(lb = exp(-6._RKG), ub = exp(+6._RKG)))

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compare integration via:")
    call disp%show("!       the Adaptive Global Gauss-Kronrod Quadrature method with automatic break point handling and extrapolation acceleration (QAGS),")
    call disp%show("! with,")
    call disp%show("!       the Adaptive Global Gauss-Kronrod Quadrature method with specified break point handling and extrapolation acceleration (QAGP).")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call test_isFailedQuad(disp, int5_type(0._RKG, 3._RKG))
    call test_isFailedQuad(disp, int5_type(-2._RKG, +5._RKG))

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compare integration via:")
    call disp%show("!       the Adaptive Global Gauss-Kronrod Quadrature method without extrapolation on a semi-infinite range `(lb, +Inf)`,")
    call disp%show("! with,")
    call disp%show("!       the Adaptive Global Gauss-Kronrod Quadrature method with extrapolation on a semi-infinite range `(lb, +Inf)`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call test_isFailedQuad(disp, intGamUpp_type())
    !call disp%show("call test_isFailedQuad(disp, intGamUpp_type(lb = 1._RKG, alpha = +1.0_RKG, beta = 3._RKG))")
    !                call test_isFailedQuad(disp, intGamUpp_type(lb = 1._RKG, alpha = +1.0_RKG, beta = 3._RKG))
    !call disp%show("call test_isFailedQuad(disp, intGamUpp_type(lb = 1._RKG, alpha = +0.5_RKG, beta = 3._RKG))")
    !                call test_isFailedQuad(disp, intGamUpp_type(lb = 1._RKG, alpha = +0.5_RKG, beta = 3._RKG))
    !call disp%show("call test_isFailedQuad(disp, intGamUpp_type(lb = 1._RKG, alpha = +0.0_RKG, beta = 3._RKG))")
    !                call test_isFailedQuad(disp, intGamUpp_type(lb = 1._RKG, alpha = +0.0_RKG, beta = 3._RKG))
    !call disp%show("call test_isFailedQuad(disp, intGamUpp_type(lb = 1._RKG, alpha = -0.5_RKG, beta = 3._RKG))")
    !                call test_isFailedQuad(disp, intGamUpp_type(lb = 1._RKG, alpha = -0.5_RKG, beta = 3._RKG))
    !call disp%show("call test_isFailedQuad(disp, intGamUpp_type(lb = 1._RKG, alpha = -1.0_RKG, beta = 3._RKG))")
    !                call test_isFailedQuad(disp, intGamUpp_type(lb = 1._RKG, alpha = -1.0_RKG, beta = 3._RKG))
    !call disp%show("call test_isFailedQuad(disp, intGamUpp_type(lb = 1._RKG, alpha = -1.5_RKG, beta = 3._RKG))")
    !                call test_isFailedQuad(disp, intGamUpp_type(lb = 1._RKG, alpha = -1.5_RKG, beta = 3._RKG))
    !call disp%show("call test_isFailedQuad(disp, intGamUpp_type(lb = 1._RKG, alpha = -2.0_RKG, beta = 3._RKG))")
    !                call test_isFailedQuad(disp, intGamUpp_type(lb = 1._RKG, alpha = -2.0_RKG, beta = 3._RKG))
    call test_isFailedQuad(disp, intLogNormPDF_type())
    call test_isFailedQuad(disp, intDoncker1_type())
    call test_isFailedQuad(disp, int6_type())
    call test_isFailedQuad(disp, int7_type())

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compare integration via:")
    call disp%show("!       the Adaptive Global Gauss-Kronrod Quadrature method without extrapolation on a semi-infinite range `(-Inf, ub)`,")
    call disp%show("! with,")
    call disp%show("!       the Adaptive Global Gauss-Kronrod Quadrature method with extrapolation on a semi-infinite range `(-Inf, ub)`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call test_isFailedQuad(disp, int8_type())
    call test_isFailedQuad(disp, intDoncker2_type())
    call test_isFailedQuad(disp, intDoncker2_type(ub = -1._RKG))

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compare integration via:")
    call disp%show("!       the Adaptive Global Gauss-Kronrod Quadrature method without extrapolation on a full-infinite range `(-Inf, +Inf)`,")
    call disp%show("! with,")
    call disp%show("!       the Adaptive Global Gauss-Kronrod Quadrature method with extrapolation on a full-infinite range `(-Inf, +Inf)`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call test_isFailedQuad(disp, int9_type())
    call test_isFailedQuad(disp, intNormPDF_type())
    call test_isFailedQuad(disp, intGenExpGammaPDF_type())
    call test_isFailedQuad(disp, intGenExpGammaPDF_type(kappa = 5._RKH, invOmega = 0.1_RKH, logSigma = 2._RKH))
    call test_isFailedQuad(disp, intPentaGammaInf_type())

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cauchy Principal Value over a finite interval.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call test_isFailedQuad(disp, intCauchy1_type())
    call test_isFailedQuad(disp, intCauchy2_type())

end program example