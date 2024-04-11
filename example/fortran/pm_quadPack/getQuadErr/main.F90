program example

    use pm_kind, only: SK, IK, RKH, RKC => RKH ! All other real kinds are also supported.
    use pm_except, only: getInfNeg, getInfPos
    use pm_quadTest, only: int1_type
    use pm_quadTest, only: int2_type
    use pm_quadTest, only: int3_type
    use pm_quadTest, only: int4_type
    use pm_quadTest, only: int5_type
    use pm_quadTest, only: int6_type
    use pm_quadTest, only: int7_type
    use pm_quadTest, only: int8_type
    use pm_quadTest, only: int9_type
    use pm_quadTest, only: intSinCos_type
    use pm_quadTest, only: intNormPDF_type
    use pm_quadTest, only: intDoncker1_type
    use pm_quadTest, only: intDoncker2_type
    use pm_quadTest, only: intLogNormPDF_type
    use pm_quadTest, only: intPentaGammaInf_type
    use pm_quadTest, only: intCauchy1_type
    use pm_quadTest, only: intCauchy2_type
    use pm_quadTest, only: test_getQuadErr
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

    call test_getQuadErr(disp, int1_type(-4._RKH, +4._RKH))
    call test_getQuadErr(disp, int2_type(+2._RKH, +3._RKH))
    call test_getQuadErr(disp, int3_type(ub = +3._RKH))
    call test_getQuadErr(disp, int4_type())
    call test_getQuadErr(disp, intSinCos_type(-4_IK, +4_IK, a = 10._RKH, b = -10._RKH))
    call test_getQuadErr(disp, intNormPDF_type(-3._RKH, +3._RKH))
    call test_getQuadErr(disp, intLogNormPDF_type(lb = exp(-6._RKC), ub = exp(+6._RKC)))

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

    call test_getQuadErr(disp, int5_type(0._RKC, 3._RKC))
    call test_getQuadErr(disp, int5_type(-2._RKC, +5._RKC))

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compare integration via:")
    call disp%show("!       the Adaptive Global Gauss-Kronrod Quadrature method without extrapolation on a semi-infinite range `(lb,+Inf)`,")
    call disp%show("! with,")
    call disp%show("!       the Adaptive Global Gauss-Kronrod Quadrature method with extrapolation on a semi-infinite range `(lb,+Inf)`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call test_getQuadErr(disp, intLogNormPDF_type())
    call test_getQuadErr(disp, intDoncker1_type())
    call test_getQuadErr(disp, int6_type())
    call test_getQuadErr(disp, int7_type())

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

    call test_getQuadErr(disp, int8_type())
    call test_getQuadErr(disp, intDoncker2_type())
    call test_getQuadErr(disp, intDoncker2_type(ub = -1._RKC))

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

    call test_getQuadErr(disp, int9_type())
    call test_getQuadErr(disp, intNormPDF_type())
    call test_getQuadErr(disp, intPentaGammaInf_type())

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cauchy Principal Value over a finite interval.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call test_getQuadErr(disp, intCauchy1_type())
    call test_getQuadErr(disp, intCauchy2_type())

end program example