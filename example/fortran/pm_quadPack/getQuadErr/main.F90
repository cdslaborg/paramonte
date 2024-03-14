program example

    use pm_kind, only: SK, IK, RKH, RKC => RKH ! All other real kinds are also supported.
    use pm_except, only: getInfNeg, getInfPos
    use pm_quadTest, only: Int1_type
    use pm_quadTest, only: Int2_type
    use pm_quadTest, only: Int3_type
    use pm_quadTest, only: Int4_type
    use pm_quadTest, only: Int5_type
    use pm_quadTest, only: Int6_type
    use pm_quadTest, only: Int7_type
    use pm_quadTest, only: Int8_type
    use pm_quadTest, only: Int9_type
    use pm_quadTest, only: IntSinCos_type
    use pm_quadTest, only: IntNormPDF_type
    use pm_quadTest, only: IntDoncker1_type
    use pm_quadTest, only: IntDoncker2_type
    use pm_quadTest, only: IntLogNormPDF_type
    use pm_quadTest, only: IntPentaGammaInf_type
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

    call test_getQuadErr(disp, Int1_type(-4._RKH, +4._RKH))
    call test_getQuadErr(disp, Int2_type(+2._RKH, +3._RKH))
    call test_getQuadErr(disp, Int3_type(ub = +3._RKH))
    call test_getQuadErr(disp, Int4_type())
    call test_getQuadErr(disp, IntSinCos_type(-4_IK, +4_IK, a = 10._RKH, b = -10._RKH))
    call test_getQuadErr(disp, IntNormPDF_type(-3._RKH, +3._RKH))
    call test_getQuadErr(disp, IntLogNormPDF_type(lb = exp(-6._RKC), ub = exp(+6._RKC)))

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

    call test_getQuadErr(disp, Int5_type(0._RKC, 3._RKC))
    call test_getQuadErr(disp, Int5_type(-2._RKC, +5._RKC))

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

    call test_getQuadErr(disp, IntLogNormPDF_type())
    call test_getQuadErr(disp, IntDoncker1_type())
    call test_getQuadErr(disp, Int6_type())
    call test_getQuadErr(disp, Int7_type())

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

    call test_getQuadErr(disp, Int8_type())
    call test_getQuadErr(disp, IntDoncker2_type())
    call test_getQuadErr(disp, IntDoncker2_type(ub = -1._RKC))

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

    call test_getQuadErr(disp, Int9_type())
    call test_getQuadErr(disp, IntNormPDF_type())
    call test_getQuadErr(disp, IntPentaGammaInf_type())

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cauchy Principal Value over a finite interval.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call test_getQuadErr(disp, intCauchy1_type())
    call test_getQuadErr(disp, intCauchy2_type())

end program example