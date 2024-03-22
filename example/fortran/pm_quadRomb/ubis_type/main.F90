program example

    use pm_kind, only: SK, IK, SP => RKS, DP => RKD, QP => RKH ! all real kinds are supported.
    use pm_mathBeta, only: getBetaInc
    use pm_distBeta, only: getBetaPDF
    use pm_quadRomb, only: getQuadRomb, lbis_type, ubis_type
    use pm_io, only: display_type

    implicit none

    integer(IK) :: neval

    real(SP)    :: quad_SP, quadref_SP, relerr_SP, alpha_SP, beta_SP
    real(DP)    :: quad_DP, quadref_DP, relerr_DP, alpha_DP, beta_DP
    real(QP)    :: quad_QP, quadref_QP, relerr_QP, alpha_QP, beta_QP

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cumulative Distribution Function (CDF) over an open interval of the Beta distribution with singular support bounds.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the numerical integration with single precision.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("alpha_SP = .5_SP")
                    alpha_SP = .5_SP
    call disp%show("beta_SP = .5_SP")
                    beta_SP = .5_SP
    call disp%show("quadref_SP = getBetaInc(1._SP, alpha = alpha_SP, beta = beta_SP) - getBetaInc(0._SP, alpha = alpha_SP, beta = beta_SP)")
                    quadref_SP = getBetaInc(1._SP, alpha = alpha_SP, beta = beta_SP) - getBetaInc(0._SP, alpha = alpha_SP, beta = beta_SP)
    call disp%show("quadref_SP")
    call disp%show( quadref_SP )
    call disp%show("quad_SP = getQuadRomb(getBetaPDF_SP, 0._SP, .5_SP, epsilon(1.) * 100, 4_IK, lbis_type(real(alpha_SP) - 1.)) + getQuadRomb(getBetaPDF_SP, .5_SP, 1._SP, epsilon(1._SP) * 100, 4_IK, ubis_type(real(beta_SP) - 1.))")
                    quad_SP = getQuadRomb(getBetaPDF_SP, 0._SP, .5_SP, epsilon(1.) * 100, 4_IK, lbis_type(real(alpha_SP) - 1.)) + getQuadRomb(getBetaPDF_SP, .5_SP, 1._SP, epsilon(1._SP) * 100, 4_IK, ubis_type(real(beta_SP) - 1.))
    call disp%show("quad_SP ! integral")
    call disp%show( quad_SP )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the numerical integration with double precision.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("alpha_DP = .5_DP")
                    alpha_DP = .5_DP
    call disp%show("beta_DP = .5_DP")
                    beta_DP = .5_DP
    call disp%show("quadref_DP = getBetaInc(1._DP, alpha = alpha_DP, beta = beta_DP) - getBetaInc(0._DP, alpha = alpha_DP, beta = beta_DP)")
                    quadref_DP = getBetaInc(1._DP, alpha = alpha_DP, beta = beta_DP) - getBetaInc(0._DP, alpha = alpha_DP, beta = beta_DP)
    call disp%show("quadref_DP")
    call disp%show( quadref_DP )
    call disp%show("quad_DP = getQuadRomb(getBetaPDF_DP, 0._DP, .5_DP, epsilon(1._DP) * 100, 10_IK, lbis_type(real(alpha_DP) - 1.)) + getQuadRomb(getBetaPDF_DP, .5_DP, 1._DP, epsilon(1._DP) * 100, 10_IK, ubis_type(real(beta_DP) - 1.))")
                    quad_DP = getQuadRomb(getBetaPDF_DP, 0._DP, .5_DP, epsilon(1._DP) * 100, 10_IK, lbis_type(real(alpha_DP) - 1.)) + getQuadRomb(getBetaPDF_DP, .5_DP, 1._DP, epsilon(1._DP) * 100, 10_IK, ubis_type(real(beta_DP) - 1.))
    call disp%show("quad_DP ! integral")
    call disp%show( quad_DP )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the integral of 1/sqrt(abs(x)) = 4 over the interval [-1, 1] with alpha singularity at x = 0.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("quad_DP = getQuadRomb(getInvSqrtAbs_DP, -1._DP, 0._DP, epsilon(1._DP) * 100, 7_IK, ubis_type(-0.5)) + getQuadRomb(getInvSqrtAbs_DP, 0._DP, 1._DP, epsilon(1._DP) * 100, 10_IK, lbis_type(-0.5))")
                    quad_DP = getQuadRomb(getInvSqrtAbs_DP, -1._DP, 0._DP, epsilon(1._DP) * 100, 7_IK, ubis_type(-0.5)) + getQuadRomb(getInvSqrtAbs_DP, 0._DP, 1._DP, epsilon(1._DP) * 100, 10_IK, lbis_type(-0.5))
    call disp%show("quad_DP ! integral")
    call disp%show( quad_DP )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the integral of 1/sqrt(abs(x)) = 4 over the interval [-1, 1] with alpha singularity at x = 0.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("quad_QP = getQuadRomb(getInvSqrtAbs_DP, -1._QP, 0._QP, epsilon(1._QP) * 100, 7_IK, ubis_type(-0.5)) + getQuadRomb(getInvSqrtAbs_QP, 0._QP, 1._QP, epsilon(1._QP) * 100, 10_IK, lbis_type(-0.5))")
                    quad_QP = getQuadRomb(getInvSqrtAbs_QP, -1._QP, 0._QP, epsilon(1._QP) * 100, 7_IK, ubis_type(-0.5)) + getQuadRomb(getInvSqrtAbs_QP, 0._QP, 1._QP, epsilon(1._QP) * 100, 10_IK, lbis_type(-0.5))
    call disp%show("quad_QP ! integral")
    call disp%show( quad_QP )
    call disp%skip()

contains

    function getBetaPDF_SP(x) result(betaPDF)
        real    , intent(in)    :: x
        real                    :: betaPDF
        betaPDF = getBetaPDF(x, alpha = alpha_SP, beta = beta_SP)
    end function

    function getBetaPDF_DP(x) result(betaPDF)
        real(DP), intent(in)    :: x
        real(DP)                :: betaPDF
        betaPDF = getBetaPDF(x, alpha = alpha_DP, beta = beta_DP)
    end function

    function getBetaPDF_QP(x) result(betaPDF)
        real(QP), intent(in)    :: x
        real(QP)                :: betaPDF
        betaPDF = getBetaPDF(x, alpha = alpha_QP, beta = beta_QP)
    end function

    pure function getInvSqrtAbs_DP(x) result(invSqrtAbs)
        real(DP), intent(in)    :: x
        real(DP)                :: invSqrtAbs
        invSqrtAbs = 1._DP / sqrt(abs(x))
    end function

    pure function getInvSqrtAbs_QP(x) result(invSqrtAbs)
        real(QP), intent(in)    :: x
        real(QP)                :: invSqrtAbs
        invSqrtAbs = 1._QP / sqrt(abs(x))
    end function

end program example