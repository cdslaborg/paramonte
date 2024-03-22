program example

    use pm_kind, only: SK, IK, SP => RKS, DP => RKD, QP => RKH ! all real kinds are supported.
    use pm_mathGamma, only: getGammaIncLow
    use pm_distExp, only: getExpCDF
    use pm_distExp, only: getExpLogPDF
    use pm_distGamma, only: getGammaLogPDF
    use pm_quadRomb, only: getQuadRomb, pexp_type
    use pm_io, only: display_type

    implicit none

    integer(IK) :: neval

    real(SP)    :: quad_SP, quadref_SP, relerr_SP, kappa_SP, invSigma_SP
    real(DP)    :: quad_DP, quadref_DP, relerr_DP, kappa_DP, invSigma_DP
    real(QP)    :: quad_QP, quadref_QP, relerr_QP, kappa_QP, invSigma_QP

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cumulative Distribution Function (CDF) over the semi-infinite interval of Exponential distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the numerical integration with single precision.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("invSigma_SP = -3._SP")
                    invSigma_SP = -3._SP
    call disp%show("quadref_SP = getExpCDF(x = log(huge(0._SP)), invSigma = -invSigma_SP)")
                    quadref_SP = getExpCDF(x = log(huge(0._SP)), invSigma = -invSigma_SP)
    call disp%show("quadref_SP")
    call disp%show( quadref_SP )
    call disp%show("quad_SP = getQuadRomb(getFunc = getExpPDF_SP, lb = -huge(0._SP), ub = 0._SP, tol = epsilon(1._SP) * 100, nref = 4_IK, interval = pexp_type(), relerr = relerr_SP, neval = neval)")
                    quad_SP = getQuadRomb(getFunc = getExpPDF_SP, lb = -huge(0._SP), ub = 0._SP, tol = epsilon(1._SP) * 100, nref = 4_IK, interval = pexp_type(), relerr = relerr_SP, neval = neval)
    call disp%show("if (relerr_SP < 0._SP) error stop 'Integration failed to converge.'")
                    if (relerr_SP < 0._SP) error stop 'Integration failed to converge.'
    call disp%show("relerr_SP ! < 0. if integration fails.")
    call disp%show( relerr_SP )
    call disp%show("neval ! # calls to the integrand.")
    call disp%show( neval )
    call disp%show("quad_SP ! integral")
    call disp%show( quad_SP )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the numerical integration with double precision.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("invSigma_DP = -3._DP")
                    invSigma_DP = -3._DP
    call disp%show("quadref_DP = getExpCDF(x = log(huge(0._DP)), invSigma = -invSigma_DP)")
                    quadref_DP = getExpCDF(x = log(huge(0._DP)), invSigma = -invSigma_DP)
    call disp%show("quadref_DP")
    call disp%show( quadref_DP )
    call disp%show("quad_DP = getQuadRomb(getFunc = getExpPDF_DP, lb = -huge(0._DP), ub = 0._DP, tol = epsilon(1._DP) * 100, nref = 4_IK, interval = pexp_type(), relerr = relerr_DP, neval = neval)")
                    quad_DP = getQuadRomb(getFunc = getExpPDF_DP, lb = -huge(0._DP), ub = 0._DP, tol = epsilon(1._DP) * 100, nref = 4_IK, interval = pexp_type(), relerr = relerr_DP, neval = neval)
    call disp%show("if (relerr_DP < 0.) error stop 'Integration failed to converge.'")
                    if (relerr_DP < 0.) error stop 'Integration failed to converge.'
    call disp%show("relerr_DP ! < 0. if integration fails.")
    call disp%show( relerr_DP )
    call disp%show("neval ! # calls to the integrand.")
    call disp%show( neval )
    call disp%show("quad_DP ! integral")
    call disp%show( quad_DP )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the numerical integration with double precision.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("invSigma_QP = -3._QP")
                    invSigma_QP = -3._QP
    call disp%show("quadref_QP = getExpCDF(x = log(huge(0._QP)), invSigma = -invSigma_QP)")
                    quadref_QP = getExpCDF(x = log(huge(0._QP)), invSigma = -invSigma_QP)
    call disp%show("quadref_QP")
    call disp%show( quadref_QP )
    call disp%show("quad_QP = getQuadRomb(getFunc = getExpPDF_QP, lb = -huge(0._QP), ub = 0._QP, tol = epsilon(1._QP) * 100, nref = 4_IK, interval = pexp_type(), relerr = relerr_QP, neval = neval)")
                    quad_QP = getQuadRomb(getFunc = getExpPDF_QP, lb = -huge(0._QP), ub = 0._QP, tol = epsilon(1._QP) * 100, nref = 4_IK, interval = pexp_type(), relerr = relerr_QP, neval = neval)
    call disp%show("if (relerr_QP < 0.) error stop 'Integration failed to converge.'")
                    if (relerr_QP < 0.) error stop 'Integration failed to converge.'
    call disp%show("relerr_QP ! < 0. if integration fails.")
    call disp%show( relerr_QP )
    call disp%show("neval ! # calls to the integrand.")
    call disp%show( neval )
    call disp%show("quad_QP ! integral")
    call disp%show( quad_QP )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cumulative Distribution Function (CDF) over the semi-infinite interval of Gamma distribution whose tail decays exponentially.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the numerical integration with single precision.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("kappa_SP = 2._SP")
                    kappa_SP = 2._SP
    call disp%show("invSigma_SP = 2._SP")
                    invSigma_SP = 2._SP
    call disp%show("quadref_SP = getGammaIncLow(log(huge(0._SP)), kappa = kappa_SP) - getGammaIncLow(0._SP, kappa = kappa_SP)")
                    quadref_SP = getGammaIncLow(log(huge(0._SP)), kappa = kappa_SP) - getGammaIncLow(0._SP, kappa = kappa_SP)
    call disp%show("quadref_SP")
    call disp%show( quadref_SP )
    call disp%show("quad_SP = getQuadRomb(getFunc = getGammaPDF_SP, lb = -huge(0._SP), ub = 0._SP, tol = epsilon(1._SP) * 100, nref = 4_IK, interval = pexp_type(), relerr = relerr_SP, neval = neval)")
                    quad_SP = getQuadRomb(getFunc = getGammaPDF_SP, lb = -huge(0._SP), ub = 0._SP, tol = epsilon(1._SP) * 100, nref = 4_IK, interval = pexp_type(), relerr = relerr_SP, neval = neval)
    call disp%show("if (relerr_SP < 0._SP) error stop 'Integration failed to converge.'")
                    if (relerr_SP < 0._SP) error stop 'Integration failed to converge.'
    call disp%show("relerr_SP ! < 0. if integration fails.")
    call disp%show( relerr_SP )
    call disp%show("neval ! # calls to the integrand.")
    call disp%show( neval )
    call disp%show("quad_SP ! integral")
    call disp%show( quad_SP )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the numerical integration with double precision.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("kappa_DP = 2._DP")
                    kappa_DP = 2._DP
    call disp%show("invSigma_DP = 2._DP")
                    invSigma_DP = 2._DP
    call disp%show("quadref_DP = getGammaIncLow(log(huge(0._DP)), kappa = kappa_DP) - getGammaIncLow(0._DP, kappa = kappa_DP)")
                    quadref_DP = getGammaIncLow(log(huge(0._DP)), kappa = kappa_DP) - getGammaIncLow(0._DP, kappa = kappa_DP)
    call disp%show("quadref_DP")
    call disp%show( quadref_DP )
    call disp%show("quad_DP = getQuadRomb(getFunc = getGammaPDF_DP, lb = -huge(0._DP), ub = 0._DP, tol = epsilon(1._DP) * 100, nref = 7_IK, interval = pexp_type(), relerr = relerr_DP, neval = neval)")
                    quad_DP = getQuadRomb(getFunc = getGammaPDF_DP, lb = -huge(0._DP), ub = 0._DP, tol = epsilon(1._DP) * 100, nref = 7_IK, interval = pexp_type(), relerr = relerr_DP, neval = neval)
    call disp%show("if (relerr_DP < 0._DP) error stop 'Integration failed to converge.'")
                    if (relerr_DP < 0._DP) error stop 'Integration failed to converge.'
    call disp%show("relerr_DP ! < 0. if integration fails.")
    call disp%show( relerr_DP )
    call disp%show("neval ! # calls to the integrand.")
    call disp%show( neval )
    call disp%show("quad_DP ! integral")
    call disp%show( quad_DP )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the numerical integration with quadro precision.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("kappa_QP = 2._QP")
                    kappa_QP = 2._QP
    call disp%show("invSigma_QP = 2._QP")
                    invSigma_QP = 2._QP
    call disp%show("quadref_QP = getGammaIncLow(log(huge(0._QP)), kappa = kappa_QP) - getGammaIncLow(0._QP, kappa = kappa_QP)")
                    quadref_QP = getGammaIncLow(log(huge(0._QP)), kappa = kappa_QP) - getGammaIncLow(0._QP, kappa = kappa_QP)
    call disp%show("quadref_QP")
    call disp%show( quadref_QP )
    call disp%show("quad_QP = getQuadRomb(getFunc = getGammaPDF_QP, lb = -huge(0._QP), ub = 0._QP, tol = sqrt(epsilon(1._QP)), nref = 10_IK, interval = pexp_type(), relerr = relerr_QP, neval = neval)")
                    quad_QP = getQuadRomb(getFunc = getGammaPDF_QP, lb = -huge(0._QP), ub = 0._QP, tol = sqrt(epsilon(1._QP)), nref = 10_IK, interval = pexp_type(), relerr = relerr_QP, neval = neval)
    call disp%show("if (relerr_QP < 0._QP) error stop 'Integration failed to converge.'")
                    if (relerr_QP < 0._QP) error stop 'Integration failed to converge.'
    call disp%show("relerr_QP ! < 0. if integration fails.")
    call disp%show( relerr_QP )
    call disp%show("neval ! # calls to the integrand.")
    call disp%show( neval )
    call disp%show("quad_QP ! integral")
    call disp%show( quad_QP )
    call disp%skip()

contains

    function getExpPDF_SP(x) result(expPDF)
        real(SP), intent(in)    :: x
        real(SP)                :: expPDF
        expPDF = exp(getExpLogPDF(-x, invSigma = -invSigma_SP))
    end function

    function getExpPDF_DP(x) result(expPDF)
        real(DP), intent(in)    :: x
        real(DP)                :: expPDF
        expPDF = exp(getExpLogPDF(-x, invSigma = -invSigma_DP))
    end function

    function getExpPDF_QP(x) result(expPDF)
        real(QP), intent(in)    :: x
        real(QP)                :: expPDF
        expPDF = exp(getExpLogPDF(-x, invSigma = -invSigma_QP))
    end function

    function getGammaPDF_SP(x) result(gammaPDF)
        real(SP), intent(in)    :: x
        real(SP)                :: gammaPDF
        gammaPDF = exp(getGammaLogPDF(-x, kappa = kappa_SP, invSigma = invSigma_SP))
    end function

    function getGammaPDF_DP(x) result(gammaPDF)
        real(DP), intent(in)    :: x
        real(DP)                :: gammaPDF
        gammaPDF = exp(getGammaLogPDF(-x, kappa = kappa_DP, invSigma = invSigma_DP))
    end function

    function getGammaPDF_QP(x) result(gammaPDF)
        real(QP), intent(in)    :: x
        real(QP)                :: gammaPDF
        gammaPDF = exp(getGammaLogPDF(-x, kappa = kappa_QP, invSigma = invSigma_QP))
    end function

end program example