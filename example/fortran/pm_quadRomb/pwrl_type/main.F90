program example

    use pm_kind, only: SK, IK, SP => RK32, DP => RK64, QP => RK128 ! all real kinds are supported.
    use pm_distPareto, only: getParetoLogCDF
    use pm_distPareto, only: getParetoLogPDF
    use pm_distNorm, only: getNormLogPDF
    use pm_quadRomb, only: getQuadRomb, pwrl_type
    use pm_io, only: display_type

    implicit none

    integer(IK) :: neval

    real(SP)    :: quad_SP, quadref_SP, relerr_SP, avg_SP, std_SP, logMinX_SP, alpha_SP
    real(DP)    :: quad_DP, quadref_DP, relerr_DP, avg_DP, std_DP, logMinX_DP, alpha_DP
    real(QP)    :: quad_QP, quadref_QP, relerr_QP, avg_QP, std_QP, logMinX_QP, alpha_QP

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cumulative Distribution Function (CDF) over the semi-infinite interval of Pareto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the numerical integration with single precision.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logMinX_SP = 3._SP; alpha_SP = -1._SP")
                    logMinX_SP = 3._SP; alpha_SP = -1._SP
    call disp%show("quadref_SP = getParetoLogCDF(logx = log(huge(0._SP)), logMinX = logMinX_SP, alpha = alpha_SP)")
                    quadref_SP = getParetoLogCDF(logx = log(huge(0._SP)), logMinX = logMinX_SP, alpha = alpha_SP)
    call disp%show("quadref_SP")
    call disp%show( quadref_SP )
    call disp%show("quad_SP = getQuadRomb(getFunc = getParetoPDF_SP, lb = logMinX_SP, ub = huge(0._SP), tol = epsilon(1._SP) * 100, nref = 4_IK, interval = pwrl_type(), relerr = relerr_SP, neval = neval)")
                    quad_SP = getQuadRomb(getFunc = getParetoPDF_SP, lb = logMinX_SP, ub = huge(0._SP), tol = epsilon(1._SP) * 100, nref = 4_IK, interval = pwrl_type(), relerr = relerr_SP, neval = neval)
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
    call disp%show("logMinX_DP = 3._DP; alpha_DP = -1._DP")
                    logMinX_DP = 3._DP; alpha_DP = -1._DP
    call disp%show("quadref_DP = getParetoLogCDF(logx = log(huge(0._DP)), logMinX = logMinX_DP, alpha = alpha_DP)")
                    quadref_DP = getParetoLogCDF(logx = log(huge(0._DP)), logMinX = logMinX_DP, alpha = alpha_DP)
    call disp%show("quadref_DP")
    call disp%show( quadref_DP )
    call disp%show("quad_DP = getQuadRomb(getFunc = getParetoPDF_DP, lb = logMinX_DP, ub = huge(0._DP), tol = epsilon(1._DP) * 100, nref = 7_IK, interval = pwrl_type(), relerr = relerr_DP, neval = neval)")
                    quad_DP = getQuadRomb(getFunc = getParetoPDF_DP, lb = logMinX_DP, ub = huge(0._DP), tol = epsilon(1._DP) * 100, nref = 7_IK, interval = pwrl_type(), relerr = relerr_DP, neval = neval)
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
    call disp%show("logMinX_QP = 1._QP; alpha_QP = -3._QP")
                    logMinX_QP = 1._QP; alpha_QP = -3._QP
    call disp%show("quadref_QP = getParetoLogCDF(logx = log(huge(0._QP)), logMinX = logMinX_QP, alpha = alpha_QP)")
                    quadref_QP = getParetoLogCDF(logx = log(huge(0._QP)), logMinX = logMinX_QP, alpha = alpha_QP)
    call disp%show("quadref_QP")
    call disp%show( quadref_QP )
    call disp%show("quad_QP = getQuadRomb(getFunc = getParetoPDF_QP, lb = logMinX_QP, ub = huge(0._QP), tol = epsilon(1._QP)*100, nref = 10_IK, interval = pwrl_type(), relerr = relerr_QP, neval = neval)")
                    quad_QP = getQuadRomb(getFunc = getParetoPDF_QP, lb = logMinX_QP, ub = huge(0._QP), tol = epsilon(1._QP)*100, nref = 10_IK, interval = pwrl_type(), relerr = relerr_QP, neval = neval)
    call disp%show("if (relerr_QP < 0.) error stop 'Integration failed to converge.'")
                    if (relerr_QP < 0.) error stop 'Integration failed to converge.'
    call disp%show("relerr_QP ! < 0. if integration fails.")
    call disp%show( relerr_QP )
    call disp%show("neval ! # calls to the integrand.")
    call disp%show( neval )
    call disp%show("quad_QP ! integral")
    call disp%show( quad_QP )
    call disp%skip()

contains

    function getParetoPDF_SP(x) result(expPDF)
        real(SP), intent(in)    :: x
        real(SP)                :: expPDF
        expPDF = exp(getParetoLogPDF(x, alpha = alpha_SP, logMinX = logMinX_SP))
    end function

    function getParetoPDF_DP(x) result(expPDF)
        real(DP), intent(in)    :: x
        real(DP)                :: expPDF
        expPDF = exp(getParetoLogPDF(x, alpha = alpha_DP, logMinX = logMinX_DP))
    end function

    function getParetoPDF_QP(x) result(expPDF)
        real(QP), intent(in)    :: x
        real(QP)                :: expPDF
        expPDF = exp(getParetoLogPDF(x, alpha = alpha_QP, logMinX = logMinX_QP))
    end function

    function getNormPDF_SP(x) result(expPDF)
        real(SP), intent(in)    :: x
        real(SP)                :: expPDF
        expPDF = exp(getNormLogPDF(x, mu = logMinX_SP))
    end function

    function getNormPDF_DP(x) result(expPDF)
        real(DP), intent(in)    :: x
        real(DP)                :: expPDF
        expPDF = exp(getNormLogPDF(x, mu = logMinX_DP))
    end function

    function getNormPDF_QP(x) result(expPDF)
        real(QP), intent(in)    :: x
        real(QP)                :: expPDF
        expPDF = exp(getNormLogPDF(x, mu = logMinX_QP))
    end function

end program example