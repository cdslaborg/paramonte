program example

    use pm_kind, only: SK, IK, QP => RKH
    use pm_mathBeta, only: getBetaInc
    use pm_distBeta, only: getBetaPDF
    use pm_quadRomb, only: getQuadRomb, open_type
    use pm_io, only: display_type

    implicit none

    integer, parameter :: SP = kind(0.e0), DP = kind(0.d0)

    integer(IK) :: neval

    real(SP)    :: quad_sp, quadref_sp, relerr_sp, alpha_sp, beta_sp
    real(DP)    :: quad_dp, quadref_dp, relerr_dp, alpha_dp, beta_dp
    real(QP)    :: quad_qp, quadref_qp, relerr_qp, alpha_qp, beta_qp

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cumulative Distribution Function (CDF) over an open interval of the Beta distribution by numerical integration.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the numerical integration with single precision.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("alpha_sp = 2.")
                    alpha_sp = 2.
    call disp%show("beta_sp = 5.")
                    beta_sp = 5.
    call disp%show("quadref_sp = getBetaInc(1., alpha = alpha_sp, beta = beta_sp) - getBetaInc(0., alpha = alpha_sp, beta = beta_sp)")
                    quadref_sp = getBetaInc(1., alpha = alpha_sp, beta = beta_sp) - getBetaInc(0., alpha = alpha_sp, beta = beta_sp)
    call disp%show("quadref_sp")
    call disp%show( quadref_sp )
    call disp%show("quad_sp = getQuadRomb(getFunc = getBetaPDF_SP, lb = 0., ub = 1., tol = epsilon(1.) * 100, nref = 4_IK, interval = open_type(), relerr = relerr_sp, neval = neval)")
                    quad_sp = getQuadRomb(getFunc = getBetaPDF_SP, lb = 0., ub = 1., tol = epsilon(1.) * 100, nref = 4_IK, interval = open_type(), relerr = relerr_sp, neval = neval)
    call disp%show("if (relerr_sp < 0.) error stop 'Integration failed to converge.'")
                    if (relerr_sp < 0.) error stop 'Integration failed to converge.'
    call disp%show("relerr_sp ! < 0. if integration fails.")
    call disp%show( relerr_sp )
    call disp%show("neval ! # calls to the integrand.")
    call disp%show( neval )
    call disp%show("quad_sp ! integral")
    call disp%show( quad_sp )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the numerical integration with double precision.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("alpha_dp = 2.d0")
                    alpha_dp = 2.d0
    call disp%show("beta_dp = 5.d0")
                    beta_dp = 5.d0
    call disp%show("quadref_dp = getBetaInc(1.d0, alpha = alpha_dp, beta = beta_dp) - getBetaInc(0.d0, alpha = alpha_dp, beta = beta_dp)")
                    quadref_dp = getBetaInc(1.d0, alpha = alpha_dp, beta = beta_dp) - getBetaInc(0.d0, alpha = alpha_dp, beta = beta_dp)
    call disp%show("quadref_dp")
    call disp%show( quadref_dp )
    call disp%show("quad_dp = getQuadRomb(getFunc = getBetaPDF_DP, lb = 0.d0, ub = 1.d0, tol = epsilon(1.d0) * 100, nref = 4_IK, interval = open_type(), relerr = relerr_dp, neval = neval)")
                    quad_dp = getQuadRomb(getFunc = getBetaPDF_DP, lb = 0.d0, ub = 1.d0, tol = epsilon(1.d0) * 100, nref = 4_IK, interval = open_type(), relerr = relerr_dp, neval = neval)
    call disp%show("if (relerr_dp < 0.) error stop 'Integration failed to converge.'")
                    if (relerr_dp < 0.) error stop 'Integration failed to converge.'
    call disp%show("relerr_dp ! < 0. if integration fails.")
    call disp%show( relerr_dp )
    call disp%show("neval ! # calls to the integrand.")
    call disp%show( neval )
    call disp%show("quad_dp ! integral")
    call disp%show( quad_dp )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the numerical integration with quadro precision.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("alpha_qp = 2._qp")
                    alpha_qp = 2._qp
    call disp%show("beta_qp = 5._qp")
                    beta_qp = 5._qp
    call disp%show("quadref_qp = getBetaInc(1._QP, alpha = alpha_qp, beta = beta_qp) - getBetaInc(0._QP, alpha = alpha_qp, beta = beta_qp)")
                    quadref_qp = getBetaInc(1._QP, alpha = alpha_qp, beta = beta_qp) - getBetaInc(0._QP, alpha = alpha_qp, beta = beta_qp)
    call disp%show("quadref_qp")
    call disp%show( quadref_qp )
    call disp%show("quad_qp = getQuadRomb(getFunc = getBetaPDF_QP, lb = 0._QP, ub = 1._QP, tol = epsilon(1._QP) * 100, nref = 4_IK, interval = open_type(), relerr = relerr_qp, neval = neval)")
                    quad_qp = getQuadRomb(getFunc = getBetaPDF_QP, lb = 0._QP, ub = 1._QP, tol = epsilon(1._QP) * 100, nref = 4_IK, interval = open_type(), relerr = relerr_qp, neval = neval)
    call disp%show("if (relerr_qp < 0.) error stop 'Integration failed to converge.'")
                    if (relerr_qp < 0.) error stop 'Integration failed to converge.'
    call disp%show("relerr_qp ! < 0. if integration fails.")
    call disp%show( relerr_qp )
    call disp%show("neval ! # calls to the integrand.")
    call disp%show( neval )
    call disp%show("quad_qp ! integral")
    call disp%show( quad_qp )
    call disp%skip()

contains

    function getBetaPDF_SP(x) result(betaPDF)
        real    , intent(in)    :: x
        real                    :: betaPDF
        betaPDF = getBetaPDF(x, alpha = alpha_sp, beta = beta_sp)
    end function

    function getBetaPDF_DP(x) result(betaPDF)
        real(DP), intent(in)    :: x
        real(DP)                :: betaPDF
        betaPDF = getBetaPDF(x, alpha = alpha_dp, beta = beta_dp)
    end function

    function getBetaPDF_QP(x) result(betaPDF)
        real(QP), intent(in)    :: x
        real(QP)                :: betaPDF
        betaPDF = getBetaPDF(x, alpha = alpha_qp, beta = beta_qp)
    end function

end program example