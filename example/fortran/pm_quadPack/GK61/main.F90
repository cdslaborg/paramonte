program example

    use pm_kind, only: SK, IK, RK => RKH ! All other real kinds are also supported.
    use pm_quadpack, only: getQuadErr, GK61
    use pm_quadpack, only: setNodeWeightGK
    use pm_distLogNorm, only: getLogNormCDF
    use pm_distNorm, only: getNormCDF
    use pm_io, only: display_type
    use pm_val2str, only: getStr

    implicit none

    integer(IK) , parameter     :: NINTMAX = 1000_IK
    real(RK)    , parameter     :: ABSTOL = 0._RK, RELTOL = epsilon(0._RK) * 10000
    real(RK)                    :: lb, ub, truth, integral, abserr, sinfo(4, NINTMAX)
    integer(IK)                 :: err, neval, nint, sindex(NINTMAX)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Gauss-Kronrod quadrature of the probability density function of the Standard Normal distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("[ABSTOL, RELTOL]")
    call disp%show( [ABSTOL, RELTOL] )
    call disp%show("lb = -100._RK; ub = +100._RK")
                    lb = -100._RK; ub = +100._RK
    call disp%show("err = getQuadErr(getNormPDF, lb, ub, ABSTOL, RELTOL, GK61, integral, abserr, sinfo, sindex, neval, nint)")
                    err = getQuadErr(getNormPDF, lb, ub, ABSTOL, RELTOL, GK61, integral, abserr, sinfo, sindex, neval, nint)
    call disp%show("if (err /= 0_IK) error stop 'integration failed: err = '//getStr(err)")
                    if (err /= 0_IK) error stop 'integration failed: err = '//getStr(err)
    call disp%show("truth = getNormCDF(ub) - getNormCDF(lb)")
                    truth = getNormCDF(ub) - getNormCDF(lb)
    call disp%show("[truth, integral, abserr, abs(integral - truth) / truth]")
    call disp%show( [truth, integral, abserr, abs(integral - truth) / truth] )
    call disp%show("[neval, nint]")
    call disp%show( [neval, nint] )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Gauss-Kronrod quadrature of the probability density function of the Standard LogNormal distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("[ABSTOL, RELTOL]")
    call disp%show( [ABSTOL, RELTOL] )
    call disp%show("lb = exp(-10._RK); ub = exp(+10._RK)")
                    lb = exp(-10._RK); ub = exp(+10._RK)
    call disp%show("err = getQuadErr(getLogNormPDF, lb, ub, ABSTOL, RELTOL, GK61, integral, abserr, sinfo, sindex, neval, nint)")
                    err = getQuadErr(getLogNormPDF, lb, ub, ABSTOL, RELTOL, GK61, integral, abserr, sinfo, sindex, neval, nint)
    call disp%show("if (err /= 0_IK) error stop 'integration failed: err = '//getStr(err)")
                    if (err /= 0_IK) error stop 'integration failed: err = '//getStr(err)
    call disp%show("truth = getLogNormCDF(ub) - getLogNormCDF(lb)")
                    truth = getLogNormCDF(ub) - getLogNormCDF(lb)
    call disp%show("[truth, integral, abserr, abs(integral - truth) / truth]")
    call disp%show( [truth, integral, abserr, abs(integral - truth) / truth] )
    call disp%show("[neval, nint]")
    call disp%show( [neval, nint] )
    call disp%skip()

contains

    function getIntSinCos(x) result(integrand)
        real(RK)    , intent(in)    :: x
        real(RK)                    :: integrand
        integrand = cos(100._RK * sin(x))
    end function

    function getNormPDF(x) result(pdf)
        use pm_distNorm, only: getNormLogPDF
        real(RK)    , intent(in)    :: x
        real(RK)                    :: pdf
        pdf = exp(getNormLogPDF(x))
    end function

    function getLogNormPDF(x) result(pdf)
        use pm_distLogNorm, only: getLogNormLogPDF
        real(RK)    , intent(in)    :: x
        real(RK)                    :: pdf
        pdf = exp(getLogNormLogPDF(x))
    end function

end program example