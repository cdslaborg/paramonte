program example

    use pm_kind, only: SK, IK, RK => RKH ! All other real kinds are also supported.
    use pm_quadpack, only: getQuadGK, ninf, pinf
    use pm_quadpack, only: setNodeWeightGK
    use pm_quadpack, only: GK15, GK21
    use pm_quadpack, only: GK31, GK41
    use pm_quadpack, only: GK51, GK61
    use pm_io, only: display_type

    implicit none

    real(RK)            :: integral, abserr, intAbsFunc, smoothness, truth

    type(display_type)  :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Gauss-Kronrod quadrature over a finite interval.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Gauss-Kronrod quadrature of the probability density function of the Beta distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    truth = 1._RK

    call disp%skip()
    call disp%show("integral = getQuadGK(getBetaPDF, 0._RK, +1._RK, GK15, abserr, intAbsFunc, smoothness)")
                    integral = getQuadGK(getBetaPDF, 0._RK, +1._RK, GK15, abserr, intAbsFunc, smoothness)
    call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
    call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
    call disp%skip()

    call disp%skip()
    call disp%show("integral = getQuadGK(getBetaPDF, 0._RK, +1._RK, GK21, abserr, intAbsFunc, smoothness)")
                    integral = getQuadGK(getBetaPDF, 0._RK, +1._RK, GK21, abserr, intAbsFunc, smoothness)
    call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
    call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
    call disp%skip()

    call disp%skip()
    call disp%show("integral = getQuadGK(getBetaPDF, 0._RK, +1._RK, GK41, abserr, intAbsFunc, smoothness)")
                    integral = getQuadGK(getBetaPDF, 0._RK, +1._RK, GK41, abserr, intAbsFunc, smoothness)
    call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
    call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
    call disp%skip()

    call disp%skip()
    call disp%show("integral = getQuadGK(getBetaPDF, 0._RK, +1._RK, GK51, abserr, intAbsFunc, smoothness)")
                    integral = getQuadGK(getBetaPDF, 0._RK, +1._RK, GK51, abserr, intAbsFunc, smoothness)
    call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
    call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
    call disp%skip()

    call disp%skip()
    call disp%show("integral = getQuadGK(getBetaPDF, 0._RK, +1._RK, GK61, abserr, intAbsFunc, smoothness)")
                    integral = getQuadGK(getBetaPDF, 0._RK, +1._RK, GK61, abserr, intAbsFunc, smoothness)
    call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
    call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Gauss-Kronrod quadrature of the probability density function of the Beta distribution for non-default node counts.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        integer(IK) , parameter :: MAX_NPG = 30_IK
        integer(IK)             :: npg
        real(RK)                :: nodeK(MAX_NPG+1), weightK(MAX_NPG+1), weightG((MAX_NPG+1)/2)
        do npg = 1_IK, MAX_NPG
            call disp%skip()
            call disp%show("npg ! The number of points in the Gauss-Legendre quadrature rule.")
            call disp%show( npg )
            call disp%show("call setNodeWeightGK(nodeK(1:npg+1), weightK(1:npg+1), weightG(1:(npg+1)/2)) ! Compute the Gauss-Kronrod nodes and weights")
                            call setNodeWeightGK(nodeK(1:npg+1), weightK(1:npg+1), weightG(1:(npg+1)/2))
            call disp%show("integral = getQuadGK(getBetaPDF, 0._RK, +1._RK, nodeK(1:npg+1), weightK(1:npg+1), weightG(1:(npg+1)/2), abserr, intAbsFunc, smoothness) ! Compute the (2*npg + 1)-points Gauss-Kronrod quadrature.")
                            integral = getQuadGK(getBetaPDF, 0._RK, +1._RK, nodeK(1:npg+1), weightK(1:npg+1), weightG(1:(npg+1)/2), abserr, intAbsFunc, smoothness)
            call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
            call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
            call disp%skip()
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Gauss-Kronrod quadrature over a positive semi-infinite interval.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Gauss-Kronrod quadrature of the probability density function of the LogNormal distribution on the semi-infinite range `(0,+infinity)`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    truth = 1._RK

    call disp%skip()
    call disp%show("integral = getQuadGK(getLogNormPDF, 0._RK, pinf, GK15, abserr, intAbsFunc, smoothness)")
                    integral = getQuadGK(getLogNormPDF, 0._RK, pinf, GK15, abserr, intAbsFunc, smoothness)
    call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
    call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
    call disp%skip()

    call disp%skip()
    call disp%show("integral = getQuadGK(getLogNormPDF, 0._RK, pinf, GK21, abserr, intAbsFunc, smoothness)")
                    integral = getQuadGK(getLogNormPDF, 0._RK, pinf, GK21, abserr, intAbsFunc, smoothness)
    call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
    call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
    call disp%skip()

    call disp%skip()
    call disp%show("integral = getQuadGK(getLogNormPDF, 0._RK, pinf, GK41, abserr, intAbsFunc, smoothness)")
                    integral = getQuadGK(getLogNormPDF, 0._RK, pinf, GK41, abserr, intAbsFunc, smoothness)
    call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
    call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
    call disp%skip()

    call disp%skip()
    call disp%show("integral = getQuadGK(getLogNormPDF, 0._RK, pinf, GK51, abserr, intAbsFunc, smoothness)")
                    integral = getQuadGK(getLogNormPDF, 0._RK, pinf, GK51, abserr, intAbsFunc, smoothness)
    call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
    call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
    call disp%skip()

    call disp%skip()
    call disp%show("integral = getQuadGK(getLogNormPDF, 0._RK, pinf, GK61, abserr, intAbsFunc, smoothness)")
                    integral = getQuadGK(getLogNormPDF, 0._RK, pinf, GK61, abserr, intAbsFunc, smoothness)
    call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
    call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Gauss-Kronrod quadrature of the probability density function of the LogNormal distribution for non-default node counts.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        integer(IK) , parameter :: MAX_NPG = 30_IK
        integer(IK)             :: npg
        real(RK)                :: nodeK(MAX_NPG+1), weightK(MAX_NPG+1), weightG((MAX_NPG+1)/2)
        do npg = 1_IK, MAX_NPG
            call disp%skip()
            call disp%show("npg ! The number of points in the Gauss-Legendre quadrature rule.")
            call disp%show( npg )
            call disp%show("call setNodeWeightGK(nodeK(1:npg+1), weightK(1:npg+1), weightG(1:(npg+1)/2)) ! Compute the Gauss-Kronrod nodes and weights")
                            call setNodeWeightGK(nodeK(1:npg+1), weightK(1:npg+1), weightG(1:(npg+1)/2))
            call disp%show("integral = getQuadGK(getLogNormPDF, 0._RK, pinf, nodeK(1:npg+1), weightK(1:npg+1), weightG(1:(npg+1)/2), abserr, intAbsFunc, smoothness) ! Compute the (2*npg + 1)-points Gauss-Kronrod quadrature.")
                            integral = getQuadGK(getLogNormPDF, 0._RK, pinf, nodeK(1:npg+1), weightK(1:npg+1), weightG(1:(npg+1)/2), abserr, intAbsFunc, smoothness)
            call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
            call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
            call disp%skip()
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Gauss-Kronrod quadrature over a negative semi-infinite interval.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Gauss-Kronrod quadrature of the probability density function of the Exponential distribution on the semi-infinite range `(-infinity,0)`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    truth = 1._RK

    call disp%skip()
    call disp%show("integral = getQuadGK(getNegExpPDF, ninf, 0._RK, GK15, abserr, intAbsFunc, smoothness)")
                    integral = getQuadGK(getNegExpPDF, ninf, 0._RK, GK15, abserr, intAbsFunc, smoothness)
    call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
    call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
    call disp%skip()

    call disp%skip()
    call disp%show("integral = getQuadGK(getNegExpPDF, ninf, 0._RK, GK21, abserr, intAbsFunc, smoothness)")
                    integral = getQuadGK(getNegExpPDF, ninf, 0._RK, GK21, abserr, intAbsFunc, smoothness)
    call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
    call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
    call disp%skip()

    call disp%skip()
    call disp%show("integral = getQuadGK(getNegExpPDF, ninf, 0._RK, GK41, abserr, intAbsFunc, smoothness)")
                    integral = getQuadGK(getNegExpPDF, ninf, 0._RK, GK41, abserr, intAbsFunc, smoothness)
    call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
    call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
    call disp%skip()

    call disp%skip()
    call disp%show("integral = getQuadGK(getNegExpPDF, ninf, 0._RK, GK51, abserr, intAbsFunc, smoothness)")
                    integral = getQuadGK(getNegExpPDF, ninf, 0._RK, GK51, abserr, intAbsFunc, smoothness)
    call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
    call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
    call disp%skip()

    call disp%skip()
    call disp%show("integral = getQuadGK(getNegExpPDF, ninf, 0._RK, GK61, abserr, intAbsFunc, smoothness)")
                    integral = getQuadGK(getNegExpPDF, ninf, 0._RK, GK61, abserr, intAbsFunc, smoothness)
    call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
    call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Gauss-Kronrod quadrature of the probability density function of the Exponential distribution for non-default node counts.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        integer(IK) , parameter :: MAX_NPG = 30_IK
        integer(IK)             :: npg
        real(RK)                :: nodeK(MAX_NPG+1), weightK(MAX_NPG+1), weightG((MAX_NPG+1)/2)
        do npg = 1_IK, MAX_NPG
            call disp%skip()
            call disp%show("npg ! The number of points in the Gauss-Legendre quadrature rule.")
            call disp%show( npg )
            call disp%show("call setNodeWeightGK(nodeK(1:npg+1), weightK(1:npg+1), weightG(1:(npg+1)/2)) ! Compute the Gauss-Kronrod nodes and weights")
                            call setNodeWeightGK(nodeK(1:npg+1), weightK(1:npg+1), weightG(1:(npg+1)/2))
            call disp%show("integral = getQuadGK(getNegExpPDF, ninf, 0._RK, nodeK(1:npg+1), weightK(1:npg+1), weightG(1:(npg+1)/2), abserr, intAbsFunc, smoothness) ! Compute the (2*npg + 1)-points Gauss-Kronrod quadrature.")
                            integral = getQuadGK(getNegExpPDF, ninf, 0._RK, nodeK(1:npg+1), weightK(1:npg+1), weightG(1:(npg+1)/2), abserr, intAbsFunc, smoothness)
            call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
            call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
            call disp%skip()
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Gauss-Kronrod quadrature over a full infinite interval (-infinity, +infinity).")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Gauss-Kronrod quadrature of the probability density function of the Normal distribution on the infinite range `(-infinity,+infinity)`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    truth = 1._RK

    call disp%skip()
    call disp%show("integral = getQuadGK(getNormPDF, ninf, pinf, GK15, abserr, intAbsFunc, smoothness)")
                    integral = getQuadGK(getNormPDF, ninf, pinf, GK15, abserr, intAbsFunc, smoothness)
    call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
    call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
    call disp%skip()

    call disp%skip()
    call disp%show("integral = getQuadGK(getNormPDF, ninf, pinf, GK21, abserr, intAbsFunc, smoothness)")
                    integral = getQuadGK(getNormPDF, ninf, pinf, GK21, abserr, intAbsFunc, smoothness)
    call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
    call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
    call disp%skip()

    call disp%skip()
    call disp%show("integral = getQuadGK(getNormPDF, ninf, pinf, GK41, abserr, intAbsFunc, smoothness)")
                    integral = getQuadGK(getNormPDF, ninf, pinf, GK41, abserr, intAbsFunc, smoothness)
    call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
    call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
    call disp%skip()

    call disp%skip()
    call disp%show("integral = getQuadGK(getNormPDF, ninf, pinf, GK51, abserr, intAbsFunc, smoothness)")
                    integral = getQuadGK(getNormPDF, ninf, pinf, GK51, abserr, intAbsFunc, smoothness)
    call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
    call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
    call disp%skip()

    call disp%skip()
    call disp%show("integral = getQuadGK(getNormPDF, ninf, pinf, GK61, abserr, intAbsFunc, smoothness)")
                    integral = getQuadGK(getNormPDF, ninf, pinf, GK61, abserr, intAbsFunc, smoothness)
    call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
    call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Gauss-Kronrod quadrature of the probability density function of the Normal distribution for non-default node counts.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        integer(IK) , parameter :: MAX_NPG = 30_IK
        integer(IK)             :: npg
        real(RK)                :: nodeK(MAX_NPG+1), weightK(MAX_NPG+1), weightG((MAX_NPG+1)/2)
        do npg = 1_IK, MAX_NPG
            call disp%skip()
            call disp%show("npg ! The number of points in the Gauss-Legendre quadrature rule.")
            call disp%show( npg )
            call disp%show("call setNodeWeightGK(nodeK(1:npg+1), weightK(1:npg+1), weightG(1:(npg+1)/2)) ! Compute the Gauss-Kronrod nodes and weights")
                            call setNodeWeightGK(nodeK(1:npg+1), weightK(1:npg+1), weightG(1:(npg+1)/2))
            call disp%show("integral = getQuadGK(getNormPDF, ninf, pinf, nodeK(1:npg+1), weightK(1:npg+1), weightG(1:(npg+1)/2), abserr, intAbsFunc, smoothness) ! Compute the (2*npg + 1)-points Gauss-Kronrod quadrature.")
                            integral = getQuadGK(getNormPDF, ninf, pinf, nodeK(1:npg+1), weightK(1:npg+1), weightG(1:(npg+1)/2), abserr, intAbsFunc, smoothness)
            call disp%show("[truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)]")
            call disp%show( [truth, integral, abserr, intAbsFunc, smoothness, abs((integral - truth) / truth)] )
            call disp%skip()
        end do
    end block

contains

    function getBetaPDF(x) result(pdf)
        use pm_distBeta, only: getBetaLogPDF
        real(RK)    , intent(in)    :: x
        real(RK)                    :: pdf
        pdf = exp(getBetaLogPDF(x, alpha = 2._RK, beta = 2._RK))
    end function

    function getLogNormPDF(x) result(pdf)
        use pm_distLogNorm, only: getLogNormLogPDF
        real(RK)    , intent(in)    :: x
        real(RK)                    :: pdf
        pdf = exp(getLogNormLogPDF(x))
    end function

    function getNegExpPDF(x) result(pdf)
        use pm_distExp, only: getExpLogPDF
        real(RK)    , intent(in)    :: x
        real(RK)                    :: pdf
        pdf = exp(getExpLogPDF(-x))
    end function

    function getNormPDF(x) result(pdf)
        use pm_distNorm, only: getNormLogPDF
        real(RK)    , intent(in)    :: x
        real(RK)                    :: pdf
        pdf = exp(getNormLogPDF(x))
    end function

end program example