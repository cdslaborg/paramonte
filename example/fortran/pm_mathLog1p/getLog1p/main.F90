program example

    use pm_kind, only: SK, IK, LK, RKH
    use pm_io, only: display_type
    use pm_mathLog1p, only: getLog1p
    use pm_arrayRank, only: getRankDense

    implicit none

    real :: x, log1p(2)
    real(RKH) :: log1p_ref
    real(RKH), allocatable :: inaccuracy(:)
    real, parameter :: EPS = epsilon(0.)
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("x = .9")
                    x = .9
    call disp%show("log1p_ref = log(1._RKH + real(x, RKH)) ! reference high-precision value for comparison")
                    log1p_ref = log(1._RKH + real(x, RKH))
    call disp%show("log1p_ref")
    call disp%show( log1p_ref )
    call disp%show("log1p = [log(1. + x), getLog1p(x)]")
                    log1p = [log(1. + x), getLog1p(x)]
    call disp%show("log1p")
    call disp%show( log1p )
    call disp%show("inaccuracy = abs(log1p_ref - log1p)")
                    inaccuracy = abs(log1p_ref - log1p)
    call disp%show("inaccuracy")
    call disp%show( inaccuracy )
    call disp%show("getRankDense(inaccuracy)")
    call disp%show( getRankDense(inaccuracy) )
    call disp%skip()

    call disp%skip()
    call disp%show("x = EPS")
                    x = EPS
    call disp%show("x")
    call disp%show( x )
    call disp%show("log1p_ref = log(1._RKH + real(x, RKH)) ! reference high-precision value for comparison")
                    log1p_ref = log(1._RKH + real(x, RKH))
    call disp%show("log1p_ref")
    call disp%show( log1p_ref )
    call disp%show("log1p = [log(1. + x), getLog1p(x)]")
                    log1p = [log(1. + x), getLog1p(x)]
    call disp%show("log1p")
    call disp%show( log1p )
    call disp%show("inaccuracy = abs(log1p_ref - log1p)")
                    inaccuracy = abs(log1p_ref - log1p)
    call disp%show("inaccuracy")
    call disp%show( inaccuracy )
    call disp%show("getRankDense(inaccuracy)")
    call disp%show( getRankDense(inaccuracy) )
    call disp%skip()

    call disp%skip()
    call disp%show("x = EPS/2")
                    x = EPS/2
    call disp%show("x")
    call disp%show( x )
    call disp%show("log1p_ref = log(1._RKH + real(x, RKH)) ! reference high-precision value for comparison")
                    log1p_ref = log(1._RKH + real(x, RKH))
    call disp%show("log1p_ref")
    call disp%show( log1p_ref )
    call disp%show("log1p = [log(1. + x), getLog1p(x)]")
                    log1p = [log(1. + x), getLog1p(x)]
    call disp%show("log1p")
    call disp%show( log1p )
    call disp%show("inaccuracy = abs(log1p_ref - log1p)")
                    inaccuracy = abs(log1p_ref - log1p)
    call disp%show("inaccuracy")
    call disp%show( inaccuracy )
    call disp%show("getRankDense(inaccuracy)")
    call disp%show( getRankDense(inaccuracy) )
    call disp%skip()

    call disp%skip()
    call disp%show("x = -EPS/2")
                    x = -EPS/2
    call disp%show("x")
    call disp%show( x )
    call disp%show("log1p_ref = log(1._RKH + real(x, RKH)) ! reference high-precision value for comparison")
                    log1p_ref = log(1._RKH + real(x, RKH))
    call disp%show("log1p_ref")
    call disp%show( log1p_ref )
    call disp%show("log1p = [log(1. + x), getLog1p(x)]")
                    log1p = [log(1. + x), getLog1p(x)]
    call disp%show("log1p")
    call disp%show( log1p )
    call disp%show("inaccuracy = abs(log1p_ref - log1p)")
                    inaccuracy = abs(log1p_ref - log1p)
    call disp%show("inaccuracy")
    call disp%show( inaccuracy )
    call disp%show("getRankDense(inaccuracy)")
    call disp%show( getRankDense(inaccuracy) )
    call disp%skip()

    call disp%skip()
    call disp%show("x = sqrt(tiny(x))")
                    x = sqrt(tiny(x))
    call disp%show("x")
    call disp%show( x )
    call disp%show("log1p_ref = log(1._RKH + real(x, RKH)) ! reference high-precision value for comparison")
                    log1p_ref = log(1._RKH + real(x, RKH))
    call disp%show("log1p_ref")
    call disp%show( log1p_ref )
    call disp%show("log1p = [log(1. + x), getLog1p(x)]")
                    log1p = [log(1. + x), getLog1p(x)]
    call disp%show("log1p")
    call disp%show( log1p )
    call disp%show("inaccuracy = abs(log1p_ref - log1p)")
                    inaccuracy = abs(log1p_ref - log1p)
    call disp%show("inaccuracy")
    call disp%show( inaccuracy )
    call disp%show("getRankDense(inaccuracy)")
    call disp%show( getRankDense(inaccuracy) )
    call disp%skip()

    call disp%skip()
    call disp%show("x = -sqrt(tiny(x))")
                    x = -sqrt(tiny(x))
    call disp%show("x")
    call disp%show( x )
    call disp%show("log1p_ref = log(1._RKH + real(x, RKH)) ! reference high-precision value for comparison")
                    log1p_ref = log(1._RKH + real(x, RKH))
    call disp%show("log1p_ref")
    call disp%show( log1p_ref )
    call disp%show("log1p = [log(1. + x), getLog1p(x)]")
                    log1p = [log(1. + x), getLog1p(x)]
    call disp%show("log1p")
    call disp%show( log1p )
    call disp%show("inaccuracy = abs(log1p_ref - log1p)")
                    inaccuracy = abs(log1p_ref - log1p)
    call disp%show("inaccuracy")
    call disp%show( inaccuracy )
    call disp%show("getRankDense(inaccuracy)")
    call disp%show( getRankDense(inaccuracy) )
    call disp%skip()

    call disp%skip()
    call disp%show("tiny(0._RKH)")
    call disp%show( tiny(0._RKH) )
    call disp%show("[log(1._RKH + tiny(0._RKH)), getLog1p(tiny(0._RKH))]")
    call disp%show( [log(1._RKH + tiny(0._RKH)), getLog1p(tiny(0._RKH))] )
    call disp%skip()

    call disp%skip()
    call disp%show("tiny(0._RKH)")
    call disp%show( tiny(0._RKH) )
    call disp%show("[log(1._RKH - tiny(0._RKH)), getLog1p(-tiny(0._RKH))]")
    call disp%show( [log(1._RKH - tiny(0._RKH)), getLog1p(-tiny(0._RKH))] )
    call disp%skip()

end program example