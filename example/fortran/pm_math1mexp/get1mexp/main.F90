program example

    use pm_kind, only: SK, IK, LK, RKH
    use pm_io, only: display_type
    use pm_arrayRank, only: getRankDense
    use pm_math1mexp, only: get1mexp, selection

    implicit none

    real :: x, result(3)
    real(RKH) :: ref
    real(RKH), allocatable :: inaccuracy(:)
    real, parameter :: EPS = epsilon(0.)
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("x = .9")
                    x = .9
    call disp%show("ref = 1._RKH - exp(real(x, RKH)) ! reference high-precision value for comparison")
                    ref = 1._RKH - exp(real(x, RKH))
    call disp%show("ref")
    call disp%show( ref )
    call disp%show("result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]")
                    result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]
    call disp%show("result")
    call disp%show( result )
    call disp%show("inaccuracy = abs(ref - result)")
                    inaccuracy = abs(ref - result)
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
    call disp%show("ref = 1._RKH - exp(real(x, RKH)) ! reference high-precision value for comparison")
                    ref = 1._RKH - exp(real(x, RKH))
    call disp%show("ref")
    call disp%show( ref )
    call disp%show("result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]")
                    result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]
    call disp%show("result")
    call disp%show( result )
    call disp%show("inaccuracy = abs(ref - result)")
                    inaccuracy = abs(ref - result)
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
    call disp%show("ref = 1._RKH - exp(real(x, RKH))) ! reference high-precision value for comparison")
                    ref = 1._RKH - exp(real(x, RKH))
    call disp%show("ref")
    call disp%show( ref )
    call disp%show("result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]")
                    result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]
    call disp%show("result")
    call disp%show( result )
    call disp%show("inaccuracy = abs(ref - result)")
                    inaccuracy = abs(ref - result)
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
    call disp%show("ref = 1._RKH - exp(real(x, RKH)) ! reference high-precision value for comparison")
                    ref = 1._RKH - exp(real(x, RKH))
    call disp%show("ref")
    call disp%show( ref )
    call disp%show("result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]")
                    result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]
    call disp%show("result")
    call disp%show( result )
    call disp%show("inaccuracy = abs(ref - result)")
                    inaccuracy = abs(ref - result)
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
    call disp%show("ref = 1._RKH - exp(real(x, RKH)) ! reference high-precision value for comparison")
                    ref = 1._RKH - exp(real(x, RKH))
    call disp%show("ref")
    call disp%show( ref )
    call disp%show("result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]")
                    result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]
    call disp%show("result")
    call disp%show( result )
    call disp%show("inaccuracy = abs(ref - result)")
                    inaccuracy = abs(ref - result)
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
    call disp%show("ref = 1._RKH - exp(real(x, RKH)) ! reference high-precision value for comparison")
                    ref = 1._RKH - exp(real(x, RKH))
    call disp%show("ref")
    call disp%show( ref )
    call disp%show("result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]")
                    result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]
    call disp%show("result")
    call disp%show( result )
    call disp%show("inaccuracy = abs(ref - result)")
                    inaccuracy = abs(ref - result)
    call disp%show("inaccuracy")
    call disp%show( inaccuracy )
    call disp%show("getRankDense(inaccuracy)")
    call disp%show( getRankDense(inaccuracy) )
    call disp%skip()

    call disp%skip()
    call disp%show("tiny(0._RKH)")
    call disp%show( tiny(0._RKH) )
    call disp%show("[(1._RKH + tiny(0._RKH)), get1mexp(tiny(0._RKH)), get1mexp(tiny(0._RKH), control = selection)]")
    call disp%show( [(1._RKH + tiny(0._RKH)), get1mexp(tiny(0._RKH)), get1mexp(tiny(0._RKH), control = selection)] )
    call disp%skip()

    call disp%skip()
    call disp%show("tiny(0._RKH)")
    call disp%show( tiny(0._RKH) )
    call disp%show("[log(1._RKH - tiny(0._RKH)), get1mexp(-tiny(0._RKH)), get1mexp(-tiny(0._RKH), control = selection)]")
    call disp%show( [log(1._RKH - tiny(0._RKH)), get1mexp(-tiny(0._RKH)), get1mexp(-tiny(0._RKH), control = selection)] )
    call disp%skip()

end program example