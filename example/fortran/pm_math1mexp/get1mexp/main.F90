program example

    use pm_kind, only: SK, IK, LK, RK128
    use pm_io, only: display_type
    use pm_arrayRank, only: getRankDense
    use pm_math1mexp, only: get1mexp, selection

    implicit none

    real :: x, Result(3)
    real(RK128) :: ref
    real(RK128), allocatable :: Inaccuracy(:)
    real, parameter :: EPS = epsilon(0.)
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("x = .9")
                    x = .9
    call disp%show("ref = 1._RK128 - exp(real(x, RK128)) ! reference high-precision value for comparison")
                    ref = 1._RK128 - exp(real(x, RK128))
    call disp%show("ref")
    call disp%show( ref )
    call disp%show("Result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]")
                    Result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]
    call disp%show("Result")
    call disp%show( Result )
    call disp%show("Inaccuracy = abs(ref - Result)")
                    Inaccuracy = abs(ref - Result)
    call disp%show("Inaccuracy")
    call disp%show( Inaccuracy )
    call disp%show("getRankDense(Inaccuracy)")
    call disp%show( getRankDense(Inaccuracy) )
    call disp%skip()

    call disp%skip()
    call disp%show("x = EPS")
                    x = EPS
    call disp%show("x")
    call disp%show( x )
    call disp%show("ref = 1._RK128 - exp(real(x, RK128)) ! reference high-precision value for comparison")
                    ref = 1._RK128 - exp(real(x, RK128))
    call disp%show("ref")
    call disp%show( ref )
    call disp%show("Result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]")
                    Result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]
    call disp%show("Result")
    call disp%show( Result )
    call disp%show("Inaccuracy = abs(ref - Result)")
                    Inaccuracy = abs(ref - Result)
    call disp%show("Inaccuracy")
    call disp%show( Inaccuracy )
    call disp%show("getRankDense(Inaccuracy)")
    call disp%show( getRankDense(Inaccuracy) )
    call disp%skip()

    call disp%skip()
    call disp%show("x = EPS/2")
                    x = EPS/2
    call disp%show("x")
    call disp%show( x )
    call disp%show("ref = 1._RK128 - exp(real(x, RK128))) ! reference high-precision value for comparison")
                    ref = 1._RK128 - exp(real(x, RK128))
    call disp%show("ref")
    call disp%show( ref )
    call disp%show("Result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]")
                    Result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]
    call disp%show("Result")
    call disp%show( Result )
    call disp%show("Inaccuracy = abs(ref - Result)")
                    Inaccuracy = abs(ref - Result)
    call disp%show("Inaccuracy")
    call disp%show( Inaccuracy )
    call disp%show("getRankDense(Inaccuracy)")
    call disp%show( getRankDense(Inaccuracy) )
    call disp%skip()

    call disp%skip()
    call disp%show("x = -EPS/2")
                    x = -EPS/2
    call disp%show("x")
    call disp%show( x )
    call disp%show("ref = 1._RK128 - exp(real(x, RK128)) ! reference high-precision value for comparison")
                    ref = 1._RK128 - exp(real(x, RK128))
    call disp%show("ref")
    call disp%show( ref )
    call disp%show("Result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]")
                    Result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]
    call disp%show("Result")
    call disp%show( Result )
    call disp%show("Inaccuracy = abs(ref - Result)")
                    Inaccuracy = abs(ref - Result)
    call disp%show("Inaccuracy")
    call disp%show( Inaccuracy )
    call disp%show("getRankDense(Inaccuracy)")
    call disp%show( getRankDense(Inaccuracy) )
    call disp%skip()

    call disp%skip()
    call disp%show("x = sqrt(tiny(x))")
                    x = sqrt(tiny(x))
    call disp%show("x")
    call disp%show( x )
    call disp%show("ref = 1._RK128 - exp(real(x, RK128)) ! reference high-precision value for comparison")
                    ref = 1._RK128 - exp(real(x, RK128))
    call disp%show("ref")
    call disp%show( ref )
    call disp%show("Result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]")
                    Result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]
    call disp%show("Result")
    call disp%show( Result )
    call disp%show("Inaccuracy = abs(ref - Result)")
                    Inaccuracy = abs(ref - Result)
    call disp%show("Inaccuracy")
    call disp%show( Inaccuracy )
    call disp%show("getRankDense(Inaccuracy)")
    call disp%show( getRankDense(Inaccuracy) )
    call disp%skip()

    call disp%skip()
    call disp%show("x = -sqrt(tiny(x))")
                    x = -sqrt(tiny(x))
    call disp%show("x")
    call disp%show( x )
    call disp%show("ref = 1._RK128 - exp(real(x, RK128)) ! reference high-precision value for comparison")
                    ref = 1._RK128 - exp(real(x, RK128))
    call disp%show("ref")
    call disp%show( ref )
    call disp%show("Result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]")
                    Result = [1. - exp(x), get1mexp(x), get1mexp(x, control = selection)]
    call disp%show("Result")
    call disp%show( Result )
    call disp%show("Inaccuracy = abs(ref - Result)")
                    Inaccuracy = abs(ref - Result)
    call disp%show("Inaccuracy")
    call disp%show( Inaccuracy )
    call disp%show("getRankDense(Inaccuracy)")
    call disp%show( getRankDense(Inaccuracy) )
    call disp%skip()

    call disp%skip()
    call disp%show("tiny(0._RK128)")
    call disp%show( tiny(0._RK128) )
    call disp%show("[(1._RK128 + tiny(0._RK128)), get1mexp(tiny(0._RK128)), get1mexp(tiny(0._RK128), control = selection)]")
    call disp%show( [(1._RK128 + tiny(0._RK128)), get1mexp(tiny(0._RK128)), get1mexp(tiny(0._RK128), control = selection)] )
    call disp%skip()

    call disp%skip()
    call disp%show("tiny(0._RK128)")
    call disp%show( tiny(0._RK128) )
    call disp%show("[log(1._RK128 - tiny(0._RK128)), get1mexp(-tiny(0._RK128)), get1mexp(-tiny(0._RK128), control = selection)]")
    call disp%show( [log(1._RK128 - tiny(0._RK128)), get1mexp(-tiny(0._RK128)), get1mexp(-tiny(0._RK128), control = selection)] )
    call disp%skip()

end program example