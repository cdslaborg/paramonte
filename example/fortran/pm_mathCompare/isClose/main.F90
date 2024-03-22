program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RKS, RKD, RKH ! all processor types and kinds are supported.
    use pm_kind, only: CKS, CKD, CKH ! all processor types and kinds are supported.
    use pm_io, only: display_type
    use pm_except, only: getInfNeg
    use pm_except, only: getInfPos
    use pm_except, only: getNAN
    use pm_mathCompare, only: isClose
    use pm_mathCompare, only: REFERENCE, STRONG, WEAK, MEAN

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("isClose(0., [0., tiny(0.), epsilon(0.), 1.])")
    call disp%show( isClose(0., [0., tiny(0.), epsilon(0.), 1.]) )
    call disp%skip

    call disp%skip
    call disp%show("isClose(0., [0., tiny(0.), epsilon(0.), 1.], MEAN)")
    call disp%show( isClose(0., [0., tiny(0.), epsilon(0.), 1.], MEAN) )
    call disp%skip

    call disp%skip
    call disp%show("isClose(0., [0., tiny(0.), epsilon(0.), 1.], WEAK)")
    call disp%show( isClose(0., [0., tiny(0.), epsilon(0.), 1.], WEAK) )
    call disp%skip

    call disp%skip
    call disp%show("isClose(0., [0., tiny(0.), epsilon(0.), 1.], STRONG)")
    call disp%show( isClose(0., [0., tiny(0.), epsilon(0.), 1.], STRONG) )
    call disp%skip

    call disp%skip
    call disp%show("isClose(0., [0., tiny(0.), epsilon(0.), 1.], REFERENCE)")
    call disp%show( isClose(0., [0., tiny(0.), epsilon(0.), 1.], REFERENCE) )
    call disp%skip

    call disp%skip
    call disp%show("isClose(1., 1. + [0., tiny(0.), epsilon(0.), 1.])")
    call disp%show( isClose(1., 1. + [0., tiny(0.), epsilon(0.), 1.]) )
    call disp%skip

    call disp%skip
    call disp%show("isClose(1., 1. + [0., tiny(0.), epsilon(0.), 1.], MEAN)")
    call disp%show( isClose(1., 1. + [0., tiny(0.), epsilon(0.), 1.], MEAN) )
    call disp%skip

    call disp%skip
    call disp%show("isClose(1., 1. + [0., tiny(0.), epsilon(0.), 1.], WEAK)")
    call disp%show( isClose(1., 1. + [0., tiny(0.), epsilon(0.), 1.], WEAK) )
    call disp%skip

    call disp%skip
    call disp%show("isClose(1., 1. + [0., tiny(0.), epsilon(0.), 1.], STRONG)")
    call disp%show( isClose(1., 1. + [0., tiny(0.), epsilon(0.), 1.], STRONG) )
    call disp%skip

    call disp%skip
    call disp%show("isClose(1., 1. + [0., tiny(0.), epsilon(0.), 1.], REFERENCE)")
    call disp%show( isClose(1., 1. + [0., tiny(0.), epsilon(0.), 1.], REFERENCE) )
    call disp%skip

    call disp%skip
    call disp%show("isClose(getInfNeg(mold = 0.), [1., getInfNeg(0.), getInfPos(0.), getNAN(0.)])")
    call disp%show( isClose(getInfNeg(mold = 0.), [1., getInfNeg(0.), getInfPos(0.), getNAN(0.)]) )
    call disp%skip

    call disp%skip
    call disp%show("isClose(getNAN(0.), getNAN(0.))")
    call disp%show( isClose(getNAN(0.), getNAN(0.)) )
    call disp%skip

    call disp%skip
    call disp%show("isClose((0., 0.), cmplx(tiny(0.), tiny(0.)))")
    call disp%show( isClose((0., 0.), cmplx(tiny(0.), tiny(0.))) )
    call disp%skip

    call disp%skip
    call disp%show("isClose((0., 0.), cmplx(epsilon(0.), epsilon(0.)))")
    call disp%show( isClose((0., 0.), cmplx(epsilon(0.), epsilon(0.))) )
    call disp%skip

    call disp%skip
    call disp%show("isClose((1., 1.), (1., 1.) + cmplx(epsilon(0.), epsilon(0.)), reltol = 2 * epsilon(0.))")
    call disp%show( isClose((1., 1.), (1., 1.) + cmplx(epsilon(0.), epsilon(0.)), reltol = 2 * epsilon(0.)) )
    call disp%skip

    call disp%skip
    call disp%show("isClose((0., 0.), cmplx(tiny(0.), tiny(0.)), abstol = 2 * tiny(0.))")
    call disp%show( isClose((0., 0.), cmplx(tiny(0.), tiny(0.)), abstol = 2 * tiny(0.)) )
    call disp%skip

end program example