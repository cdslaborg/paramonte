program example

    use pm_kind, only: LK, SK
    use pm_kind, only: CKC => CK32 ! all complex kinds are supported.
    use pm_complexDiv, only: getDiv
    use pm_io, only: display_type

    implicit none

    complex(CKC)    :: dividend, divisor, quotient
    real(CKC)       :: large

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("dividend = (-1., -2.); divisor = (+1., +2.)")
                    dividend = (-1., -2.); divisor = (+1., +2.)
    call disp%show("quotient = getDiv(dividend, divisor)")
                    quotient = getDiv(dividend, divisor)
    call disp%show("[quotient, dividend / divisor]")
    call disp%show( [quotient, dividend / divisor] )
    call disp%skip()

    call disp%skip()
    call disp%show("large = (huge(0._CKC))/2; dividend = cmplx(large, large, CKC); divisor = cmplx(large, large, CKC)")
                    large = (huge(0._CKC))/2; dividend = cmplx(large, large, CKC); divisor = cmplx(large, large, CKC)
    call disp%show("quotient = getDiv(dividend, divisor)")
                    quotient = getDiv(dividend, divisor)
    call disp%show("[quotient, dividend / divisor]")
    call disp%show( [quotient, dividend / divisor] )
    call disp%skip()

end program example