program example

    use pm_kind, only: LK, SK
    use pm_kind, only: CKG => CKS ! all complex kinds are supported.
    use pm_complexDiv, only: getDiv
    use pm_io, only: display_type

    implicit none

    complex(CKG)    :: dividend, divisor, quotient
    real(CKG)       :: large

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
    call disp%show("large = (huge(0._CKG))/2; dividend = cmplx(large, large, CKG); divisor = cmplx(large, large, CKG)")
                    large = (huge(0._CKG))/2; dividend = cmplx(large, large, CKG); divisor = cmplx(large, large, CKG)
    call disp%show("quotient = getDiv(dividend, divisor)")
                    quotient = getDiv(dividend, divisor)
    call disp%show("[quotient, dividend / divisor]")
    call disp%show( [quotient, dividend / divisor] )
    call disp%skip()

end program example