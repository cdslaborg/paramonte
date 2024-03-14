program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: isStrUpper

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isStrUpper('A')")
    call disp%show( isStrUpper('A') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpper('ABCD')")
    call disp%show( isStrUpper('ABCD') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpper('A B C')")
    call disp%show( isStrUpper('A B C') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpper(' ')")
    call disp%show( isStrUpper(' ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpper('1')")
    call disp%show( isStrUpper('1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpper('Aa')")
    call disp%show( isStrUpper('Aa') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpper('AA')")
    call disp%show( isStrUpper('AA') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpper('aa')")
    call disp%show( isStrUpper('aa') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpper('?')")
    call disp%show( isStrUpper('?') )
    call disp%skip()

end program example