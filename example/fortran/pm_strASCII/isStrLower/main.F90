program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: isStrLower

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isStrLower('a')")
    call disp%show( isStrLower('a') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLower('abcd')")
    call disp%show( isStrLower('abcd') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLower('a b c')")
    call disp%show( isStrLower('a b c') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLower(' ')")
    call disp%show( isStrLower(' ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLower('1')")
    call disp%show( isStrLower('1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLower('Aa')")
    call disp%show( isStrLower('Aa') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLower('AA')")
    call disp%show( isStrLower('AA') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLower('aa')")
    call disp%show( isStrLower('aa') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLower('?')")
    call disp%show( isStrLower('?') )
    call disp%skip()

end program example