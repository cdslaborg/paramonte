program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: isStrUpperAny

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isStrUpperAny('A')")
    call disp%show( isStrUpperAny('A') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAny('ABCD')")
    call disp%show( isStrUpperAny('ABCD') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAny([character(1) :: 'A', 'Z', 'L'])")
    call disp%show( isStrUpperAny([character(1) :: 'A', 'Z', 'L']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAny([character(3) :: 'AAA', 'ZZZ', 'LPQ'])")
    call disp%show( isStrUpperAny([character(3) :: 'AAA', 'ZZZ', 'LPQ']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAny([character(3) :: 'A', 'ZZ', 'LPQ'])")
    call disp%show( isStrUpperAny([character(3) :: 'A', 'ZZ', 'LPQ']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAny([character(2) :: 'A ', ' Z', 'LPQ'])")
    call disp%show( isStrUpperAny([character(2) :: 'A ', ' Z', 'LPQ']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAny('A B C')")
    call disp%show( isStrUpperAny('A B C') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAny(' ')")
    call disp%show( isStrUpperAny(' ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAny('1')")
    call disp%show( isStrUpperAny('1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAny('Aa')")
    call disp%show( isStrUpperAny('Aa') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAny('AA')")
    call disp%show( isStrUpperAny('AA') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAny('aa')")
    call disp%show( isStrUpperAny('aa') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAny('?')")
    call disp%show( isStrUpperAny('?') )
    call disp%skip()

end program example