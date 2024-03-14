program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: isStrUpperAll

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isStrUpperAll('A')")
    call disp%show( isStrUpperAll('A') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAll('ABCD')")
    call disp%show( isStrUpperAll('ABCD') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAll([character(1) :: 'A', 'Z', 'L'])")
    call disp%show( isStrUpperAll([character(1) :: 'A', 'Z', 'L']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAll([character(3) :: 'AAA', 'ZZZ', 'LPQ'])")
    call disp%show( isStrUpperAll([character(3) :: 'AAA', 'ZZZ', 'LPQ']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAll([character(3) :: 'A', 'ZZ', 'LPQ'])")
    call disp%show( isStrUpperAll([character(3) :: 'A', 'ZZ', 'LPQ']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAll([character(2) :: 'A ', ' Z', 'LPQ'])")
    call disp%show( isStrUpperAll([character(2) :: 'A ', ' Z', 'LPQ']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAll('A B C')")
    call disp%show( isStrUpperAll('A B C') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAll(' ')")
    call disp%show( isStrUpperAll(' ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAll('1')")
    call disp%show( isStrUpperAll('1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAll('Aa')")
    call disp%show( isStrUpperAll('Aa') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAll('AA')")
    call disp%show( isStrUpperAll('AA') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAll('aa')")
    call disp%show( isStrUpperAll('aa') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrUpperAll('?')")
    call disp%show( isStrUpperAll('?') )
    call disp%skip()

end program example