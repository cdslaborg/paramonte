program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: isStrAlphaAll

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isStrAlphaAll('A')")
    call disp%show( isStrAlphaAll('A') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaAll('a')")
    call disp%show( isStrAlphaAll('a') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaAll('aBcD')")
    call disp%show( isStrAlphaAll('aBcD') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaAll('ParaMonte')")
    call disp%show( isStrAlphaAll('ParaMonte') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaAll([character(4) :: 'abcd', 'ParaMonte', 'DCBA'])")
    call disp%show( isStrAlphaAll([character(4) :: 'abcd', 'ParaMonte', 'DCBA']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaAll([character(9) :: 'abcd', 'ParaMonte', 'DCBA'])")
    call disp%show( isStrAlphaAll([character(9) :: 'abcd', 'ParaMonte', 'DCBA']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaAll([character(4) :: ' ', '1aBc2', '#@!?'])")
    call disp%show( isStrAlphaAll([character(4) :: ' ', '1aBc2', '#@!?']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaAll('-a')")
    call disp%show( isStrAlphaAll('-a') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaAll('_A')")
    call disp%show( isStrAlphaAll('_A') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaAll('_')")
    call disp%show( isStrAlphaAll('_') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaAll(' ')")
    call disp%show( isStrAlphaAll(' ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlphaAll('')")
    call disp%show( isStrAlphaAll('') )
    call disp%skip()

end program example