program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: isStrLowerAny

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isStrLowerAny('a')")
    call disp%show( isStrLowerAny('a') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAny('abcd')")
    call disp%show( isStrLowerAny('abcd') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAny([character(1) :: 'a', 'z', 'l'])")
    call disp%show( isStrLowerAny([character(1) :: 'a', 'z', 'l']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAny([character(3) :: 'aaa', 'zzz', 'lpq'])")
    call disp%show( isStrLowerAny([character(3) :: 'aaa', 'zzz', 'lpq']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAny([character(3) :: 'a', 'zz', 'lpq'])")
    call disp%show( isStrLowerAny([character(3) :: 'a', 'zz', 'lpq']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAny([character(2) :: 'a ', ' z', 'lpq'])")
    call disp%show( isStrLowerAny([character(2) :: 'a ', ' z', 'lpq']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAny('a b c')")
    call disp%show( isStrLowerAny('a b c') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAny(' ')")
    call disp%show( isStrLowerAny(' ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAny('1')")
    call disp%show( isStrLowerAny('1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAny('Aa')")
    call disp%show( isStrLowerAny('Aa') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAny('AA')")
    call disp%show( isStrLowerAny('AA') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAny('aa')")
    call disp%show( isStrLowerAny('aa') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAny('?')")
    call disp%show( isStrLowerAny('?') )
    call disp%skip()

end program example