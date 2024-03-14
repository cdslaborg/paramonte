program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: isStrLowerAll

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isStrLowerAll('a')")
    call disp%show( isStrLowerAll('a') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAll('abcd')")
    call disp%show( isStrLowerAll('abcd') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAll([character(1,SK) :: 'a', 'z', 'l'])")
    call disp%show( isStrLowerAll([character(1,SK) :: 'a', 'z', 'l']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAll([character(3) :: 'aaa', 'zzz', 'lpq'])")
    call disp%show( isStrLowerAll([character(3) :: 'aaa', 'zzz', 'lpq']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAll([character(3) :: 'a', 'zz', 'lpq'])")
    call disp%show( isStrLowerAll([character(3) :: 'a', 'zz', 'lpq']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAll([character(2) :: 'a ', ' z', 'lpq'])")
    call disp%show( isStrLowerAll([character(2) :: 'a ', ' z', 'lpq']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAll('a b c')")
    call disp%show( isStrLowerAll('a b c') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAll(' ')")
    call disp%show( isStrLowerAll(' ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAll('1')")
    call disp%show( isStrLowerAll('1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAll('Aa')")
    call disp%show( isStrLowerAll('Aa') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAll('AA')")
    call disp%show( isStrLowerAll('AA') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAll('aa')")
    call disp%show( isStrLowerAll('aa') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrLowerAll('?')")
    call disp%show( isStrLowerAll('?') )
    call disp%skip()

end program example