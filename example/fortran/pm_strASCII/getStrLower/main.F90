program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: getStrLower

    implicit none

    character(9), allocatable :: Lower9(:)
    character(5), allocatable :: Lower5(:)
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getStrLower('A')")
    call disp%show( getStrLower('A') , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStrLower('a')")
    call disp%show( getStrLower('a') , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStrLower('PARAMONTE')")
    call disp%show( getStrLower('PARAMONTE') , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("Lower9 = getStrLower([character(9) :: 'abcd', 'ParaMonte', 'DCBA'])")
                    Lower9 = getStrLower([character(9) :: 'abcd', 'ParaMonte', 'DCBA'])
    call disp%show("Lower9")
    call disp%show( Lower9 , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("Lower5 = getStrLower([character(5) :: ' ', '1aBc2', '#@!?'])")
                    Lower5 = getStrLower([character(5) :: ' ', '1aBc2', '#@!?'])
    call disp%show("Lower5")
    call disp%show( Lower5 , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStrLower('-a')")
    call disp%show( getStrLower('-a') , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStrLower('_A')")
    call disp%show( getStrLower('_A') , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStrLower('_')")
    call disp%show( getStrLower('_') , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStrLower(' ')")
    call disp%show( getStrLower(' ') , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStrLower('')")
    call disp%show( getStrLower('') , deliml = """" )
    call disp%skip()

end program example