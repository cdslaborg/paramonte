program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: getStrUpper

    implicit none

    character(9), allocatable :: Upper9(:)
    character(5), allocatable :: Upper5(:)
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getStrUpper('A')")
    call disp%show( getStrUpper('A') , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStrUpper('a')")
    call disp%show( getStrUpper('a') , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStrUpper('paramonte')")
    call disp%show( getStrUpper('paramonte') , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("Upper9 = getStrUpper([character(9) :: 'abcd', 'ParaMonte', 'DCBA'])")
                    Upper9 = getStrUpper([character(9) :: 'abcd', 'ParaMonte', 'DCBA'])
    call disp%show("Upper9")
    call disp%show( Upper9 , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("Upper5 = getStrUpper([character(5) :: ' ', '1aBc2', '#@!?'])")
                    Upper5 = getStrUpper([character(5) :: ' ', '1aBc2', '#@!?'])
    call disp%show("Upper5")
    call disp%show( Upper5 , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStrUpper('-a')")
    call disp%show( getStrUpper('-a') , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStrUpper('_A')")
    call disp%show( getStrUpper('_A') , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStrUpper('_')")
    call disp%show( getStrUpper('_') , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStrUpper(' ')")
    call disp%show( getStrUpper(' ') , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStrUpper('')")
    call disp%show( getStrUpper('') , deliml = """" )
    call disp%skip()

end program example