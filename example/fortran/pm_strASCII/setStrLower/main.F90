program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: setStrLower

    implicit none

    character(:, SK), allocatable :: str, StrVec(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    str = SK_"A"
    call disp%skip()
    call disp%show("str")
    call disp%show( str , deliml = SK_"""" )
    call disp%show("call setStrLower(str)")
                    call setStrLower(str)
    call disp%show("str")
    call disp%show( str , deliml = SK_"""" )
    call disp%skip()

    str = SK_"a"
    call disp%skip()
    call disp%show("str")
    call disp%show( str , deliml = SK_"""" )
    call disp%show("call setStrLower(str)")
                    call setStrLower(str)
    call disp%show("str")
    call disp%show( str , deliml = SK_"""" )
    call disp%skip()

    str = SK_"PARAMONTE"
    call disp%skip()
    call disp%show("str")
    call disp%show( str , deliml = SK_"""" )
    call disp%show("call setStrLower(str)")
                    call setStrLower(str)
    call disp%show("str")
    call disp%show( str , deliml = SK_"""" )
    call disp%skip()

    StrVec = [character(9,SK) :: 'abcd', 'ParaMonte', 'DCBA']
    call disp%skip()
    call disp%show("StrVec")
    call disp%show( StrVec , deliml = SK_"""" )
    call disp%show("call setStrLower(StrVec)")
                    call setStrLower(StrVec)
    call disp%show("StrVec")
    call disp%show( StrVec , deliml = SK_"""" )
    call disp%skip()

    StrVec = [character(5,SK) :: '', ' ', '1aB2_', '#@!?']
    call disp%skip()
    call disp%show("StrVec")
    call disp%show( StrVec , deliml = SK_"""" )
    call disp%show("call setStrLower(StrVec)")
                    call setStrLower(StrVec)
    call disp%show("StrVec")
    call disp%show( StrVec , deliml = SK_"""" )
    call disp%skip()

end program example