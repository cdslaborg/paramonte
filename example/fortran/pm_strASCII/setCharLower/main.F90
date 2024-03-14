program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: setCharLower

    implicit none

    character(1,SK), allocatable :: chr, StrVec(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    chr = SK_"A"
    call disp%skip()
    call disp%show("chr")
    call disp%show( chr , deliml = SK_"""" )
    call disp%show("call setCharLower(chr)")
                    call setCharLower(chr)
    call disp%show("chr")
    call disp%show( chr , deliml = SK_"""" )
    call disp%skip()

    chr = SK_"a"
    call disp%skip()
    call disp%show("chr")
    call disp%show( chr , deliml = SK_"""" )
    call disp%show("call setCharLower(chr)")
                    call setCharLower(chr)
    call disp%show("chr")
    call disp%show( chr , deliml = SK_"""" )
    call disp%skip()

    chr = SK_" "
    call disp%skip()
    call disp%show("chr")
    call disp%show( chr , deliml = SK_"""" )
    call disp%show("call setCharLower(chr)")
                    call setCharLower(chr)
    call disp%show("chr")
    call disp%show( chr , deliml = SK_"""" )
    call disp%skip()

    StrVec = ['a', 'A', '?', ' ']
    call disp%skip()
    call disp%show("StrVec")
    call disp%show( StrVec , deliml = SK_"""" )
    call disp%show("call setCharLower(StrVec)")
                    call setCharLower(StrVec)
    call disp%show("StrVec")
    call disp%show( StrVec , deliml = SK_"""" )
    call disp%skip()

end program example