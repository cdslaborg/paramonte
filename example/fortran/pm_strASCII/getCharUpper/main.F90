program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: getCharUpper

    implicit none

    character(1,SK), allocatable :: chr, charVec(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    chr = SK_"A"
    call disp%skip()
    call disp%show("chr")
    call disp%show( chr , deliml = SK_"""" )
    call disp%show("chr = getCharUpper(chr)")
                    chr = getCharUpper(chr)
    call disp%show("chr")
    call disp%show( chr , deliml = SK_"""" )
    call disp%skip()

    chr = SK_"a"
    call disp%skip()
    call disp%show("chr")
    call disp%show( chr , deliml = SK_"""" )
    call disp%show("chr = getCharUpper(chr)")
                    chr = getCharUpper(chr)
    call disp%show("chr")
    call disp%show( chr , deliml = SK_"""" )
    call disp%skip()

    chr = SK_" "
    call disp%skip()
    call disp%show("chr")
    call disp%show( chr , deliml = SK_"""" )
    call disp%show("chr = getCharUpper(chr)")
                    chr = getCharUpper(chr)
    call disp%show("chr")
    call disp%show( chr , deliml = SK_"""" )
    call disp%skip()

    charVec = ['a', 'A', '?', ' ']
    call disp%skip()
    call disp%show("charVec")
    call disp%show( charVec , deliml = SK_"""" )
    call disp%show("charVec = getCharUpper(charVec)")
                    charVec = getCharUpper(charVec)
    call disp%show("charVec")
    call disp%show( charVec , deliml = SK_"""" )
    call disp%skip()

end program example