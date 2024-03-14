program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: getLocNonSpace

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getLocNonSpace(SK_'')")
    call disp%show( getLocNonSpace(SK_'') )
    call disp%skip()

    call disp%skip()
    call disp%show("getLocNonSpace(SK_' ')")
    call disp%show( getLocNonSpace(SK_' ') )
    call disp%skip()

    call disp%skip()
    call disp%show("getLocNonSpace(SK_'  ')")
    call disp%show( getLocNonSpace(SK_'  ') )
    call disp%skip()

    call disp%skip()
    call disp%show("getLocNonSpace(SK_'a ')")
    call disp%show( getLocNonSpace(SK_'a ') )
    call disp%skip()

    call disp%skip()
    call disp%show("getLocNonSpace(SK_'aa')")
    call disp%show( getLocNonSpace(SK_'aa') )
    call disp%skip()

    call disp%skip()
    call disp%show("getLocNonSpace([character(2,SK) :: 'a ', ' a', 'aa'])")
    call disp%show( getLocNonSpace([character(2,SK) :: 'a ', ' a', 'aa']) )
    call disp%skip()

end program example