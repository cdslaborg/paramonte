program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: isStrAlpha

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isStrAlpha(SK_'A')")
    call disp%show( isStrAlpha(SK_'A') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlpha(SK_'a')")
    call disp%show( isStrAlpha(SK_'a') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlpha(SK_'aBcD')")
    call disp%show( isStrAlpha(SK_'aBcD') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlpha(SK_'ParaMonte')")
    call disp%show( isStrAlpha(SK_'ParaMonte') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlpha(SK_'-a')")
    call disp%show( isStrAlpha(SK_'-a') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlpha(SK_'_A')")
    call disp%show( isStrAlpha(SK_'_A') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlpha(SK_'_')")
    call disp%show( isStrAlpha(SK_'_') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlpha(SK_' ')")
    call disp%show( isStrAlpha(SK_' ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrAlpha(SK_'')")
    call disp%show( isStrAlpha(SK_'') )
    call disp%skip()

end program example