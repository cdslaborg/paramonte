program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: isCharAlphaNum

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isCharAlphaNum('1')")
    call disp%show( isCharAlphaNum('1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharAlphaNum(['a','-',' ','3','A'])")
    call disp%show( isCharAlphaNum(['a','-',' ','3','A']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharAlphaNum('.')")
    call disp%show( isCharAlphaNum('.') )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharAlphaNum('+')")
    call disp%show( isCharAlphaNum('+') )
    call disp%skip()

end program example