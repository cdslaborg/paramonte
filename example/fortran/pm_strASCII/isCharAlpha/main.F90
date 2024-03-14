program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: isCharAlpha

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isCharAlpha('A')")
    call disp%show( isCharAlpha('A') )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharAlpha('a')")
    call disp%show( isCharAlpha('a') )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharAlpha(['.','B','c','D'])")
    call disp%show( isCharAlpha(['.','B','c','D']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharAlpha('-')")
    call disp%show( isCharAlpha('-') )
    call disp%skip()

    call disp%skip()
    call disp%show("isCharAlpha(' ')")
    call disp%show( isCharAlpha(' ') )
    call disp%skip()

end program example