program example

    use pm_kind, only: IK, LK
    use pm_kind, only: SK ! All kinds are supported.
    use pm_str, only: getLenIndent
    use pm_io, only: display_type
    use pm_strASCII, only: LF ! linefeed for better display

    implicit none

    character(:, SK), allocatable   :: str

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    str = ".-+.-+ParaMonte"

    call disp%show("str")
    call disp%show( str, deliml = SK_"""" )
    call disp%show("getLenIndent(str, SK_'.-+')")
    call disp%show( getLenIndent(str, SK_'.-+') )
    call disp%skip()

    str = "ParaMonte"

    call disp%show("str")
    call disp%show( str, deliml = SK_"""" )
    call disp%show("getLenIndent(str, SK_'.-+')")
    call disp%show( getLenIndent(str, SK_'.-+') )
    call disp%skip()

    call disp%show("getLenIndent([character(50,SK) :: 'ParaMonte', ' ParaMonte', '  ParaMonte'], SK_'  ')")
    call disp%show( getLenIndent([character(50,SK) :: 'ParaMonte', ' ParaMonte', '  ParaMonte'], SK_'  ') )
    call disp%skip()

end program example