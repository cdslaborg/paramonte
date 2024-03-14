program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: getStrQuoted

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getStrQuoted('')")
    call disp%show( getStrQuoted('') )
    call disp%skip()

    call disp%skip()
    call disp%show("getStrQuoted(' ')")
    call disp%show( getStrQuoted(' ') )
    call disp%skip()

    call disp%skip()
    call disp%show("getStrQuoted('"" ""')")
    call disp%show( getStrQuoted('" "') )
    call disp%skip()

    call disp%skip()
    call disp%show("getStrQuoted(''' ''')")
    call disp%show( getStrQuoted(''' ''') )
    call disp%skip()

    call disp%skip()
    call disp%show("getStrQuoted('A')")
    call disp%show( getStrQuoted('A') )
    call disp%skip()

    call disp%skip()
    call disp%show("getStrQuoted('a""')")
    call disp%show( getStrQuoted('a"' ) )
    call disp%skip()

end program example