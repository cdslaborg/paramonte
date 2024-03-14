program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: setStrQuoted

    implicit none

    character(:), allocatable :: strQuoted
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("call setStrQuoted(strQuoted, '')")
                    call setStrQuoted(strQuoted, '')
    call disp%show("strQuoted")
    call disp%show( strQuoted )
    call disp%skip()

    call disp%skip()
    call disp%show("call setStrQuoted(strQuoted, ' ')")
                    call setStrQuoted(strQuoted, ' ')
    call disp%show("strQuoted")
    call disp%show( strQuoted )
    call disp%skip()

    call disp%skip()
    call disp%show("call setStrQuoted(strQuoted, '"" ""')")
                    call setStrQuoted(strQuoted, '" "'  )
    call disp%show("strQuoted")
    call disp%show( strQuoted )
    call disp%skip()

    call disp%skip()
    call disp%show("call setStrQuoted(strQuoted, ''' ''')")
                    call setStrQuoted(strQuoted, ''' ''')
    call disp%show("strQuoted")
    call disp%show( strQuoted )
    call disp%skip()

    call disp%skip()
    call disp%show("call setStrQuoted(strQuoted, 'A')")
                    call setStrQuoted(strQuoted, 'A')
    call disp%show("strQuoted")
    call disp%show( strQuoted )
    call disp%skip()

    call disp%skip()
    call disp%show("call setStrQuoted(strQuoted, 'a""')")
                    call setStrQuoted(strQuoted, 'a""')
    call disp%show("strQuoted")
    call disp%show( strQuoted )
    call disp%skip()

    call disp%skip()
    call disp%show("call setStrQuoted(strQuoted, 'a""')")
                    call setStrQuoted(strQuoted, 'a""')
    call disp%show("strQuoted")
    call disp%show( strQuoted )
    call disp%skip()

end program example