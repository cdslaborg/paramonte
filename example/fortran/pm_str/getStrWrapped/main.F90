program example

    use pm_strASCII, only: LF ! linefeed for better display
    use pm_kind, only: IK, LK
    use pm_kind, only: SK ! All kinds are supported.
    use pm_str, only: getStrWrapped
    use pm_io, only: display_type
    use pm_arrayReplace, only: getReplaced

    implicit none

    character(:, SK), allocatable   :: str, strWrapped

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    str =   "A man asked my friend Jaime Cohen: `What is the human being's funniest characteristic?` &
            &Cohen said: ‘Our contradictoriness. We are in such a hurry to grow up, and then we long for our lost childhood. &
            &We make ourselves ill earning money, and then spend all our money on getting well again. &
            &We think so much about the future that we neglect the present, and thus experience neither the present nor the future. &
            &We live as if we were never going to die, and die as if we had never lived."

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Wrap string within the default 132 characters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("str")
    call disp%show( str, deliml = SK_"""" )
    call disp%show("strWrapped = getStrWrapped(str)")
                    strWrapped = getStrWrapped(str)
    call disp%show("strWrapped")
    call disp%show( strWrapped, deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Wrap string with an indent.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("'    '//str")
    call disp%show( '    '//str, deliml = SK_"""" )
    call disp%show("strWrapped = getStrWrapped('    '//str)")
                    strWrapped = getStrWrapped('    '//str)
    call disp%show("LF//strWrapped")
    call disp%show( LF//strWrapped, deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Wrap string with an indent and a prefix for a given non-default width.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("'....'//str")
    call disp%show( '....'//str, deliml = SK_"""" )
    call disp%show("strWrapped = getStrWrapped('....'//str, prefix = SK_'NOTE - ', indent = SK_'.', width = 100_IK)")
                    strWrapped = getStrWrapped('....'//str, prefix = SK_'NOTE - ', indent = SK_'.', width = 100_IK)
    call disp%show("LF//strWrapped")
    call disp%show( LF//strWrapped, deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Wrap string with an indent and a prefix for a given non-default width at any character by specifying an empty `break` to ensure a precise width.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("'    '//str")
    call disp%show( '    '//str, deliml = SK_"""" )
    call disp%show("strWrapped = getStrWrapped('    '//str, prefix = SK_'NOTE - ', indent = SK_'    ', break = SK_'', width = 50_IK)")
                    strWrapped = getStrWrapped('    '//str, prefix = SK_'NOTE - ', indent = SK_'    ', break = SK_'', width = 50_IK)
    call disp%show("LF//strWrapped")
    call disp%show( LF//strWrapped, deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Setting a too small `width` (e.g., less than the initial indentation)can lead to awkward wrapping.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("'....'//str")
    call disp%show( '....'//str, deliml = SK_"""" )
    call disp%show("strWrapped = getStrWrapped(repeat(SK_'.', 8)//str, prefix = SK_'NOTE - ', indent = repeat(SK_'.', 4), width = 5_IK)")
                    strWrapped = getStrWrapped(repeat(SK_'.', 8)//str, prefix = SK_'NOTE - ', indent = repeat(SK_'.', 4), width = 5_IK)
    call disp%show("LF//strWrapped")
    call disp%show( LF//strWrapped, deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Wrap string while converting all C-style newline instances (\n) to actual linefeed characters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    str =   "The ParaMonte serial/parallel library contains Monte Carlo routines for \n\n&
            &    +   optimization, \n&
                &+   sampling, and \n&
                &+   integration \n\n&
            &of mathematical objective functions of arbitrary-dimensions in particular, the posterior distributions of Bayesian models in \n\n&
            &    +   data science, \n&
                &+   Machine Learning, and \n&
                &+   scientific inference, \n\n&
            &with the design goal of unifying the \n&
            &    +   automation (of Monte Carlo simulations), \n&
            &    +   user-friendliness (of the library), \n&
            &    +   accessibility (from multiple programming environments), \n&
            &    +   high-performance (at runtime), and \n&
            &    +   scalability (across many parallel processors)."

    call disp%show("str")
    call disp%show( str, deliml = SK_"""" )
    call disp%show("strWrapped = getStrWrapped(str, prefix = SK_'NOTE - ', newline = '\n', linefeed = LF, width = 80_IK)")
                    strWrapped = getStrWrapped(str, prefix = SK_'NOTE - ', newline = '\n', linefeed = LF, width = 80_IK)
    call disp%show("LF//strWrapped")
    call disp%show( LF//strWrapped, deliml = SK_"""" )
    call disp%skip()
    call disp%show("strWrapped = getReplaced(getStrWrapped(str, prefix = SK_'NOTE - ', newline = '\n', width = 80_IK), pattern = '\n', replacement = LF) ! The slower verbose way of achieving the above.")
                    strWrapped = getReplaced(getStrWrapped(str, prefix = SK_'NOTE - ', newline = '\n', width = 80_IK), pattern = '\n', replacement = LF)
    call disp%show("LF//strWrapped")
    call disp%show( LF//strWrapped, deliml = SK_"""" )
    call disp%skip()

    str = SK_'    '//str

    call disp%show("str")
    call disp%show( str, deliml = SK_"""" )
    call disp%show("strWrapped = getStrWrapped(str, prefix = SK_'NOTE - ', newline = '\n', linefeed = LF, width = 80_IK)")
                    strWrapped = getStrWrapped(str, prefix = SK_'NOTE - ', newline = '\n', linefeed = LF, width = 80_IK)
    call disp%show("LF//strWrapped")
    call disp%show( LF//strWrapped, deliml = SK_"""" )
    call disp%skip()
    call disp%show("strWrapped = getReplaced(getStrWrapped(str, prefix = SK_'NOTE - ', newline = '\n', width = 80_IK), pattern = '\n', replacement = LF) ! The slower verbose way of achieving the above.")
                    strWrapped = getReplaced(getStrWrapped(str, prefix = SK_'NOTE - ', newline = '\n', width = 80_IK), pattern = '\n', replacement = LF)
    call disp%show("LF//strWrapped")
    call disp%show( LF//strWrapped, deliml = SK_"""" )
    call disp%skip()

end program example