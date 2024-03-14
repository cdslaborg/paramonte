program example

    use pm_kind, only: IK, LK
    use pm_kind, only: SK ! All kinds are supported.
    use pm_err, only: setAborted
    use pm_io, only: display_type

    implicit none

    character(:, SK), allocatable   :: str, strWrapped

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    str =   "ParaMonte is a serial/parallel library of Monte Carlo routines for \n&
            &   +   optimization, \n&
            &   +   sampling, and \n&
            &   +   integration \n&
            &of mathematical objective functions of arbitrary-dimensions in particular, the posterior distributions of Bayesian models in \n&
            &   +   data science, \n&
            &   +   Machine Learning, and \n&
            &   +   scientific inference, \n&
            &with the design goal of unifying the \n&
            &   +   automation (of Monte Carlo simulations), \n&
            &   +   user-friendliness (of the library), \n&
            &   +   accessibility (from multiple programming environments), \n&
            &   +   high-performance (at runtime), and \n&
            &   +   scalability (across many parallel processors)."

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Notify the user about a fatal error occurrence.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("str")
    call disp%show( str, deliml = SK_"""" )
    call disp%show("call setAborted(str, newline = SK_'\n', renabled = .true._LK, unit = disp%unit) ! To halt the program, omit `renabled`.")
                    call setAborted(str, newline = SK_'\n', renabled = .true._LK, unit = disp%unit)
    call disp%skip()
    call disp%show("call setAborted(str, prefix = SK_'    ParaMonte', newline = SK_'\n', renabled = .true._LK, unit = disp%unit)")
                    call setAborted(str, prefix = SK_'    ParaMonte', newline = SK_'\n', renabled = .true._LK, unit = disp%unit)
    call disp%skip()
    call disp%show("call setAborted(str, prefix = SK_'    ParaMonte', newline = SK_'\n', width = 72_IK, renabled = .true._LK, unit = disp%unit)")
                    call setAborted(str, prefix = SK_'    ParaMonte', newline = SK_'\n', width = 72_IK, renabled = .true._LK, unit = disp%unit)
    call disp%skip()
    call disp%show("call setAborted(SK_'    '//str, help = SK_'Contact cdslab.org for help.', prefix = SK_'    ParaMonte', newline = SK_'\n', width = 72_IK, renabled = .true._LK, unit = disp%unit)")
                    call setAborted(SK_'    '//str, help = SK_'Contact cdslab.org for help.', prefix = SK_'    ParaMonte', newline = SK_'\n', width = 72_IK, renabled = .true._LK, unit = disp%unit)
    call disp%skip()

end program example