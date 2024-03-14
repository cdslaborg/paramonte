program example

    use pm_kind, only: IK, LK
    use pm_kind, only: SK ! All kinds are supported.
    use pm_err, only: setWarned
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
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Warn the user about an exceptional event.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("str")
    call disp%show( str, deliml = SK_"""" )
    call disp%show("call setWarned(str, newline = SK_'\n', unit = disp%unit)")
                    call setWarned(str, newline = SK_'\n', unit = disp%unit)
    call disp%skip()
    call disp%show("call setWarned(str, prefix = SK_'    ParaMonte', newline = SK_'\n', unit = disp%unit)")
                    call setWarned(str, prefix = SK_'    ParaMonte', newline = SK_'\n', unit = disp%unit)
    call disp%skip()
    call disp%show("call setWarned(str, prefix = SK_'    ParaMonte', newline = SK_'\n', width = 72_IK, unit = disp%unit)")
                    call setWarned(str, prefix = SK_'    ParaMonte', newline = SK_'\n', width = 72_IK, unit = disp%unit)
    call disp%skip()

end program example