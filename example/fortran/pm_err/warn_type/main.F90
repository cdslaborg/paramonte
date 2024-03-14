program example

    use pm_kind, only: IK, LK
    use pm_io, only: display_type
    use pm_err, only: warn_type
    use pm_kind, only: SK ! All kinds are supported.
    use iso_fortran_env, only: output_unit

    implicit none

    character(:, SK), allocatable :: str, strWrapped

    type(display_type) :: disp

    disp = display_type(file = SK_"main.out.F90")

    str =   "ParaMonte is a serial/parallel library of Monte Carlo routines for \n&
            &   +   optimization, \n&
            &   +   sampling, and \n&
            &   +   integration \n&
                    &of mathematical objective functions of arbitrary-dimensions, &
                    &in particular, the posterior distributions of Bayesian models in \n&
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
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Notify the user about an important message.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        type(warn_type) :: warn
        call disp%show("warn = warn_type(unit = disp%unit, newline = SK_'\n', prefix = SK_'    ParaMonte')")
                        warn = warn_type(unit = disp%unit, newline = SK_'\n', prefix = SK_'    ParaMonte')
        call disp%show("str")
        call disp%show( str, deliml = SK_"""" )
        call disp%skip()
        call disp%show("call warn%show(str)")
                        call warn%show(str)
        call disp%skip()
        call disp%show("call warn%show(str, width = 72_IK, sticky = .true._LK)")
                        call warn%show(str, width = 72_IK, sticky = .true._LK)
        call disp%skip()
        call disp%show("call warn%show(str) ! same width as before.")
                        call warn%show(str) ! same width as before.
        call disp%skip()
    end block

end program example