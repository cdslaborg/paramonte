program example

    use pm_kind, only: IK, LK
    use pm_io, only: display_type
    use pm_err, only: note_type
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
        type(note_type) :: note
        call disp%show("note = note_type(unit = disp%unit, newline = SK_'\n', prefix = SK_'    ParaMonte')")
                        note = note_type(unit = disp%unit, newline = SK_'\n', prefix = SK_'    ParaMonte')
        call disp%show("str")
        call disp%show( str, deliml = SK_"""" )
        call disp%skip()
        call disp%show("call note%show(str)")
                        call note%show(str)
        call disp%skip()
        call disp%show("call note%show(str, width = 72_IK, sticky = .true._LK)")
                        call note%show(str, width = 72_IK, sticky = .true._LK)
        call disp%skip()
        call disp%show("call note%show(str) ! same width as before.")
                        call note%show(str) ! same width as before.
        call disp%skip()
    end block

end program example