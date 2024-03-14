program example

    use pm_io, only: wrap_type
    use pm_io, only: display_type
    use pm_str, only: getStrWrapped
    use pm_kind, only: SK, IK, LK ! All kinds are supported.
    use pm_strASCII, only: LF ! linefeed for better display.

    implicit none

    character(:, SK), allocatable   :: str, strWrapped

    type(wrap_type) :: text
    type(display_type) :: disp

    disp = display_type(file = SK_"main.out.F90")

    str =   "ParaMonte is a serial/parallel library of Monte Carlo routines for sampling mathematical density functions of arbitrary-dimensions&
            &and Machine Learning algorithms for scientific inference, with the design goal of unifying **automation** (of simulations and tasks),&
            &**user-friendliness** (of algorithms), **accessibility** (from any platform or programming environment),&
            &**high-performance** (at runtime), and **scalability** (across many parallel processors)."

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Decorate the input 72-characters wrapped string within the default width.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("str")
    call disp%show( str, deliml = SK_"""" )
    call disp%show("strWrapped = getStrWrapped(str, width = 72_IK)")
                    strWrapped = getStrWrapped(str, width = 72_IK)
    call disp%show("text = wrap_type(unit = disp%unit)")
                    text = wrap_type(unit = disp%unit)
    call disp%show("call text%wrap(LF//strWrapped//LF)")
                    call text%wrap(LF//strWrapped//LF)
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Decorate the input 72-characters wrapped string within non-default settings.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("strWrapped")
    call disp%show( strWrapped, deliml = SK_"""" )
    call disp%show("call text%wrap(LF//strWrapped//LF, lwsize = 1_IK, rwsize = 1_IK, lwfill = SK_'|', rwfill = SK_'|', twfill = SK_'-', bwfill = SK_'-')")
                    call text%wrap(LF//strWrapped//LF, lwsize = 1_IK, rwsize = 1_IK, lwfill = SK_'|', rwfill = SK_'|', twfill = SK_'-', bwfill = SK_'-')
    call disp%skip()

    call disp%skip()
    call disp%show("strWrapped")
    call disp%show( strWrapped, deliml = SK_"""" )
    call disp%show("call text%wrap(LF//strWrapped//LF, fill = SK_'.', lwsize = 1_IK, rwsize = 1_IK, lwfill = SK_'\', rwfill = SK_'/', twfill = SK_'_', bwfill = SK_'_')")
                    call text%wrap(LF//strWrapped//LF, fill = SK_'.', lwsize = 1_IK, rwsize = 1_IK, lwfill = SK_'\', rwfill = SK_'/', twfill = SK_'_', bwfill = SK_'_')
    call disp%skip()

    str =   "ParaMonte is a serial/parallel library of Monte Carlo routines for sampling mathematical objective functions of arbitrary-dimensions."

    call disp%skip()
    call disp%show("str")
    call disp%show( str, deliml = SK_"""" )
    call disp%show("strWrapped = getStrWrapped(str)")
                    strWrapped = getStrWrapped(str)
    call disp%show("call text%wrap(strWrapped) ! Note what happens to the overflowed text.")
                    call text%wrap(strWrapped) ! Note what happens to the overflowed text.
    call disp%skip()

    call disp%skip()
    call disp%show("call text%wrap(LF//SK_'ParaMonte'//LF)")
                    call text%wrap(LF//SK_'ParaMonte'//LF)
    call disp%skip()

    call disp%skip()
    call disp%show("call text%wrap('ParaMonte')")
                    call text%wrap('ParaMonte')
    call disp%skip()

    call disp%skip()
    call disp%show("call text%wrap(LF//SK_''//LF)")
                    call text%wrap(LF//SK_''//LF)
    call disp%skip()

    call disp%skip()
    call disp%show("call text%wrap(LF)")
                    call text%wrap(LF)
    call disp%skip()

    call disp%skip()
    call disp%show("call text%wrap('')")
                    call text%wrap('')
    call disp%skip()

end program example