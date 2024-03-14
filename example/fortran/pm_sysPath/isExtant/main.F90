program example

    use iso_fortran_env, only: output_unit, input_unit, error_unit
    use pm_kind, only: SK, IK, LK
    use pm_sysPath, only: isExtant
    use pm_io, only: display_type

    implicit none

    integer(IK) :: unit
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("isExtant('temp.tmp')")
    call disp%show( isExtant('temp.tmp') )
    call disp%skip

    call disp%skip
    call disp%show("open(newunit = unit, file = 'temp.tmp')")
                    open(newunit = unit, file = 'temp.tmp')
    call disp%show("isExtant('temp.tmp')")
    call disp%show( isExtant('temp.tmp') )
    call disp%show("close(unit)")
                    close(unit)
    call disp%show("isExtant('temp.tmp')")
    call disp%show( isExtant('temp.tmp') )
    call disp%skip

    call disp%skip
    call disp%show("isExtant('main.F90')")
    call disp%show( isExtant('main.F90') )
    call disp%skip

    call disp%skip
    call disp%show("isExtant('.') ! Current directory.")
    call disp%show( isExtant('.') )
    call disp%skip

    call disp%skip
    call disp%show("isExtant('./') ! Current directory.")
    call disp%show( isExtant('./') )
    call disp%skip

    call disp%skip
    call disp%show("isExtant('..') ! Parent directory.")
    call disp%show( isExtant('..') )
    call disp%skip

    call disp%skip
    call disp%show("isExtant('../') ! Parent directory.")
    call disp%show( isExtant('../') )
    call disp%skip

end program example