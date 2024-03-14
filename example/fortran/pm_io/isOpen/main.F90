program example

    use iso_fortran_env, only: output_unit, input_unit, error_unit
    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_io, only: isOpen

    implicit none

    integer(IK) :: unit
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("isOpen('temp.tmp')")
    call disp%show( isOpen('temp.tmp') )
    call disp%skip

    call disp%skip
    call disp%show("open(newunit = unit, file = 'temp.tmp')")
                    open(newunit = unit, file = 'temp.tmp')
    call disp%show("isOpen('temp.tmp')")
    call disp%show( isOpen('temp.tmp') )
    call disp%show("isOpen(unit)")
    call disp%show( isOpen(unit) )
    call disp%skip

    call disp%skip
    call disp%show("isOpen([output_unit, input_unit, error_unit])")
    call disp%show( isOpen([output_unit, input_unit, error_unit]) )
    call disp%skip

end program example