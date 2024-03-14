program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use iso_fortran_env, only: input_unit, output_unit, error_unit
    use pm_io, only: isPreconnected

    implicit none

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("open(unit = 11, file = 'temp.tmp')")
                    open(unit = 11, file = 'temp.tmp')
    call disp%show("isPreconnected(11_IK)")
    call disp%show( isPreconnected(11_IK) )
    call disp%show("close(unit = 11)")
                    close(unit = 11)
    call disp%skip

    call disp%skip
    call disp%show("isPreconnected([integer(IK) :: input_unit, output_unit, error_unit])")
    call disp%show( isPreconnected([integer(IK) :: input_unit, output_unit, error_unit]) )
    call disp%skip

end program example