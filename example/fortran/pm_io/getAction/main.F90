program example

    use iso_fortran_env, only: output_unit, input_unit, error_unit
    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_sysPath, only: getPathNew
    use pm_io, only: getAction

    implicit none

    integer(IK) :: unit
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("getAction('temp.tmp') ! file is not connected.")
    call disp%show( getAction('temp.tmp') , deliml = SK_"""" )
    call disp%skip

    call disp%skip
    call disp%show("open(newunit = unit, file = 'temp.tmp', action = 'write')")
                    open(newunit = unit, file = 'temp.tmp', action = 'write')
    call disp%show("getAction('temp.tmp')")
    call disp%show( getAction('temp.tmp') , deliml = SK_"""" )
    call disp%show("getAction(unit)")
    call disp%show( getAction(unit) , deliml = SK_"""" )
    call disp%show("close(unit)")
                    close(unit)
    call disp%skip

    call disp%skip
    call disp%show("open(newunit = unit, file = 'temp.tmp', action = 'read')")
                    open(newunit = unit, file = 'temp.tmp', action = 'read')
    call disp%show("getAction('temp.tmp')")
    call disp%show( getAction('temp.tmp') , deliml = SK_"""" )
    call disp%show("getAction(unit)")
    call disp%show( getAction(unit) , deliml = SK_"""" )
    call disp%show("close(unit)")
                    close(unit)
    call disp%skip

    call disp%skip
    call disp%show("open(newunit = unit, file = 'temp.tmp', action = 'readwrite')")
                    open(newunit = unit, file = 'temp.tmp', action = 'readwrite')
    call disp%show("getAction('temp.tmp')")
    call disp%show( getAction('temp.tmp') , deliml = SK_"""" )
    call disp%show("getAction(unit)")
    call disp%show( getAction(unit) , deliml = SK_"""" )
    call disp%show("close(unit)")
                    close(unit)
    call disp%skip

    call disp%skip
    call disp%show("getAction([output_unit, input_unit, error_unit])")
    call disp%show( getAction([output_unit, input_unit, error_unit]) , deliml = SK_"""" )
    call disp%skip

end program example