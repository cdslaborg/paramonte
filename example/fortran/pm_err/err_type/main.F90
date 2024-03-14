program example

    use pm_kind, only: IK, LK
    use pm_kind, only: SK ! All kinds are supported.
    use pm_err, only: err_type
    use pm_io, only: display_type
    use iso_fortran_env, only: output_unit

    implicit none

    integer :: fileunit
    type(err_type) :: err
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%show("err = err_type()")
                    err = err_type()
    call disp%show("err%msg")
    call disp%show( err%msg , deliml = SK_"""" )
    call disp%show("err%stat")
    call disp%show( err%stat )
    call disp%show("open(newunit = fileunit, file = 'nonesense.io', status = 'old', iostat = err%stat, iomsg = err%msg)")
                    open(newunit = fileunit, file = 'nonesense.io', status = 'old', iostat = err%stat, iomsg = err%msg)
    call disp%skip()
    call disp%show("err%stat")
    call disp%show( err%stat )
    call disp%skip()
    call disp%show("trim(err%msg)")
    call disp%show( trim(err%msg) )
    call disp%skip()

end program example