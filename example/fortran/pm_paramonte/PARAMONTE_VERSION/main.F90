program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_paramonte, only: PARAMONTE_VERSION

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("PARAMONTE_VERSION")
    call disp%show( PARAMONTE_VERSION , deliml = SK_"""" )
    call disp%skip()

end program example