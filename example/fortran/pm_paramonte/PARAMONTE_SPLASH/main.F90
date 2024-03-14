program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_paramonte, only: PARAMONTE_SPLASH

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("PARAMONTE_SPLASH")
    call disp%show( PARAMONTE_SPLASH , deliml = SK_"""" )
    call disp%skip()

end program example