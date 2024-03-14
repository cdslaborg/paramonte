program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_paramonte, only: PARAMONTE_BUILD_DATE

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("PARAMONTE_BUILD_DATE")
    call disp%show( PARAMONTE_BUILD_DATE , deliml = SK_"""" )
    call disp%skip()

end program example