program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_io, only: openArg_type

    implicit none

    integer(IK)                     :: iostat
    type(openArg_type)              :: openArg

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("openArg = openArg_type(file = 'foo.bar', status = 'new')")
                    openArg = openArg_type(file = 'foo.bar', status = 'new')
    call disp%skip

    call disp%skip
    call disp%show("openArg = openArg_type(file = 'foo.bar', status = 'scratch', iostat = iostat)")
                    openArg = openArg_type(file = 'foo.bar', status = 'scratch', iostat = iostat)
    call disp%show("iostat")
    call disp%show( iostat )
    call disp%show("openArg%iomsg")
    call disp%show( openArg%iomsg , deliml = SK_"""" )
    call disp%skip

end program example