program example

    use pm_kind, only: LK, IK
    use pm_kind, only: SK ! all processor types and kinds are supported.
    use pm_io, only: display_type
    use pm_sysPath, only: getPathTemp

    implicit none

    logical(LK) :: failed

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getPathTemp()")
    call disp%show( getPathTemp() , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getPathTemp(prefix = SK_'paramonte')")
    call disp%show( getPathTemp(prefix = SK_'paramonte') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getPathTemp(prefix = SK_'paramonte', ext = SK_'.txt')")
    call disp%show( getPathTemp(prefix = SK_'paramonte', ext = SK_'.txt') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getPathTemp(prefix = SK_'paramonte', ext = SK_'.txt', pid = 3_IK)")
    call disp%show( getPathTemp(prefix = SK_'paramonte', ext = SK_'.txt', pid = 3_IK) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getPathTemp(prefix = SK_'paramonte', pid = 3_IK, failed = failed)")
    call disp%show( getPathTemp(prefix = SK_'paramonte', pid = 3_IK, failed = failed) , deliml = SK_"""" )
    call disp%show("failed")
    call disp%show( failed )
    call disp%skip()

    call disp%skip()
    call disp%show("getPathTemp(prefix = SK_'paramonte', sep = SK_'.', pid = 3_IK, failed = failed)")
    call disp%show( getPathTemp(prefix = SK_'paramonte', sep = SK_'.', pid = 3_IK, failed = failed) , deliml = SK_"""" )
    call disp%show("failed")
    call disp%show( failed )
    call disp%skip()

end program example