program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! all processor types and kinds are supported.
    use pm_io, only: display_type
    use pm_sysPath, only: getPathSep

    implicit none

    character(255, SK)  :: errmsg = SK_""
    logical(LK)         :: failed

    type(display_type)  :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getPathSep()")
    call disp%show( getPathSep() , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("getPathSep(mold = SK_'a')")
    call disp%show( getPathSep(mold = SK_'a') , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("getPathSep(failed)")
    call disp%show( getPathSep(failed) , deliml = """" )
    call disp%show("failed")
    call disp%show( failed )
    call disp%skip()

    call disp%skip()
    call disp%show("getPathSep(failed, errmsg)")
    call disp%show( getPathSep(failed, errmsg) , deliml = """" )
    call disp%show("failed")
    call disp%show( failed )
    call disp%show("trim(errmsg)")
    call disp%show( trim(errmsg) , deliml = SK_"""" )
    call disp%skip()

end program example