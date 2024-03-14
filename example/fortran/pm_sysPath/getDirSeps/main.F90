program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! all processor types and kinds are supported.
    use pm_io, only: display_type
    use pm_sysPath, only: getDirSeps

    implicit none

    character(255, SK)  :: errmsg = SK_""
    logical(LK)         :: failed

    type(display_type)  :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getDirSeps()")
    call disp%show( getDirSeps() , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDirSeps(mold = SK_'a')")
    call disp%show( getDirSeps(mold = SK_'a') , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDirSeps(failed)")
    call disp%show( getDirSeps(failed) , deliml = """" )
    call disp%show("failed")
    call disp%show( failed )
    call disp%skip()

    call disp%skip()
    call disp%show("getDirSeps(failed, errmsg)")
    call disp%show( getDirSeps(failed, errmsg) , deliml = """" )
    call disp%show("failed")
    call disp%show( failed )
    call disp%show("trim(errmsg)")
    call disp%show( trim(errmsg) , deliml = SK_"""" )
    call disp%skip()

end program example