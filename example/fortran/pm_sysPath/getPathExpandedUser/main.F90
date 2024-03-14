program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysPath, only: getPathExpandedUser
    use pm_sysPath, only: getDirHome

    implicit none

    character(:, SK), allocatable   :: epath
    character(2047, SK)             :: errmsg
    logical(LK)                     :: failed

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("epath = getPathExpandedUser(SK_'~')")
                    epath = getPathExpandedUser(SK_'~')
    call disp%show("epath")
    call disp%show( epath , deliml = SK_"""" )
    call disp%show("getDirHome()")
    call disp%show( getDirHome() , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("epath = getPathExpandedUser(SK_'~paramonte')")
                    epath = getPathExpandedUser(SK_'~paramonte')
    call disp%show("epath")
    call disp%show( epath , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("epath = getPathExpandedUser(SK_'~/paramonte')")
                    epath = getPathExpandedUser(SK_'~/paramonte')
    call disp%show("epath")
    call disp%show( epath , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("epath = getPathExpandedUser(SK_'./~/paramonte', failed = failed) ! Gracefully capture any runtime error while retrieving the current working directory.")
                    epath = getPathExpandedUser(SK_'./~/paramonte', failed = failed)
    call disp%show("failed")
    call disp%show( failed )
    call disp%show("if (.not. failed) call disp%show(epath)")
                    if (.not. failed) call disp%show(epath, deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("epath = getPathExpandedUser(SK_'', failed = failed) ! Gracefully capture any runtime error while retrieving the current working directory.")
                    epath = getPathExpandedUser(SK_'', failed = failed)
    call disp%show("failed")
    call disp%show( failed )
    call disp%show("if (.not. failed) call disp%show(epath)")
                    if (.not. failed) call disp%show(epath, deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("epath = getPathExpandedUser(SK_'~~/paramonte', failed = failed, errmsg = errmsg) ! Gracefully capture any runtime error while retrieving the current working directory.")
                    epath = getPathExpandedUser(SK_'~~/paramonte', failed = failed, errmsg = errmsg)
    call disp%show("failed")
    call disp%show( failed )
    call disp%show("if (failed) then; call disp%show(trim(errmsg), deliml = SK_''''); else; call disp%show(epath, deliml = SK_''''); end if )")
                    if (failed) then; call disp%show(trim(errmsg), deliml = SK_''''); else; call disp%show(epath, deliml = SK_''''); end if
    call disp%skip()

end program example