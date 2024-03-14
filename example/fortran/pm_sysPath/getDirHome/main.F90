program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysPath, only: getDirHome, isDir

    implicit none

    character(:, SK), allocatable   :: hpath
    character(2047, SK)             :: errmsg
    logical(LK)                     :: failed

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("hpath = getDirHome() ! get the current user home directory.")
                    hpath = getDirHome()
    call disp%show("hpath")
    call disp%show( hpath , deliml = SK_"""" )
    call disp%show("isDir(hpath)")
    call disp%show( isDir(hpath) )
    call disp%skip()

    call disp%skip()
    call disp%show("hpath = getDirHome(user = SK_'') ! empty user returns the current user home.")
                    hpath = getDirHome(user = SK_'')
    call disp%show("hpath")
    call disp%show( hpath , deliml = SK_"""" )
    call disp%show("isDir(hpath)")
    call disp%show( isDir(hpath) )
    call disp%skip()

    call disp%skip()
    call disp%show("hpath = getDirHome(user = SK_'paramonte') ! get the home directory of user `paramonte`.")
                    hpath = getDirHome(user = SK_'paramonte')
    call disp%show("hpath")
    call disp%show( hpath , deliml = SK_"""" )
    call disp%show("isDir(hpath)")
    call disp%show( isDir(hpath) )
    call disp%skip()

    call disp%skip()
    call disp%show("hpath = getDirHome(failed = failed) ! Gracefully capture any runtime error while retrieving the home directory.")
                    hpath = getDirHome(failed = failed)
    call disp%show("failed")
    call disp%show( failed )
    call disp%show("if (.not. failed) call disp%show(hpath)")
                    if (.not. failed) call disp%show(hpath, deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("hpath = getDirHome(failed = failed, errmsg = errmsg) ! Gracefully capture any runtime error while retrieving the home directory.")
                    hpath = getDirHome(failed = failed, errmsg = errmsg)
    call disp%show("failed")
    call disp%show( failed )
    call disp%show("if (failed) then; call disp%show(trim(errmsg), deliml = SK_''''); else; call disp%show(hpath, deliml = SK_''''); end if )")
                    if (failed) then; call disp%show(trim(errmsg), deliml = SK_''''); else; call disp%show(hpath, deliml = SK_''''); end if
    call disp%skip()

end program example