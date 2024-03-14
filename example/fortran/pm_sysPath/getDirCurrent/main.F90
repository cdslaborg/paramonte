program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysPath, only: getDirCurrent, isDir

    implicit none

    character(:, SK), allocatable   :: cpath
    character(2047, SK)             :: errmsg
    logical(LK)                     :: failed

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("cpath = getDirCurrent() ! get the current working directory.")
                    cpath = getDirCurrent()
    call disp%show("cpath")
    call disp%show( cpath , deliml = SK_"""" )
    call disp%show("isDir(cpath)")
    call disp%show( isDir(cpath) )
    call disp%skip()

    call disp%skip()
    call disp%show("cpath = getDirCurrent(failed) ! Gracefully capture any runtime error while retrieving the current working directory.")
                    cpath = getDirCurrent(failed)
    call disp%show("failed")
    call disp%show( failed )
    call disp%show("if (.not. failed) call disp%show(cpath)")
                    if (.not. failed) call disp%show(cpath, deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("cpath = getDirCurrent(failed, errmsg) ! Gracefully capture any runtime error while retrieving the current working directory.")
                    cpath = getDirCurrent(failed, errmsg)
    call disp%show("failed")
    call disp%show( failed )
    call disp%show("if (failed) then; call disp%show(trim(errmsg), deliml = SK_''''); else; call disp%show(cpath, deliml = SK_''''); end if )")
                    if (failed) then; call disp%show(trim(errmsg), deliml = SK_''''); else; call disp%show(cpath, deliml = SK_''''); end if
    call disp%skip()

end program example