program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysPath, only: isFailedMakeDirTemp, isDir

    implicit none

    character(:, SK), allocatable :: tempdir

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isFailedMakeDirTemp(tempdir)")
    call disp%show( isFailedMakeDirTemp(tempdir) )
    call disp%show("tempdir")
    call disp%show( tempdir , deliml = SK_"""" )
    call disp%show("isDir(tempdir)")
    call disp%show( isDir(tempdir) )
    call disp%skip()

end program example