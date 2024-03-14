program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_container, only: css_type
    use pm_arrayResize, only: setResized
    use pm_sysPath, only: isFailedRemove, isFailedMakeDir, isExtant, isFile, getPathNew, getDirSep, ls, glob, isDir

    implicit none

    integer :: i, unit

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%")
    call disp%show("! Remove file.")
    call disp%show("!%%%%%%%%%%%%%")
    call disp%skip()

    block
        character(:, SK), allocatable :: path, dir
        call disp%skip()
        call disp%show("dir = getPathNew(prefix = SK_'tempdir_') ! define a temp directory name.")
                        dir = getPathNew(prefix = SK_'tempdir_')
        call disp%show("dir")
        call disp%show( dir , deliml = SK_"""" )
        call disp%show("if (isFailedMakeDir(dir)) error stop 'directory creation failed.'")
                        if (isFailedMakeDir(dir)) error stop 'directory creation failed.'
        call disp%show("path = getPathNew(dir) ! define a temp path name.")
                        path = getPathNew(dir)
        call disp%show("path")
        call disp%show( path , deliml = SK_"""" )
        call disp%show("open(newunit = unit, file = path, status = 'replace'); close(unit) ! Create a file in a randomly-named directory.")
                        open(newunit = unit, file = path, status = 'replace'); close(unit) ! Create a file in a randomly-named directory.
        call disp%show("isExtant(path)")
        call disp%show( isExtant(path) )
        call disp%show("if (isFailedRemove(path)) error stop 'removal failed.'")
                        if (isFailedRemove(path)) error stop 'removal failed.'
        call disp%show("[isFile(path), isExtant(path)]")
        call disp%show( [isFile(path), isExtant(path)] )
        call disp%show("if (.not. isFailedRemove(path)) error stop 'removal of non existing path fails unless `forced`.'")
                        if (.not. isFailedRemove(path)) error stop 'removal of non existing path fails unless `forced`.'
        call disp%show("if (isFailedRemove(path, forced = .true._LK)) error stop 'removal of non existing path fails unless `forced`.'")
                        if (isFailedRemove(path, forced = .true._LK)) error stop 'removal of non existing path fails unless `forced`.'
        call disp%show("isDir(dir)")
        call disp%show( isDir(dir) )
        call disp%show("if (isFailedRemove(dir, forced = .true._LK, recursive = .true._LK)) error stop 'removal of directory failed.'")
                        if (isFailedRemove(dir, forced = .true._LK, recursive = .true._LK)) error stop 'removal of directory failed.'
        call disp%show("isDir(dir)")
        call disp%show( isDir(dir) )
        call disp%skip()
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove directory or files matching pattern.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        type(css_type), allocatable :: list(:)
        character(:, SK), allocatable :: dir, pattern
        call disp%skip()
        call disp%show("dir = getPathNew(prefix = SK_'tempdir_') ! define a temp directory name.")
                        dir = getPathNew(prefix = SK_'tempdir_')
        call disp%show("dir")
        call disp%show( dir , deliml = SK_"""" )
        call disp%show("if (isFailedMakeDir(dir)) error stop 'directory creation failed.'")
                        if (isFailedMakeDir(dir)) error stop 'directory creation failed.'
        call disp%show("call setResized(list, 5_IK)")
                        call setResized(list, 5_IK)
        call disp%show("do i = 1, 3; list(i)%val = getPathNew(dir, prefix = SK_'tempfile_'); end do ! define temp path names.")
                        do i = 1, 3; list(i)%val = getPathNew(dir, prefix = SK_'tempfile_'); end do
        call disp%show("do i = 4, size(list); list(i)%val = getPathNew(dir, prefix = SK_'dumpfile_'); end do ! define temp path names.")
                        do i = 4, size(list); list(i)%val = getPathNew(dir, prefix = SK_'dumpfile_'); end do
        call disp%show("list")
        call disp%show( list , deliml = SK_"""" )
        call disp%show("ls(dir)")
        call disp%show( ls(dir) , deliml = SK_"""" )
        call disp%show("do i = 1, size(list); open(newunit = unit, file = list(i)%val, status = 'replace'); close(unit); end do ! Create random files with specified patterns.")
                        do i = 1, size(list); open(newunit = unit, file = list(i)%val, status = 'replace'); close(unit); end do ! Create random files with specified patterns.
        call disp%show("ls(dir)")
        call disp%show( ls(dir) , deliml = SK_"""" )
        call disp%show("pattern = dir//getDirSep()//SK_'*ump*'")
                        pattern = dir//getDirSep()//SK_'*ump*'
        call disp%show("pattern")
        call disp%show( pattern , deliml = SK_"""" )
        call disp%show("glob(pattern)")
        call disp%show( glob(pattern) , deliml = SK_"""" )
        call disp%show("if (isFailedRemove(pattern, forced = .true._LK, recursive = .true._LK)) error stop 'removal of file(s) matching pattern '''//pattern//''' failed.'")
                        if (isFailedRemove(pattern, forced = .true._LK, recursive = .true._LK)) error stop 'removal of file(s) matching pattern '''//pattern//''' failed.'
        call disp%show("glob(pattern)")
        call disp%show( glob(pattern) , deliml = SK_"""" )
        call disp%skip()
        call disp%show("pattern = SK_'./*emp*'")
                        pattern = SK_'./*emp*'
        call disp%show("pattern")
        call disp%show( pattern , deliml = SK_"""" )
        call disp%show("glob(pattern)")
        call disp%show( glob(pattern) , deliml = SK_"""" )
        call disp%show("if (isFailedRemove(pattern, forced = .true._LK, recursive = .true._LK)) error stop 'removal of directories and files matching pattern '''//pattern//''' failed.'")
                        if (isFailedRemove(pattern, forced = .true._LK, recursive = .true._LK)) error stop 'removal of directories and files matching pattern '''//pattern//''' failed.'
        call disp%show("glob(pattern)")
        call disp%show( glob(pattern) , deliml = SK_"""" )
        call disp%skip()
    end block

end program example