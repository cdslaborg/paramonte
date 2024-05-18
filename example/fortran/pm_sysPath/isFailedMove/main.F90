program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysPath, only: isFailedMove, isFile, isDir, getPathNew, getDirSep

    implicit none

    character(:, SK), allocatable :: from, to

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Move file to a new non-existing destination file.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("from = SK_'main.F90' ! define the source path name.")
                    from = SK_'main.F90'
    call disp%show("to = getPathNew(ext = SK_'.main.F90') ! define a random unique path name.")
                    to = getPathNew(ext = SK_'.main.F90')
    call disp%show("to")
    call disp%show( to , deliml = SK_"""" )
    call disp%show("[isFile(from), isFile(to)]")
    call disp%show( [isFile(from), isFile(to)] )
    call disp%show("if (isFailedMove(from, to)) error stop 'movement failed.'")
                    if (isFailedMove(from, to)) error stop 'movement failed.'
    call disp%show("[isFile(from), isFile(to)]")
    call disp%show( [isFile(from), isFile(to)] )
    call disp%show("if (isFailedMove(to, from)) error stop 'movement failed.'")
                    if (isFailedMove(to, from)) error stop 'movement failed.'
    call disp%show("[isFile(from), isFile(to)]")
    call disp%show( [isFile(from), isFile(to)] )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Move file to a new non-existing doubly-nested destination folder.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("from = SK_'main.F90' ! define the source path name.")
                    from = SK_'main.F90'
    call disp%show("to = getPathNew(dir = getPathNew(), prefix = SK_'sub')//getDirSep() ! Ending the directory name with directory separator is crucial if it does not exit already, otherwise it is interpreted as file name.")
                    to = getPathNew(dir = getPathNew(), prefix = SK_'sub')//getDirSep()
    call disp%show("to")
    call disp%show( to , deliml = SK_"""" )
    call disp%show("[isDir(from), isFile(from), isDir(to), isFile(to)]")
    call disp%show( [isDir(from), isFile(from), isDir(to), isFile(to)] )
    call disp%show("if (isFailedMove(from, to)) error stop 'movement failed.'")
                    if (isFailedMove(from, to)) error stop 'movement failed.'
    call disp%show("[isFile(from), isFile(to//getDirSep()//from)]")
    call disp%show( [isFile(from), isFile(to//getDirSep()//from)] )
    call disp%show("if (isFailedMove(to//from, from)) error stop 'movement failed.'")
                    if (isFailedMove(to//from, from)) error stop 'movement failed.'
    call disp%show("[isFile(from), isFile(to)]")
    call disp%show( [isFile(from), isFile(to)] )
    call disp%skip()

end program example