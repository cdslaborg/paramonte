program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysPath, only: isFailedCopy, isFile, isDir, getPathNew, getDirSep

    implicit none

    character(:, SK), allocatable :: cpath, dir, subdir, destination

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Copy file to a new non-existing destination file.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("cpath = getPathNew(ext = SK_'.main.out.F90') ! define a random unique path name.")
                    cpath = getPathNew(ext = SK_'.main.out.F90')
    call disp%show("cpath")
    call disp%show( cpath , deliml = SK_"""" )
    call disp%show("isFile(cpath)")
    call disp%show( isFile(cpath) )
    call disp%show("isFailedCopy(from = SK_'main.out.F90', to = cpath)")
    call disp%show( isFailedCopy(from = SK_'main.out.F90', to = cpath) )
    call disp%show("isFile(cpath)")
    call disp%show( isFile(cpath) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Copy file to a new non-existing doubly-nested destination folder.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("dir = getPathNew(dir = getPathNew(), prefix = SK_'sub')//getDirSep() ! Ending the directory name with directory separator is crucial if it does not exit already, otherwise it is interpreted as file name.")
                    dir = getPathNew(dir = getPathNew(), prefix = SK_'sub')//getDirSep()
    call disp%show("dir")
    call disp%show( dir , deliml = SK_"""" )
    call disp%show("isDir(dir)")
    call disp%show( isDir(dir) )
    call disp%show("isFile(dir)")
    call disp%show( isFile(dir) )
    call disp%show("isFailedCopy(from = SK_'main.out.F90', to = dir)")
    call disp%show( isFailedCopy(from = SK_'main.out.F90', to = dir) )
    call disp%show("isDir(dir)")
    call disp%show( isDir(dir) )
    call disp%show("isFile(dir//SK_'main.out.F90')")
    call disp%show( isFile(dir//SK_'main.out.F90') )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Copy file to a new location and overwrite the existing file.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("destination = dir//SK_'main.out.F90'")
                    destination = dir//SK_'main.out.F90'
    call disp%show("destination")
    call disp%show( destination , deliml = SK_"""" )
    call disp%show("isFile(destination)")
    call disp%show( isFile(destination) )
    call disp%show("isFailedCopy(from = SK_'main.out.F90', to = destination, forced = .true.)")
    call disp%show( isFailedCopy(from = SK_'main.out.F90', to = destination, forced = .true.) )
    call disp%show("isDir(dir)")
    call disp%show( isDir(dir) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Copy file to a new non-existing destination file name in a non-existing subfolder.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("subdir = getPathNew(dir = dir, prefix = SK_'subsub')")
                    subdir = getPathNew(dir = dir, prefix = SK_'subsub')
    call disp%show("subdir")
    call disp%show( subdir , deliml = SK_"""" )
    call disp%show("isDir(subdir)")
    call disp%show( isDir(subdir) )
    call disp%show("cpath = subdir//getDirSep()//SK_'main.out.F90'")
                    cpath = subdir//getDirSep()//SK_'main.out.F90'
    call disp%show("cpath")
    call disp%show( cpath , deliml = SK_"""" )
    call disp%show("isFile(cpath)")
    call disp%show( isFile(cpath) )
    call disp%show("isFailedCopy(from = SK_'main.out.F90', to = cpath)")
    call disp%show( isFailedCopy(from = SK_'main.out.F90', to = cpath) )
    call disp%show("isFile(cpath)")
    call disp%show( isFile(cpath) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Copy 'the contents of' a folder recursively into another (potentially non-existing) destination folder.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("dir ! Note the ending directory separator which is essential for copying the contents of the folder and not the folder itself.")
    call disp%show( dir , deliml = SK_"""" )
    call disp%show("isDir(dir)")
    call disp%show( isDir(dir) )
    call disp%show("cpath = getPathNew(prefix = SK_'recursive_copy')")
                    cpath = getPathNew(prefix = SK_'recursive_copy')
    call disp%show("isDir(cpath)")
    call disp%show( isDir(cpath) )
    call disp%show("isFailedCopy(from = dir, to = cpath, recursive = .true._LK)")
    call disp%show( isFailedCopy(from = dir, to = cpath, recursive = .true._LK) )
    call disp%show("isDir(cpath//getDirSep()//dir)")
    call disp%show( isDir(cpath//getDirSep()//dir) )
    call disp%skip()

end program example