program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysPath, only: getPathSep
    use pm_sysPath, only: setPathMatch
    use pm_arrayReplace, only: getReplaced

    implicit none

    integer(IK) :: lenList
    type(display_type) :: disp
    character(:, SK), allocatable :: list
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("lenList = 9999")
                    lenList = 9999
    call disp%show("list = repeat(' ', lenList)")
                    list = repeat(' ', lenList)
    call disp%show("call setPathMatch(key = SK_'', inc = SK_'', sep = getPathSep(), list = list, lenList = lenList)")
                    call setPathMatch(key = SK_'', inc = SK_'', sep = getPathSep(), list = list, lenList = lenList)
    call disp%show("lenList")
    call disp%show( lenList )
    call disp%show("getReplaced(trim(list), getPathSep(), new_line(SK_'a'))")
    call disp%show( getReplaced(trim(list), getPathSep(), new_line(SK_'a')) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("lenList = 9999")
                    lenList = 9999
    call disp%show("list = repeat(' ', lenList)")
                    list = repeat(' ', lenList)
    call disp%show("call setPathMatch(key = SK_'intel:oneapi:mpi:', inc = SK_'', sep = getPathSep(), list = list, lenList = lenList)")
                    call setPathMatch(key = SK_'intel:oneapi:mpi:', inc = SK_'', sep = getPathSep(), list = list, lenList = lenList)
    call disp%show("lenList")
    call disp%show( lenList )
    call disp%show("getReplaced(trim(list), SK_':', new_line(SK_'a'))")
    call disp%show( getReplaced(trim(list), SK_':', new_line(SK_'a')) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("lenList = 9999")
                    lenList = 9999
    call disp%show("list = repeat(' ', lenList)")
                    list = repeat(' ', lenList)
    call disp%show("call setPathMatch(key = SK_'intel:oneapi:mpi', inc = SK_'mpiexec', sep = getPathSep(), list = list, lenList = lenList)")
                    call setPathMatch(key = SK_'intel:oneapi:mpi', inc = SK_'mpiexec', sep = getPathSep(), list = list, lenList = lenList)
    call disp%show("lenList")
    call disp%show( lenList )
    call disp%show("getReplaced(trim(list), getPathSep(), new_line(SK_'a'))")
    call disp%show( getReplaced(trim(list), getPathSep(), new_line(SK_'a')) , deliml = SK_"""" )
    call disp%skip()

end program example