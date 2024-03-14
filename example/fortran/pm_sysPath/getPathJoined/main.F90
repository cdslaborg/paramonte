program example

    use pm_kind, only: SK, LK
    use pm_io, only: display_type
    use pm_sysPath, only: getPathJoined
    use pm_sysShell, only: isShellWindows
    use pm_sysShell, only: isShellPosix

    implicit none

    character(:, SK), allocatable :: pathJoined

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Note the runtime shell type below before checking the example outputs.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("isShellPosix()")
    call disp%show( isShellPosix() )
    call disp%show("isShellWindows()")
    call disp%show( isShellWindows() )
    call disp%skip()

    call disp%skip()
    call disp%show("pathJoined = getPathJoined(SK_'', SK_'')")
                    pathJoined = getPathJoined(SK_'', SK_'')
    call disp%show("pathJoined")
    call disp%show( pathJoined , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("pathJoined = getPathJoined(SK_'C:', SK_'\')")
                    pathJoined = getPathJoined(SK_'C:', SK_'\')
    call disp%show("pathJoined")
    call disp%show( pathJoined , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("pathJoined = getPathJoined(SK_'C:', SK_'/')")
                    pathJoined = getPathJoined(SK_'C:', SK_'/')
    call disp%show("pathJoined")
    call disp%show( pathJoined , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("pathJoined = getPathJoined(SK_'C:', SK_'./')")
                    pathJoined = getPathJoined(SK_'C:', SK_'./')
    call disp%show("pathJoined")
    call disp%show( pathJoined , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("pathJoined = getPathJoined(SK_'C:', SK_'paramonte')")
                    pathJoined = getPathJoined(SK_'C:', SK_'paramonte')
    call disp%show("pathJoined")
    call disp%show( pathJoined , deliml = SK_'''' )
    call disp%skip()

    call disp%skip()
    call disp%show("pathJoined = getPathJoined(SK_'./', SK_'paramonte')")
                    pathJoined = getPathJoined(SK_'./', SK_'paramonte')
    call disp%show("pathJoined")
    call disp%show( pathJoined , deliml = SK_'''' )
    call disp%skip()

end program example