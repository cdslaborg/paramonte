program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysShell, only: isFailedExec
    use pm_sysPath, only: isFile, getPathNew

    implicit none

    character(:, SK), allocatable  :: file
    character(255, SK)  :: cmdmsg = SK_""
    integer(IK)         :: cmdstat
    integer(IK)         :: exitstat = -huge(0_IK)
    logical(LK)         :: cmdExecFailed

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("cmdExecFailed = isFailedExec(SK_'')")
                    cmdExecFailed = isFailedExec(SK_'')
    call disp%show("cmdExecFailed")
    call disp%show( cmdExecFailed )
    call disp%skip()

    call disp%skip()
    call disp%show("cmdExecFailed = isFailedExec(SK_'ls')")
                    cmdExecFailed = isFailedExec(SK_'ls')
    call disp%show("cmdExecFailed")
    call disp%show( cmdExecFailed )
    call disp%skip()

    call disp%skip()
    call disp%show("file = getPathNew() ! Create a file name unique in the current directory.")
                    file = getPathNew()
    call disp%show("file")
    call disp%show( file , deliml = SK_"""" )
    call disp%show("isFile(file)")
    call disp%show( isFile(file) )
    call disp%show("cmdExecFailed = isFailedExec(SK_'echo 2>'//file) ! Create an empty file in the current directory.")
                    cmdExecFailed = isFailedExec(SK_'echo 2>'//file)
    call disp%show("cmdExecFailed")
    call disp%show( cmdExecFailed )
    call disp%show("isFile(file)")
    call disp%show( isFile(file) )
    call disp%skip()

    call disp%skip()
    call disp%show("cmdExecFailed = isFailedExec(SK_'ls', wait = .false._LK) ! run the Bash command asynchronously.")
                    cmdExecFailed = isFailedExec(SK_'ls', wait = .false._LK)
    call disp%show("cmdExecFailed")
    call disp%show( cmdExecFailed )
    call disp%skip()

    call disp%skip()
    call disp%show("cmdExecFailed = isFailedExec(SK_'', exitstat = exitstat, cmdstat = cmdstat, cmdmsg = cmdmsg)")
                    cmdExecFailed = isFailedExec(SK_'', exitstat = exitstat, cmdstat = cmdstat, cmdmsg = cmdmsg)
    call disp%show("cmdExecFailed")
    call disp%show( cmdExecFailed )
    call disp%show("exitstat")
    call disp%show( exitstat )
    call disp%show("cmdstat")
    call disp%show( cmdstat )
    call disp%show("trim(cmdmsg)")
    call disp%show( trim(cmdmsg) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("cmdExecFailed = isFailedExec(SK_'dir', exitstat = exitstat, cmdstat = cmdstat, cmdmsg = cmdmsg)")
                    cmdExecFailed = isFailedExec(SK_'dir', exitstat = exitstat, cmdstat = cmdstat, cmdmsg = cmdmsg)
    call disp%show("cmdExecFailed")
    call disp%show( cmdExecFailed )
    call disp%show("exitstat")
    call disp%show( exitstat )
    call disp%show("cmdstat")
    call disp%show( cmdstat )
    call disp%show("trim(cmdmsg)")
    call disp%show( trim(cmdmsg) , deliml = SK_"""" )
    call disp%skip()

end program example