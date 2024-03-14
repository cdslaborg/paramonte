program example

    use pm_kind, only: SK, IK, LK, RKD
    use pm_io, only: display_type
    use pm_timer, only: getTimeSYS
    use pm_timer, only: setIdleSYS

    implicit none

    real(RKD)           :: since

    type(display_type)  :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("since = getTimeSYS()")
                    since = getTimeSYS()
    call disp%show("call setIdleSYS(seconds = 0.2_RKD)")
                    call setIdleSYS(seconds = 0.2_RKD)
    call disp%show("getTimeSYS(since)")
    call disp%show( getTimeSYS(since) )
    call disp%skip()

end program example