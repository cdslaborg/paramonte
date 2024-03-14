program example

    use pm_kind, only: SK, IK, LK, RKD
    use pm_io, only: display_type
    use pm_timer, only: getTime
    use pm_timer, only: setIdle

    implicit none

    real(RKD)           :: since

    type(display_type)  :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("since = getTime()")
                    since = getTime()
    call disp%show("call setIdle(seconds = 0.2_RKD)")
                    call setIdle(seconds = 0.2_RKD)
    call disp%show("getTime(since)")
    call disp%show( getTime(since) )
    call disp%skip()

end program example