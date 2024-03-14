program example

    use pm_kind, only: SK, IK, LK, RKD
    use pm_io, only: display_type
    use pm_timer, only: getTimeCPU
    use pm_timer, only: setIdleCPU

    implicit none

    real(RKD)           :: since

    type(display_type)  :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("since = getTimeCPU()")
                    since = getTimeCPU()
    call disp%show("call setIdleCPU(seconds = 0.1_RKD)")
                    call setIdleCPU(seconds = 0.1_RKD)
    call disp%show("getTimeCPU(since)")
    call disp%show( getTimeCPU(since) )
    call disp%skip()

end program example