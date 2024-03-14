program example

    use pm_kind, only: SK, IK, LK, RKD
    use pm_io, only: display_type
    use pm_timer, only: getResTimerCPU

    implicit none

    type(display_type)  :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getResTimerCPU()")
    call disp%show( getResTimerCPU() )
    call disp%skip()

end program example