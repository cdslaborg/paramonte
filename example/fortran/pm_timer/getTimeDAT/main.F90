program example

    use pm_kind, only: SK, IK, LK, RKD
    use pm_io, only: display_type
    use pm_timer, only: getTimeDAT
    use pm_timer, only: setIdleDAT

    implicit none

    real(RKD)           :: since

    type(display_type)  :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("since = getTimeDAT()")
                    since = getTimeDAT()
    call disp%show("call setIdleDAT(seconds = 0.1_RKD)")
                    call setIdleDAT(seconds = 0.1_RKD)
    call disp%show("getTimeDAT(since)")
    call disp%show( getTimeDAT(since) )
    call disp%skip()

end program example