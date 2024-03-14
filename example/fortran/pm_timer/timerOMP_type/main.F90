program example

    use pm_kind, only: IK, LK, SK, RKD
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

#if !defined _OPENMP
    call disp%show("The OMP timer is unavailable, because this example has not been linked with an OpenMP library.")
#else
    block
        use pm_timer, only: timerOMP_type
        type(timerOMP_type) :: timer
        logical(LK)         :: failed
        integer             :: i
    
        call disp%skip()
        call disp%show("timer = timerOMP_type()")
                        timer = timerOMP_type()
        call showTimerComponents()
        call disp%skip()
    
        call disp%skip()
        call disp%show("timer%clock = timer%time()")
                        timer%clock = timer%time()
        call showTimerComponents()
        call disp%skip()
    
        call disp%skip()
        call disp%show("timer%delta = timer%time(since = timer%start)")
                        timer%delta = timer%time(since = timer%start)
        call showTimerComponents()
        call disp%skip()
    
        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!Idle for a certain number of seconds.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()
    
        call disp%skip()
        call disp%show("timer = timerOMP_type()")
                        timer = timerOMP_type()
        call disp%show("call timer%wait(seconds = 0.1_RKD)")
                        call timer%wait(seconds = 0.1_RKD)
        call disp%show("timer%delta = timer%time() - timer%start")
                        timer%delta = timer%time() - timer%start
        call showTimerComponents()
        call disp%skip()
    end block
contains
    impure subroutine showTimerComponents()
        call disp%show("timer%start")
        call disp%show( timer%start )
        call disp%show("timer%clock")
        call disp%show( timer%clock )
        call disp%show("timer%delta")
        call disp%show( timer%delta )
        call disp%show("timer%resol")
        call disp%show( timer%resol )
    end subroutine
#endif

end program example