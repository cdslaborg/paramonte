program example

    use pm_kind, only: SK, IK, LK, RKD
    use pm_io, only: display_type

    implicit none

    type(display_type)  :: disp
    disp = display_type(file = "main.out.F90")

#if defined _OPENMP
    block
        use pm_timer, only: getResTimerOMP
        call disp%skip()
        call disp%show("getResTimerOMP()")
        call disp%show( getResTimerOMP() )
        call disp%skip()
    end block
#else
    call disp%show("The OMP timer is unavailable, because this example has not been linked with an OpenMP library.")
#endif

end program example