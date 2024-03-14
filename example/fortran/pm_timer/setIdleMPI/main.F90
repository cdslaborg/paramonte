program example

    use pm_kind, only: SK, IK, LK, RKD
    use pm_io, only: display_type

    implicit none


    type(display_type)  :: disp
    disp = display_type(file = "main.out.F90")

#if defined MPI_VERSION
    block
        use pm_timer, only: getTimeMPI
        use pm_timer, only: setIdleMPI
        real(RKD) :: since
        call disp%skip()
        call disp%show("since = getTimeMPI()")
                        since = getTimeMPI()
        call disp%show("call setIdleMPI(seconds = 0.2_RKD)")
                        call setIdleMPI(seconds = 0.2_RKD)
        call disp%show("getTimeMPI(since)")
        call disp%show( getTimeMPI(since) )
        call disp%skip()
    end block
#else
    call disp%show("The MPI timer is unavailable, because this example has not been linked with an MPI library.")
#endif

end program example