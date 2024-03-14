program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RK, RK32, RK64, RK128
    use pm_io, only: display_type
    use pm_distBern, only: isHead

    implicit none

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Flip a potentially biased coin (as determined by the success rate `p`).")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("isHead() ! flip a fair coin.")
    call disp%show( isHead() )
    call disp%skip()

    call disp%show("isHead(size = 5_IK) ! flip a fair coin 5 times.")
    call disp%show( isHead(size = 5_IK) )
    call disp%skip()

    call disp%show("isHead(p = 0.5_RK32) ! flip a fair coin with real-32bit precision.")
    call disp%show( isHead(p = 0.5_RK32) )
    call disp%skip()

    call disp%show("isHead(p = 0.5_RK64) ! flip a fair coin with real-64bit precision.")
    call disp%show( isHead(p = 0.5_RK64) )
    call disp%skip()

    call disp%show("isHead(p = 0.5_RK128) ! flip a fair coin with real-128bit precision.")
    call disp%show( isHead(p = 0.5_RK128) )
    call disp%skip()

    call disp%show("isHead(0.25) ! flip an unfair coin with 0.25 odds of success.")
    call disp%show( isHead(0.25) )
    call disp%skip()

    call disp%show("isHead(0.25, size = 5_IK) ! flip an unfair coin with 0.25 odds of success 5 times.")
    call disp%show( isHead(0.25, size = 5_IK) )
    call disp%skip()

    call disp%show("isHead(p = [0.05, 0.25, 0.75, 0.95]) ! flip four different unfair coins.")
    call disp%show( isHead(p = [0.05, 0.25, 0.75, 0.95]) )
    call disp%skip()

end program example