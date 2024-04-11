program example

    use pm_kind, only: SK, IK, IKL
    use pm_mathUnsigned, only: operator(.uadd.)
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    disp = display_type(file = SK_"main.out.F90")

    call disp%skip
    call disp%show("1 .uadd. 1")
    call disp%show( 1 .uadd. 1 )
    call disp%skip

    call disp%skip
    call disp%show("huge(1) .uadd. [huge(1), huge(1) / 2, 1, 0]")
    call disp%show( huge(1) .uadd. [huge(1), huge(1) / 2, 1, 0] )
    call disp%skip

    call disp%skip
    call disp%show("127_IKL .uadd. [huge(0_IKL), 127_IKL, 1_IKL]")
    call disp%show( 127_IKL .uadd. [huge(0_IKL), 127_IKL, 1_IKL] )
    call disp%skip

end program example