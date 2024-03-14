program example

    use pm_kind, only: SK, IK
    use iso_fortran_env, only: int8
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
    call disp%show("127_int8 .uadd. [huge(0_int8), 127_int8, 1_int8]")
    call disp%show( 127_int8 .uadd. [huge(0_int8), 127_int8, 1_int8] )
    call disp%skip

end program example