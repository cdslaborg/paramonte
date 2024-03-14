program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! All kinds are supported.
    use pm_kind, only: IK ! All kinds are supported.
    use pm_kind, only: CK ! All kinds are supported.
    use pm_kind, only: RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayFill, only: getFilled

    implicit none

    integer(IK) :: i
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%show("getFilled(0., 5)")
    call disp%show( getFilled(0., 5) , deliml = SK_"""" )
    call disp%skip()

    call disp%show("getFilled(0.d0, 5, 4, 3)")
    call disp%show( getFilled(0.d0, 5, 4, 3) , deliml = SK_"""" )
    call disp%skip()

    call disp%show("getFilled('AZ', 4, 5)")
    call disp%show( getFilled('AZ', 4, 5) , deliml = SK_"""" )
    call disp%skip()

    call disp%show("getFilled(.false., 4, 5, 0)")
    call disp%show( getFilled(.false., 4, 5, 0) , deliml = SK_"""" )
    call disp%skip()

    call disp%show("getFilled(.false., 4, 5, 1)")
    call disp%show( getFilled(.false., 4, 5, 1) , deliml = SK_"""" )
    call disp%skip()

    call disp%show("getFilled(333, 2, 5)")
    call disp%show( getFilled(333, 2, 5) , deliml = SK_"""" )
    call disp%skip()

    call disp%show("getFilled((0., 1.), 5, 3)")
    call disp%show( getFilled((0., 1.), 5, 3) , deliml = SK_"""" )
    call disp%skip()

end program example