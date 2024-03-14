program example

    use pm_kind, only: SK, IK
    use pm_distBand, only: getBandEpeak
    use pm_distBand, only: getBandEbreak
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getBandEpeak(alpha = -1., beta = -3., ebreak = 300.)")
    call disp%show( getBandEpeak(alpha = -1., beta = -3., ebreak = 300.) )
    call disp%skip()

    call disp%skip()
    call disp%show("getBandEbreak(alpha = -1., beta = -3., epeak = getBandEpeak(alpha = 1., beta = -3., ebreak = 300.))")
    call disp%show( getBandEbreak(alpha = -1., beta = -3., epeak = getBandEpeak(alpha = 1., beta = -3., ebreak = 300.)) )
    call disp%skip()

    call disp%skip()
    call disp%show("getBandEpeak(alpha = -1., beta = -2., ebreak = [real :: 100, 200, 300])")
    call disp%show( getBandEpeak(alpha = -1., beta = -2., ebreak = [real :: 100, 200, 300]) )
    call disp%skip()

end program example