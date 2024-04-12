program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_dateTime, only: dateTimeInt_type

    implicit none

    type(dateTimeInt_type) :: dti

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct a Gregorian calendar date and time.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("dti = dateTimeInt_type(year = 2020_IK, zone = -5 * 60_IK, minute = 48_IK)")
                    dti = dateTimeInt_type(year = 2020_IK, zone = -5 * 60_IK, minute = 48_IK)
    ! Intel compiler bug on Windows OS: forrtl: severe (32): invalid logical unit number, unit -129, file unknown
#if !__INTEL_COMPILER || !_WIN32
    call disp%show("write(disp%unit, ""(*(g0,:,', '))"") dti")
                    write(disp%unit, "(*(g0,:,', '))") dti
#endif
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct the current date and time and moment in the Gregorian calendar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("dti = dateTimeInt_type()")
                    dti = dateTimeInt_type()
    ! Intel compiler bug on Windows OS: forrtl: severe (32): invalid logical unit number, unit -129, file unknown
#if !__INTEL_COMPILER || !_WIN32
    call disp%show("write(disp%unit, ""(*(g0,:,', '))"") dti")
                    write(disp%unit, "(*(g0,:,', '))") dti
#endif
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Get the object components as a vector of `Values` as returned by the Fortran instrinsic `date_and_time()`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("dti = dateTimeInt_type(year = 2020_IK, zone = -5 * 60_IK, minute = 48_IK)")
                    dti = dateTimeInt_type(year = 2020_IK, zone = -5 * 60_IK, minute = 48_IK)
    ! Intel compiler bug on Windows OS: forrtl: severe (32): invalid logical unit number, unit -129, file unknown
#if !__INTEL_COMPILER || !_WIN32
    call disp%show("write(disp%unit, ""(*(g0,:,', '))"") dti")
                    write(disp%unit, "(*(g0,:,', '))") dti
#endif
    call disp%show("dti%getValues()")
    call disp%show( dti%getValues() )
    call disp%show("size(dti%getValues())")
    call disp%show( size(dti%getValues()) )
    call disp%skip()

end program example