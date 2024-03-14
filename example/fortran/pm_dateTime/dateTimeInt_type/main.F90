program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_dateTime, only: dateTimeInt_type

    implicit none

    type(dateTimeInt_type) :: DTI

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct a Gregorian calendar date and time.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("DTI = dateTimeInt_type(year = 2020_IK, zone = -5 * 60_IK, minute = 48_IK)")
                    DTI = dateTimeInt_type(year = 2020_IK, zone = -5 * 60_IK, minute = 48_IK)
    call disp%show("write(disp%unit, ""(*(g0,:,', '))"") DTI")
                    write(disp%unit, "(*(g0,:,', '))") DTI
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct the current date and time and moment in the Gregorian calendar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("DTI = dateTimeInt_type()")
                    DTI = dateTimeInt_type()
    call disp%show("write(disp%unit, ""(*(g0,:,', '))"") DTI")
                    write(disp%unit, "(*(g0,:,', '))") DTI
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Get the object components as a vector of `Values` as returned by the Fortran instrinsic `date_and_time()`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("DTI = dateTimeInt_type(year = 2020_IK, zone = -5 * 60_IK, minute = 48_IK)")
                    DTI = dateTimeInt_type(year = 2020_IK, zone = -5 * 60_IK, minute = 48_IK)
    call disp%show("write(disp%unit, ""(*(g0,:,', '))"") DTI")
                    write(disp%unit, "(*(g0,:,', '))") DTI
    call disp%show("DTI%getValues()")
    call disp%show( DTI%getValues() )
    call disp%show("size(DTI%getValues())")
    call disp%show( size(DTI%getValues()) )
    call disp%skip()

end program example