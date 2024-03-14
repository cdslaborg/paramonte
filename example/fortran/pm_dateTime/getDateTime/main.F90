program example

    use pm_kind, only: IK, RK, SK
    use pm_io, only: display_type
    use pm_dateTime, only: getDateTime

    implicit none

    integer(IK) :: Values(8)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate and return the current local Gregorian date and time.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime()")
    call disp%show( getDateTime() )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate and return the specified Gregorian date and time as vector.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(2000_IK)")
    call disp%show( getDateTime(2000_IK) )
    call disp%show("getDateTime(2000_IK, 12_IK)")
    call disp%show( getDateTime(2000_IK, 12_IK) )
    call disp%show("getDateTime(2000_IK, 12_IK, 13_IK)")
    call disp%show( getDateTime(2000_IK, 12_IK, 13_IK) )
    call disp%show("getDateTime(2000_IK, 12_IK, 13_IK, -420_IK)")
    call disp%show( getDateTime(2000_IK, 12_IK, 13_IK, -420_IK) )
    call disp%show("getDateTime(2000_IK, 12_IK, 13_IK, -420_IK, 18)")
    call disp%show( getDateTime(2000_IK, 12_IK, 13_IK, -420_IK, 18) )
    call disp%show("getDateTime(2000_IK, 12_IK, 13_IK, -420_IK, 18, 21)")
    call disp%show( getDateTime(2000_IK, 12_IK, 13_IK, -420_IK, 18, 21) )
    call disp%show("getDateTime(2000_IK, 12_IK, 13_IK, -420_IK, 18, 21, 35)")
    call disp%show( getDateTime(2000_IK, 12_IK, 13_IK, -420_IK, 18, 21, 35) )
    call disp%show("getDateTime(2000_IK, 12_IK, 13_IK, -420_IK, 18, 21, 35, 847)")
    call disp%show( getDateTime(2000_IK, 12_IK, 13_IK, -420_IK, 18, 21, 35, 847) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate and return the UTC Gregorian date and time as a vector corresponding to the input Julian day.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("getDateTime(1720694.5_RK) ! January 1, -1, (2 BC) start of day.")
    call disp%show( getDateTime(1720694.5_RK) )

    call disp%show("getDateTime(1720694.5_RK + 365._RK) ! January 1, -1, (1 BC) start of day.")
    call disp%show( getDateTime(1720694.5_RK + 365._RK) )

    call disp%show("getDateTime(1720694.5_RK + 365._RK + 366._RK) ! January 1, 1, (AD 1) start of day.")
    call disp%show( getDateTime(1720694.5_RK + 365._RK + 366._RK) )

    call disp%show("getDateTime(2299160.5_RK) ! October 15, 1582, start of day (first day of Gregorian reform).")
    call disp%show( getDateTime(2299160.5_RK) )

    call disp%show("getDateTime(2415385.5_RK) ! January 1, 1901, start of day (start of the 20th century).")
    call disp%show( getDateTime(2415385.5_RK) )

    call disp%show("getDateTime(2440587.5_RK) ! January 1, 1970, start of day (Unix reference date).")
    call disp%show( getDateTime(2440587.5_RK) )

    call disp%show("getDateTime(2444239.5_RK) ! January 1, 1980, start of day (Unix reference date).")
    call disp%show( getDateTime(2444239.5_RK) )

    call disp%show("getDateTime(2444240.0_RK) ! January 1, 1980, noon (Unix reference date).")
    call disp%show( getDateTime(2444240.0_RK) )

    call disp%show("getDateTime(2451545.25_RK) ! January 1, 2000, 18:00 UTC.")
    call disp%show( getDateTime(2451545.25_RK) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate and return the Gregorian date and time as a string.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'')")
    call disp%show( getDateTime(SK_'') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'Today is %A.')")
    call disp%show( getDateTime(SK_'Today is %A.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'%A is abbreviated as %a.')")
    call disp%show( getDateTime(SK_'%A is abbreviated as %a.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'Today is a %A of month %B.')")
    call disp%show( getDateTime(SK_'Today is a %A of month %B.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'%B is abbreviated as %b.')")
    call disp%show( getDateTime(SK_'%B is abbreviated as %b.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'Current Date and Time is %c')")
    call disp%show( getDateTime(SK_'Current Date and Time is %c') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The Julian Period starts at %c of the Gregorian Calendar.', getDateTime(-4713_IK, 1_IK, 1_IK, 0_IK, 12_IK))")
    call disp%show( getDateTime(SK_'The Julian Period starts at %c of the Gregorian Calendar.', getDateTime(-4713_IK, 1_IK, 1_IK, 0_IK, 12_IK)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The last century was %C.')")
    call disp%show( getDateTime(SK_'The last century was %C.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The last century BC was %C.', getDateTime(year = -1_IK))")
    call disp%show( getDateTime(SK_'The last century BC was %C.', getDateTime(year = -1_IK)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The Holocene geological epoch roughly began at century %C.', getDateTime(year = -12000_IK))")
    call disp%show( getDateTime(SK_'The Holocene geological epoch roughly began at century %C.', getDateTime(year = -12000_IK)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The current day of month is %d.')")
    call disp%show( getDateTime(SK_'The current day of month is %d.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The first day of every month (zero-padded) is %d.', getDateTime(0_IK))")
    call disp%show( getDateTime(SK_'The first day of every month (zero-padded) is %d.', getDateTime(0_IK)) , deliml = SK_"""" )
    call disp%skip()
    call disp%skip()
    call disp%show("getDateTime(SK_'The first day of every month (blank-padded) is %e.', getDateTime(0_IK))")
    call disp%show( getDateTime(SK_'The first day of every month (blank-padded) is %e.', getDateTime(0_IK)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'March 21, 1987 in short (slashed) format is %D.', getDateTime(1987_IK, 3_IK, 21_IK))")
    call disp%show( getDateTime(SK_'March 21, 1987 in short (slashed) format is %D.', getDateTime(1987_IK, 3_IK, 21_IK)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The millisecond of the moment now (possibly padded with leading zeros) %f.')")
    call disp%show( getDateTime(SK_'The millisecond of the moment now (possibly padded with leading zeros) %f.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'March 21, 1987 in short (dashed) format is %F.', getDateTime(-1987_IK, 3_IK, 21_IK))")
    call disp%show( getDateTime(SK_'March 21, 1987 in short (dashed) format is %F.', getDateTime(-1987_IK, 3_IK, 21_IK)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The last two digits of the week year of the Gregorian Calendar date 1 Jan 1977 are %g.', getDateTime(1977_IK, 1_IK, 1_IK))")
    call disp%show( getDateTime(SK_'The last two digits of the week year of the Gregorian Calendar date 1 Jan 1977 are %g.', getDateTime(1977_IK, 1_IK, 1_IK)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The full week year of the Gregorian Calendar date 1 Jan 1977 is %G.', getDateTime(1977_IK, 1_IK, 1_IK))")
    call disp%show( getDateTime(SK_'The full week year of the Gregorian Calendar date 1 Jan 1977 is %G.', getDateTime(1977_IK, 1_IK, 1_IK)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The locale abbreviation for month %B is %h.')")
    call disp%show( getDateTime(SK_'The locale abbreviation for month %B is %h.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The hour 6 pm in 24-hour format is %H.', getDateTime(1_IK, 1_IK, 1_IK, 0_IK, 18_IK))")
    call disp%show( getDateTime(SK_'The hour 6 pm in 24-hour format is %H.', getDateTime(1_IK, 1_IK, 1_IK, 0_IK, 18_IK)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The hour 6 pm in 12-hour format is %I.', getDateTime(1_IK, 1_IK, 1_IK, 0_IK, 18_IK))")
    call disp%show( getDateTime(SK_'The hour 6 pm in 12-hour format is %I.', getDateTime(1_IK, 1_IK, 1_IK, 0_IK, 18_IK)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The current day of year (ordinal day) is %j.')")
    call disp%show( getDateTime(SK_'The current day of year (ordinal day) is %j.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The day of year (ordinal day) corresponding to March 21, 1987 is %j.', getDateTime(1987_IK, 3_IK, 21_IK))")
    call disp%show( getDateTime(SK_'The day of year (ordinal day) corresponding to March 21, 1987 is %j.', getDateTime(1987_IK, 3_IK, 21_IK)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The month of March 21, 1987 as a decimal number is %m.', getDateTime(1987_IK, 3_IK, 21_IK))")
    call disp%show( getDateTime(SK_'The month of March 21, 1987 as a decimal number is %m.', getDateTime(1987_IK, 3_IK, 21_IK)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The current minute of time is %M.')")
    call disp%show( getDateTime(SK_'The current minute of time is %M.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The specifier %%n yields a new line character: %n.')")
    call disp%show( getDateTime(SK_'The specifier %%n yields a new line character: %n.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The current moment is in %p period of day.')")
    call disp%show( getDateTime(SK_'The current moment is in %p period of day.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The current 12-hour clock time is %r.')")
    call disp%show( getDateTime(SK_'The current 12-hour clock time is %r.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The current 24-hour clock time is %r.')")
    call disp%show( getDateTime(SK_'The current 24-hour clock time is %r.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The initial hour and minute setting of the Doomsday Clock was %R.', getDateTime(1947, 1_IK, 1_IK, 0_IK, 23_IK, 53_IK))")
    call disp%show( getDateTime(SK_'The initial hour and minute setting of the Doomsday Clock was %R.', getDateTime(1947, 1_IK, 1_IK, 0_IK, 23_IK, 53_IK)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The second of the moment now is %S.')")
    call disp%show( getDateTime(SK_'The second of the moment now is %S.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The specifier %%t yields a horizontal tab character: %t.')")
    call disp%show( getDateTime(SK_'The specifier %%t yields a horizontal tab character: %t.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The current local time in ISO 8601 format is %T.')")
    call disp%show( getDateTime(SK_'The current local time in ISO 8601 format is %T.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The current ISO 8601 weekday is %u.')")
    call disp%show( getDateTime(SK_'The current ISO 8601 weekday is %u.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The current ISO 8601 week number is %V.')")
    call disp%show( getDateTime(SK_'The current ISO 8601 week number is %V.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The current weekday as a decimal number with Sunday as 0 is %w.')")
    call disp%show( getDateTime(SK_'The current weekday as a decimal number with Sunday as 0 is %w.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'March 21, 1987 in locale format is %x.', getDateTime(1987_IK, 3_IK, 21_IK))")
    call disp%show( getDateTime(SK_'March 21, 1987 in locale format is %x.', getDateTime(1987_IK, 3_IK, 21_IK)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The current time in locale format is %X.')")
    call disp%show( getDateTime(SK_'The current time in locale format is %X.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The last two digits of the current year are %y.')")
    call disp%show( getDateTime(SK_'The last two digits of the current year are %y.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The last two digits of year 1987 are %y.', getDateTime(1987_IK))")
    call disp%show( getDateTime(SK_'The last two digits of year 1987 are %y.', getDateTime(1987_IK)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The current local ISO 8601 offset from UTC in timezone in units of minutes is %z.')")
    call disp%show( getDateTime(SK_'The current local ISO 8601 offset from UTC in timezone in units of minutes is %z.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The current local timezone abbreviation is %Z.')")
    call disp%show( getDateTime(SK_'The current local timezone abbreviation is %Z.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The timezone abbreviations corresponding to ISO timezone UTC-12 is %Z.', getDateTime(1_IK, 1_IK, 1_IK, -12_IK * 60_IK))")
    call disp%show( getDateTime(SK_'The timezone abbreviations corresponding to ISO timezone UTC-12 is %Z.', getDateTime(1_IK, 1_IK, 1_IK, -12_IK * 60_IK)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The timezone abbreviations corresponding to ISO timezone UTC-11 is %Z.', getDateTime(1_IK, 1_IK, 1_IK, -11_IK * 60_IK))")
    call disp%show( getDateTime(SK_'The timezone abbreviations corresponding to ISO timezone UTC-11 is %Z.', getDateTime(1_IK, 1_IK, 1_IK, -11_IK * 60_IK)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The timezone abbreviations corresponding to ISO timezone UTC-10 is %Z.', getDateTime(1_IK, 1_IK, 1_IK, -10_IK * 60_IK))")
    call disp%show( getDateTime(SK_'The timezone abbreviations corresponding to ISO timezone UTC-10 is %Z.', getDateTime(1_IK, 1_IK, 1_IK, -10_IK * 60_IK)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The timezone abbreviations corresponding to ISO timezone UTC-9:30 is %Z.', getDateTime(1_IK, 1_IK, 1_IK, -9_IK * 60_IK - 30_IK))")
    call disp%show( getDateTime(SK_'The timezone abbreviations corresponding to ISO timezone UTC-9:30 is %Z.', getDateTime(1_IK, 1_IK, 1_IK, -9_IK * 60_IK - 30_IK)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The timezone abbreviations corresponding to ISO timezone UTC-9 is %Z.', getDateTime(1_IK, 1_IK, 1_IK, -9_IK * 60_IK))")
    call disp%show( getDateTime(SK_'The timezone abbreviations corresponding to ISO timezone UTC-9 is %Z.', getDateTime(1_IK, 1_IK, 1_IK, -9_IK * 60_IK)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The timezone abbreviations corresponding to ISO timezone UTC-8:30 is %Z.', getDateTime(1_IK, 1_IK, 1_IK, -8_IK * 60_IK - 30_IK)) ! non-existent timezone yields empty abbreviation.")
    call disp%show( getDateTime(SK_'The timezone abbreviations corresponding to ISO timezone UTC-8:30 is %Z.', getDateTime(1_IK, 1_IK, 1_IK, -8_IK * 60_IK - 30_IK)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The timezone abbreviations corresponding to ISO timezone UTC+0 is %Z.', getDateTime(1_IK, 1_IK, 1_IK, +0_IK * 60_IK))")
    call disp%show( getDateTime(SK_'The timezone abbreviations corresponding to ISO timezone UTC+0 is %Z.', getDateTime(1_IK, 1_IK, 1_IK, +0_IK * 60_IK)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The timezone abbreviations corresponding to ISO timezone UTC+8 is %Z.', getDateTime(1_IK, 1_IK, 1_IK, +8_IK * 60_IK))")
    call disp%show( getDateTime(SK_'The timezone abbreviations corresponding to ISO timezone UTC+8 is %Z.', getDateTime(1_IK, 1_IK, 1_IK, +8_IK * 60_IK)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'The current year is %Y.')")
    call disp%show( getDateTime(SK_'The current year is %Y.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'%Q : An undefined specifier has no effects.')")
    call disp%show( getDateTime(SK_'%Q : An undefined specifier has no effects.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'%%A : two sequential percentage characters translate to a single percentage character.')")
    call disp%show( getDateTime(SK_'%%A : two sequential percentage characters translate to a single percentage character.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDateTime(SK_'% : dangling percentage character does nothing.')")
    call disp%show( getDateTime(SK_'% : dangling percentage character does nothing.') , deliml = SK_"""" )
    call disp%skip()

end program example