!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief This file contains the implementations of the tests of module [pm_dateTime](@ref pm_dateTime).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, March 22, 2012, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_dateTime) routines

    use pm_val2str, only: getStr
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_constructDateTimeInt

        type(dateTimeInt_type) :: DTI
        integer(IK) :: i, values_ref(8)
        integer(IK), allocatable :: Values(:)
        call date_and_time(values = values_ref)
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        values_ref = [integer(IK) :: 1, 1, 1, 0, 0, 0, 0, 0]
        assertion = assertion .and. logical(DTI%year        == values_ref(1), LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor with no input arguments must correctly initialize the `year` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTI%month       == values_ref(2), LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor with no input arguments must correctly initialize the `month` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTI%day         == values_ref(3), LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor with no input arguments must correctly initialize the `day` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTI%zone        == values_ref(4), LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor with no input arguments must correctly initialize the `zone` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTI%hour        == values_ref(5), LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor with no input arguments must correctly initialize the `hour` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTI%minute      == values_ref(6), LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor with no input arguments must correctly initialize the `minute` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTI%second      == values_ref(7), LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor with no input arguments must correctly initialize the `second` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTI%millisecond == values_ref(8), LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor with no input arguments must correctly initialize the `millisecond` component.", int(__LINE__, IK))
        Values = DTI%getValues()
        call test%assert(assertion, SK_"@dateTimeInt_type(): The size of the returned vector by `getValues()` method must be 8. size(DTI%getValues()) = "//getStr(size(DTI%getValues())), int(__LINE__, IK))
        do i = 1, size(values_ref)
            assertion = assertion .and. logical(Values(i) == values_ref(i), LK)
            call test%assert(assertion, SK_"@dateTimeInt_type(): The type constructor must correctly extract the `i`th component value. i, Values(i), values_ref(i)"//getStr([i, Values(i), values_ref(i)]), int(__LINE__, IK))
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        DTI = dateTimeInt_type  ( year          = values_ref(1) &
                                , month         = values_ref(2) &
                                , day           = values_ref(3) &
                                , zone          = values_ref(4) &
                                , hour          = values_ref(5) &
                                , minute        = values_ref(6) &
                                , second        = values_ref(7) &
                                , millisecond   = values_ref(8) &
                                )
        assertion = assertion .and. logical(DTI%year        == values_ref(1), LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor must correctly assign the `year` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTI%month       == values_ref(2), LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor must correctly assign the `month` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTI%day         == values_ref(3), LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor must correctly assign the `day` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTI%zone        == values_ref(4), LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor must correctly assign the `zone` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTI%hour        == values_ref(5), LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor must correctly assign the `hour` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTI%minute      == values_ref(6), LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor must correctly assign the `minute` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTI%second      == values_ref(7), LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor must correctly assign the `second` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTI%millisecond == values_ref(8), LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor must correctly assign the `millisecond` component.", int(__LINE__, IK))
        Values = DTI%getValues()
        assertion = assertion .and. logical(size(Values) == 8, LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The size of the returned vector by `getValues()` method must be 8. size(DTI%getValues()) = "//getStr(size(DTI%getValues())), int(__LINE__, IK))
        do i = 1, size(values_ref)
            assertion = assertion .and. logical(Values(i) == values_ref(i), LK)
            call test%assert(assertion, SK_"@dateTimeInt_type(): The type constructor must correctly extract the `i`th component value. i, Values(i), values_ref(i)"//getStr([i, Values(i), values_ref(i)]), int(__LINE__, IK))
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 10
            values_ref = getDateTime()
            DTI = dateTimeInt_type()
            if (DTI%year        == values_ref(1) .and. & ! LCOV_EXCL_LINE
                DTI%month       == values_ref(2) .and. & ! LCOV_EXCL_LINE
                DTI%day         == values_ref(3) .and. & ! LCOV_EXCL_LINE
                DTI%zone        == values_ref(4) .and. & ! LCOV_EXCL_LINE
                DTI%hour        == values_ref(5) .and. & ! LCOV_EXCL_LINE
                DTI%minute      == values_ref(6) .and. & ! LCOV_EXCL_LINE
                DTI%second      == values_ref(7)) exit
        end do
        assertion = assertion .and. logical(DTI%year        == values_ref(1), LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor must correctly assign the `year` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTI%month       == values_ref(2), LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor must correctly assign the `month` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTI%day         == values_ref(3), LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor must correctly assign the `day` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTI%zone        == values_ref(4), LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor must correctly assign the `zone` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTI%hour        == values_ref(5), LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor must correctly assign the `hour` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTI%minute      == values_ref(6), LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor must correctly assign the `minute` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTI%second      == values_ref(7), LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor must correctly assign the `second` component.", int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_constructDateTimeStr

        use pm_container, only: strc => css_type
        type(dateTimeStr_type) :: DTS
        type(strc), allocatable :: values_ref(:)
        character( 8, SK) :: date
        character(10, SK) :: time
        character( 5, SK) :: zone
        integer :: i
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        values_ref = [strc(SK_"0001"), strc(SK_"01"), strc(SK_"01"), strc(SK_"+0000"), strc(SK_"00"), strc(SK_"00"), strc(SK_"00"), strc(SK_"000")]
        assertion = assertion .and. logical(DTS%year        == values_ref(1)%val, LK)
        call test%assert(assertion, SK_"@dateTimeStr_type(): The default constructor with no input arguments must correctly initialize the `year` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTS%month       == values_ref(2)%val, LK)
        call test%assert(assertion, SK_"@dateTimeStr_type(): The default constructor with no input arguments must correctly initialize the `month` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTS%day         == values_ref(3)%val, LK)
        call test%assert(assertion, SK_"@dateTimeStr_type(): The default constructor with no input arguments must correctly initialize the `day` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTS%zone        == values_ref(4)%val, LK)
        call test%assert(assertion, SK_"@dateTimeStr_type(): The default constructor with no input arguments must correctly initialize the `zone` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTS%hour        == values_ref(5)%val, LK)
        call test%assert(assertion, SK_"@dateTimeStr_type(): The default constructor with no input arguments must correctly initialize the `hour` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTS%minute      == values_ref(6)%val, LK)
        call test%assert(assertion, SK_"@dateTimeStr_type(): The default constructor with no input arguments must correctly initialize the `minute` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTS%second      == values_ref(7)%val, LK)
        call test%assert(assertion, SK_"@dateTimeStr_type(): The default constructor with no input arguments must correctly initialize the `second` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTS%millisecond == values_ref(8)%val, LK)
        call test%assert(assertion, SK_"@dateTimeStr_type(): The default constructor with no input arguments must correctly initialize the `millisecond` component.", int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call date_and_time(date = date, zone = zone, time = time)
        values_ref(1)%val = date(1:4)
        values_ref(2)%val = date(5:6)
        values_ref(3)%val = date(7:8)
        values_ref(4)%val = zone
        values_ref(5)%val = time(1:2)
        values_ref(6)%val = time(3:4)
        values_ref(7)%val = time(5:6)
        values_ref(8)%val = time(8:10)
        DTS = dateTimeStr_type  ( year          = values_ref(1)%val &
                                , month         = values_ref(2)%val &
                                , day           = values_ref(3)%val &
                                , zone          = values_ref(4)%val &
                                , hour          = values_ref(5)%val &
                                , minute        = values_ref(6)%val &
                                , second        = values_ref(7)%val &
                                , millisecond   = values_ref(8)%val &
                                )
        assertion = assertion .and. logical(DTS%year        == values_ref(1)%val, LK)
        call test%assert(assertion, SK_"@dateTimeStr_type(): The default constructor must correctly assign the `year` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTS%month       == values_ref(2)%val, LK)
        call test%assert(assertion, SK_"@dateTimeStr_type(): The default constructor must correctly assign the `month` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTS%day         == values_ref(3)%val, LK)
        call test%assert(assertion, SK_"@dateTimeStr_type(): The default constructor must correctly assign the `day` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTS%zone        == values_ref(4)%val, LK)
        call test%assert(assertion, SK_"@dateTimeStr_type(): The default constructor must correctly assign the `zone` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTS%hour        == values_ref(5)%val, LK)
        call test%assert(assertion, SK_"@dateTimeStr_type(): The default constructor must correctly assign the `hour` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTS%minute      == values_ref(6)%val, LK)
        call test%assert(assertion, SK_"@dateTimeStr_type(): The default constructor must correctly assign the `minute` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTS%second      == values_ref(7)%val, LK)
        call test%assert(assertion, SK_"@dateTimeStr_type(): The default constructor must correctly assign the `second` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTS%millisecond == values_ref(8)%val, LK)
        call test%assert(assertion, SK_"@dateTimeStr_type(): The default constructor must correctly assign the `millisecond` component.", int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 10
            call date_and_time(date = date, zone = zone, time = time)
            DTS = dateTimeStr_type()
            values_ref(1)%val = date(1:4)
            values_ref(2)%val = date(5:6)
            values_ref(3)%val = date(7:8)
            values_ref(4)%val = zone
            values_ref(5)%val = time(1:2)
            values_ref(6)%val = time(3:4)
            values_ref(7)%val = time(5:6)
            values_ref(8)%val = time(8:10)
            if (DTS%year        == values_ref(1)%val .and. & ! LCOV_EXCL_LINE
                DTS%month       == values_ref(2)%val .and. & ! LCOV_EXCL_LINE
                DTS%day         == values_ref(3)%val .and. & ! LCOV_EXCL_LINE
                DTS%zone        == values_ref(4)%val .and. & ! LCOV_EXCL_LINE
                DTS%hour        == values_ref(5)%val .and. & ! LCOV_EXCL_LINE
                DTS%minute      == values_ref(6)%val .and. & ! LCOV_EXCL_LINE
                DTS%second      == values_ref(7)%val) exit
        end do
        assertion = assertion .and. logical(DTS%year        == values_ref(1)%val, LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor must correctly assign the `year` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTS%month       == values_ref(2)%val, LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor must correctly assign the `month` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTS%day         == values_ref(3)%val, LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor must correctly assign the `day` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTS%zone        == values_ref(4)%val, LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor must correctly assign the `zone` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTS%hour        == values_ref(5)%val, LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor must correctly assign the `hour` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTS%minute      == values_ref(6)%val, LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor must correctly assign the `minute` component.", int(__LINE__, IK))
        assertion = assertion .and. logical(DTS%second      == values_ref(7)%val, LK)
        call test%assert(assertion, SK_"@dateTimeInt_type(): The default constructor must correctly assign the `second` component.", int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_getJulianDay
        use pm_kind, only: RKG => RK, IKG => IK
        real(RKG) :: tol, julianDay, julianDay_ref
        real(RKG), parameter :: ABSTOL = epsilon(0._RKG) * 100
        integer(IKG) :: year, month, day, zone, hour, minute, second, millisecond
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        assertion = assertion .and. logical(abs(getJulianDay() - getJulianDay()) < 1._RKG / SECONDS_PER_DAY, LK)
        call test%assert(assertion, SK_"@test_getJulianDay(): The condition `getJulianDay() - getJulianDay() < 1._RKG / SECONDS_PER_DAY` must hold, unless the computation takes more than a few seconds!", int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tol = ABSTOL
        julianDay_ref = 1720694.5_RKG

        year        = -1_IKG
        month       = +1_IKG
        day         = +1_IKG
        zone        = +0_IKG
        hour        = +0_IKG
        minute      = +0_IKG
        second      = +0_IKG
        millisecond = +0_IKG

        call setAssertion(__LINE__, 1)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tol = ABSTOL
        julianDay_ref = -.5_RKG

        year        = -4713_IKG
        month       = +11_IKG
        day         = +24_IKG
        zone        = +0_IKG
        hour        = +0_IKG
        minute      = +0_IKG
        second      = +0_IKG
        millisecond = +0_IKG

        call setAssertion(__LINE__, 3)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tol = ABSTOL
        julianDay_ref = 0._RKG

        year        = -4713_IKG
        month       = +11_IKG
        day         = +24_IKG
        zone        = +0_IKG
        hour        = +12_IKG
        minute      = +0_IKG
        second      = +0_IKG
        millisecond = +0_IKG

        call setAssertion(__LINE__, 5)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tol = ABSTOL
        julianDay_ref = +.5_RKG

        year        = -4713_IKG
        month       = +11_IKG
        day         = +25_IKG
        zone        = +0_IKG
        hour        = +0_IKG
        minute      = +0_IKG
        second      = +0_IKG
        millisecond = +0_IKG

        call setAssertion(__LINE__, 3)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tol = ABSTOL
        julianDay_ref = +2299160.5_RKG

        year        = +1582_IKG
        month       = +10_IKG
        day         = +15_IKG
        zone        = +0_IKG
        hour        = +0_IKG
        minute      = +0_IKG
        second      = +0_IKG
        millisecond = +0_IKG

        call setAssertion(__LINE__, 3)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tol = ABSTOL
        julianDay_ref = 2415385.5_RKG

        year        = +1901_IKG
        month       = +1_IKG
        day         = +1_IKG
        zone        = +0_IKG
        hour        = +0_IKG
        minute      = +0_IKG
        second      = +0_IKG
        millisecond = +0_IKG

        call setAssertion(__LINE__, 1)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tol = ABSTOL
        julianDay_ref = 2440587.5_RKG

        year        = +1970_IKG
        month       = +1_IKG
        day         = +1_IKG
        zone        = +0_IKG
        hour        = +0_IKG
        minute      = +0_IKG
        second      = +0_IKG
        millisecond = +0_IKG

        call setAssertion(__LINE__, 1)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tol = ABSTOL
        julianDay_ref = 2444239.5_RKG

        year        = +1980_IKG
        month       = +1_IKG
        day         = +1_IKG
        zone        = +0_IKG
        hour        = +0_IKG
        minute      = +0_IKG
        second      = +0_IKG
        millisecond = +0_IKG

        call setAssertion(__LINE__, 1)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tol = ABSTOL
        julianDay_ref = 2444240.5_RKG

        year        = +1980_IKG
        month       = +1_IKG
        day         = +2_IKG
        zone        = +0_IKG
        hour        = +0_IKG
        minute      = +0_IKG
        second      = +0_IKG
        millisecond = +0_IKG

        call setAssertion(__LINE__, 3)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tol = ABSTOL
        julianDay_ref = 2444239.25_RKG

        year        = +1979_IKG
        month       = +12_IKG
        day         = +31_IKG
        zone        = +0_IKG
        hour        = +18_IKG
        minute      = +0_IKG
        second      = +0_IKG
        millisecond = +0_IKG

        call setAssertion(__LINE__, 5)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tol = ABSTOL
        julianDay_ref = 2444239._RKG

        year        = +1979_IKG
        month       = +12_IKG
        day         = +31_IKG
        zone        = +360_IKG
        hour        = +18_IKG
        minute      = +0_IKG
        second      = +0_IKG
        millisecond = +0_IKG

        call setAssertion(__LINE__, 5)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tol = ABSTOL
        julianDay_ref = 2444239._RKG

        year        = +1979_IKG
        month       = +12_IKG
        day         = +31_IKG
        zone        = -360_IKG
        hour        = +6_IKG
        minute      = +0_IKG
        second      = +0_IKG
        millisecond = +0_IKG

        call setAssertion(__LINE__, 5)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tol = ABSTOL
        julianDay_ref = 38245308.75_RKG

        year        = +99999_IKG
        month       = +12_IKG
        day         = +31_IKG
        zone        = +0_IKG
        hour        = +6_IKG
        minute      = +0_IKG
        second      = +0_IKG
        millisecond = +0_IKG

        call setAssertion(__LINE__, 5)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine setAssertion(line, minarg)

            integer, intent(in) :: line, minarg

            if (minarg > 8) return
            julianDay = getJulianDay(year, month, day, zone, hour, minute, second, millisecond)
            call report(line, 8)
            julianDay = getJulianDay([year, month, day, zone, hour, minute, second, millisecond])
            call report(line, 8)

            if (minarg > 7) return
            julianDay = getJulianDay(year, month, day, zone, hour, minute, second)
            call report(line, 7)
            julianDay = getJulianDay([year, month, day, zone, hour, minute, second])
            call report(line, 7)

            if (minarg > 6) return
            julianDay = getJulianDay(year, month, day, zone, hour, minute)
            call report(line, 6)
            julianDay = getJulianDay([year, month, day, zone, hour, minute])
            call report(line, 6)

            if (minarg > 5) return
            julianDay = getJulianDay(year, month, day, zone, hour)
            call report(line, 5)
            julianDay = getJulianDay([year, month, day, zone, hour])
            call report(line, 5)

            if (minarg > 4) return
            julianDay = getJulianDay(year, month, day, zone)
            call report(line, 4)
            julianDay = getJulianDay([year, month, day, zone])
            call report(line, 4)

            if (minarg > 3) return
            julianDay = getJulianDay(year, month, day)
            call report(line, 3)
            julianDay = getJulianDay([year, month, day])
            call report(line, 3)

            if (minarg > 2) return
            julianDay = getJulianDay(year, month)
            call report(line, 2)
            julianDay = getJulianDay([year, month])
            call report(line, 2)

            if (minarg > 1) return
            julianDay = getJulianDay(year)
            call report(line, 1)
            julianDay = getJulianDay([year])
            call report(line, 1)

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, narg)
            integer, intent(in) :: line, narg
            real(RKG) :: diff
            diff = abs(julianDay - julianDay_ref)
            assertion = assertion .and. logical(diff <= tol * abs(julianDay_ref) + epsilon(0._RKG), LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "year           ", narg
                write(test%disp%unit,"(*(g0,:,', '))") "year           ", year
                write(test%disp%unit,"(*(g0,:,', '))") "month          ", month
                write(test%disp%unit,"(*(g0,:,', '))") "day            ", day
                write(test%disp%unit,"(*(g0,:,', '))") "zone           ", zone
                write(test%disp%unit,"(*(g0,:,', '))") "hour           ", hour
                write(test%disp%unit,"(*(g0,:,', '))") "minute         ", minute
                write(test%disp%unit,"(*(g0,:,', '))") "second         ", second
                write(test%disp%unit,"(*(g0,:,', '))") "millisecond    ", millisecond
                write(test%disp%unit,"(*(g0,:,', '))") "julianDay_ref  ", julianDay_ref
                write(test%disp%unit,"(*(g0,:,', '))") "julianDay      ", julianDay
                write(test%disp%unit,"(*(g0,:,', '))") "diff           ", diff
                write(test%disp%unit,"(*(g0,:,', '))") "tol            ", tol
                write(test%disp%unit,"(*(g0,:,', '))") "reltol         ", tol * abs(julianDay_ref)
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"@test_getJulianDay(): The Julian day must be computed correctly for the specified date and time with narg = "//getStr(narg), int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_getDateTimeShifted

        use pm_kind, only: RKG => RK, IKG => IK
        integer(IKG), allocatable :: dateTimeShifted(:), dateTimeShifted_ref(:)
        integer(IKG) :: year, month, day, zone, hour, minute, second, millisecond
        !real(RKG), parameter :: ABSTOL = epsilon(0._RKG) * 100
        real(RKG) :: amount
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dateTimeShifted_ref = getDateTime()
        dateTimeShifted = getDateTimeShifted(1._RKG)
        dateTimeShifted_ref(3) = dateTimeShifted_ref(3) + 1_IKG
        assertion = assertion .and. logical(all(dateTimeShifted(1:6) == dateTimeShifted_ref(1:6)), LK) .or. any([dateTimeShifted_ref(1:3), dateTimeShifted_ref(5:6)] == 1_IK)
        call test%assert(assertion, SK_"@test_getDateTimeShifted(): The input time must be shifted correctly. dateTimeShifted = "//getStr(dateTimeShifted)//SK_", dateTimeShifted_ref = "//getStr(dateTimeShifted_ref), int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dateTimeShifted_ref = [integer(IKG) :: +2000, +1, +1, +0, +0, +0, +0, +0]
        amount      = +365._RKG
        year        = +1999_IKG
        month       = +1_IKG
        day         = +1_IKG
        zone        = +0_IKG
        hour        = +0_IKG
        minute      = +0_IKG
        second      = +0_IKG
        millisecond = +0_IKG

        call setAssertion(__LINE__, 1)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dateTimeShifted_ref = [integer(IKG) :: +0, +12, +30, +660, +19, +12, +0, +0]
        amount      = -1.2_RKG
        year        = +1_IKG
        month       = +1_IKG
        day         = +1_IKG
        zone        = +660_IKG
        hour        = +0_IKG
        minute      = +0_IKG
        second      = +0_IKG
        millisecond = +0_IKG

        call setAssertion(__LINE__, 4)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dateTimeShifted_ref = [integer(IKG) :: +0, +1, +3, -660, +8, +0, +0, +0]
        amount      = +2.5_RKG
        year        = -1_IKG
        month       = +12_IKG
        day         = +31_IKG
        zone        = -660_IKG
        hour        = +20_IKG
        minute      = +0_IKG
        second      = +0_IKG
        millisecond = +0_IKG

        call setAssertion(__LINE__, 5)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dateTimeShifted_ref = [integer(IKG) :: +1999, +12, +31, +0, +0, +0, +0, +0]
        amount      = -366._RKG
        year        = +2000_IKG
        month       = +12_IKG
        day         = +31_IKG
        zone        = +0_IKG
        hour        = +0_IKG
        minute      = +0_IKG
        second      = +0_IKG
        millisecond = +0_IKG

        call setAssertion(__LINE__, 3)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dateTimeShifted_ref = [integer(IKG) :: +1998, +11, +21, +660, +8, +21, +35, +847]
        amount      = -100._RKG
        year        = +1999_IKG
        month       = +3_IKG
        day         = +1_IKG
        zone        = +660_IKG
        hour        = +8_IKG
        minute      = +21_IKG
        second      = +35_IKG
        millisecond = +847_IKG

        call setAssertion(__LINE__, 8)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dateTimeShifted_ref = [integer(IKG) :: +2002, +12, +31, -660, +0, +0, +0, +0]
        amount      = +2*365._RKG
        year        = +2000_IKG
        month       = +12_IKG
        day         = +31_IKG
        zone        = -660_IKG
        hour        = +0_IKG
        minute      = +0_IKG
        second      = +0_IKG
        millisecond = +0_IKG

        call setAssertion(__LINE__, 4)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine setAssertion(line, minarg)

            integer, intent(in) :: line, minarg

            if (minarg > 8) return
            dateTimeShifted = getDateTimeShifted(amount, year, month, day, zone, hour, minute, second, millisecond)
            call report(line, 8)
            dateTimeShifted = getDateTimeShifted(amount, [year, month, day, zone, hour, minute, second, millisecond])
            call report(line, 8)

            if (minarg > 7) return
            dateTimeShifted = getDateTimeShifted(amount, year, month, day, zone, hour, minute, second)
            call report(line, 7)
            dateTimeShifted = getDateTimeShifted(amount, [year, month, day, zone, hour, minute, second])
            call report(line, 7)

            if (minarg > 6) return
            dateTimeShifted = getDateTimeShifted(amount, year, month, day, zone, hour, minute)
            call report(line, 6)
            dateTimeShifted = getDateTimeShifted(amount, [year, month, day, zone, hour, minute])
            call report(line, 6)

            if (minarg > 5) return
            dateTimeShifted = getDateTimeShifted(amount, year, month, day, zone, hour)
            call report(line, 5)
            dateTimeShifted = getDateTimeShifted(amount, [year, month, day, zone, hour])
            call report(line, 5)

            if (minarg > 4) return
            dateTimeShifted = getDateTimeShifted(amount, year, month, day, zone)
            call report(line, 4)
            dateTimeShifted = getDateTimeShifted(amount, [year, month, day, zone])
            call report(line, 4)

            if (minarg > 3) return
            dateTimeShifted = getDateTimeShifted(amount, year, month, day)
            call report(line, 3)
            dateTimeShifted = getDateTimeShifted(amount, [year, month, day])
            call report(line, 3)

            if (minarg > 2) return
            dateTimeShifted = getDateTimeShifted(amount, year, month)
            call report(line, 2)
            dateTimeShifted = getDateTimeShifted(amount, [year, month])
            call report(line, 2)

            if (minarg > 1) return
            dateTimeShifted = getDateTimeShifted(amount, year)
            call report(line, 1)
            dateTimeShifted = getDateTimeShifted(amount, [year])
            call report(line, 1)

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, narg)
            integer, intent(in) :: line, narg
            assertion = assertion .and. logical(all(dateTimeShifted == dateTimeShifted_ref), LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "narg               ", narg
                write(test%disp%unit,"(*(g0,:,', '))") "year               ", year
                write(test%disp%unit,"(*(g0,:,', '))") "month              ", month
                write(test%disp%unit,"(*(g0,:,', '))") "day                ", day
                write(test%disp%unit,"(*(g0,:,', '))") "zone               ", zone
                write(test%disp%unit,"(*(g0,:,', '))") "hour               ", hour
                write(test%disp%unit,"(*(g0,:,', '))") "minute             ", minute
                write(test%disp%unit,"(*(g0,:,', '))") "second             ", second
                write(test%disp%unit,"(*(g0,:,', '))") "millisecond        ", millisecond
                write(test%disp%unit,"(*(g0,:,', '))") "dateTimeShifted    ", dateTimeShifted
                write(test%disp%unit,"(*(g0,:,', '))") "dateTimeShifted_ref = ", dateTimeShifted_ref
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"@test_getDateTimeShifted(): The input time must be shifted correctly with narg = "//getStr(narg), int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end  procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_getDateTimeDiff

        use pm_kind, only: RKG => RK, IKG => IK
        integer(IKG), allocatable :: DateTime1(:), DateTime2(:)
        !   \warning
        !   Note that the conversion of time difference to integer days has a roughly single precision accuracy.
        !   As such, the following tests will fail if ABSTOL is set to anything significantly larger than single precision 32-bit accuracy.
        real(RKG), parameter :: ABSTOL = real(epsilon(0.), RKG)
        real(RKG) :: dateTimeDiff, dateTimeDiff_ref, reltol
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dateTimeDiff_ref = 1._RKG
        DateTime1 = getDateTime(1999_IKG, 3_IKG, 1_IKG)
        DateTime2 = getDateTime(1999_IKG, 2_IKG, 28_IKG)
        dateTimeDiff = getDateTimeDiff(DateTime1, DateTime2)
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dateTimeDiff_ref = 2._RKG
        DateTime1 = getDateTime(2000_IKG, 3_IKG, 1_IKG)
        DateTime2 = getDateTime(2000_IKG, 2_IKG, 28_IKG)
        dateTimeDiff = getDateTimeDiff(DateTime1, DateTime2)
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dateTimeDiff_ref = -364._RKG
        DateTime1 = getDateTime(2019_IKG, 1_IKG, 1_IKG)
        DateTime2 = getDateTime(2019_IKG, 12_IKG, 31_IKG)
        dateTimeDiff = getDateTimeDiff(DateTime1, DateTime2)
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dateTimeDiff_ref = -365._RKG
        DateTime1 = getDateTime(2020_IKG, 1_IKG, 1_IKG)
        DateTime2 = getDateTime(2020_IKG, 12_IKG, 31_IKG)
        dateTimeDiff = getDateTimeDiff(DateTime1, DateTime2)
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dateTimeDiff_ref = -1._RKG
        DateTime1 = getDateTime(1_IKG, 1_IKG, 1_IKG, zone = +12_IKG * 60_IKG)
        DateTime2 = getDateTime(1_IKG, 1_IKG, 1_IKG, zone = -12_IKG * 60_IKG)
        dateTimeDiff = getDateTimeDiff(DateTime1, DateTime2)
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dateTimeDiff_ref = +1._RKG
        DateTime1 = getDateTime(-1_IKG, 1_IKG, 1_IKG, zone = -12_IKG * 60_IKG)
        DateTime2 = getDateTime(-1_IKG, 1_IKG, 1_IKG, zone = +12_IKG * 60_IKG)
        dateTimeDiff = getDateTimeDiff(DateTime1, DateTime2)
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dateTimeDiff_ref = +8334.7544704633765_RKG
        DateTime1 = getDateTime(+2022_IKG, +10_IKG, +26_IKG, +0_IKG, +18_IKG, +6_IKG, +26_IKG, +248_IKG)
        DateTime2 = getDateTime(2000_IK)
        dateTimeDiff = getDateTimeDiff(DateTime1, DateTime2)
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        block
            use pm_distUnif, only: getUnifRand
            integer :: i
            do i = 1, 5000
                !   \warning
                !   Note that the conversion of time difference to integer days has a roughly single precision accuracy.
                !   As such, the following tests will fail if ABSTOL is set to anything significantly larger than single precision 32-bit accuracy.
                dateTimeDiff_ref = getUnifRand(-300000._RKG, +300000._RKG)
                DateTime1 = getDateTime()
                DateTime2 = getDateTimeShifted(-dateTimeDiff_ref, DateTime1)
                dateTimeDiff = getDateTimeDiff(DateTime1, DateTime2)
                call report(__LINE__)
            end do
        end block

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line)
            integer, intent(in) :: line
            real(RKG) :: diff
            diff = abs(dateTimeDiff - dateTimeDiff_ref)
            reltol = abs(dateTimeDiff_ref) * ABSTOL + epsilon(0._RKG)
            assertion = assertion .and. logical(diff <= reltol, LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "DateTime1          ", DateTime1
                write(test%disp%unit,"(*(g0,:,', '))") "DateTime2          ", DateTime2
                write(test%disp%unit,"(*(g0,:,', '))") "ABSTOL             ", ABSTOL
                write(test%disp%unit,"(*(g0,:,', '))") "reltol             ", reltol
                write(test%disp%unit,"(*(g0,:,', '))") "diff               ", diff
                write(test%disp%unit,"(*(g0,:,', '))") "dateTimeDiff       ", dateTimeDiff
                write(test%disp%unit,"(*(g0,:,', '))") "dateTimeDiff_ref   ", dateTimeDiff_ref
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"@test_getDateTimeShifted(): The input time must be shifted correctly.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end  procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_getDateTimeUTC

        use pm_distUnif, only: getUnifRand
        use pm_kind, only: RKG => RK, IKG => IK
        integer(IKG), allocatable :: DateTime(:), DateTimeUTC(:), DateTimeUTC_ref(:)
        !!   \warning
        !!   Note that the conversion of time difference to integer days has a roughly single precision accuracy.
        !!   As such, the following tests will fail if ABSTOL is set to anything significantly larger than single precision 32-bit accuracy.
        !real(RKG), parameter :: ABSTOL = real(epsilon(0.), RKG)
        !real(RKG) :: dateTimeDiff, dateTimeDiff_ref, reltol
        integer :: i, ub
        integer(IKG) :: zone
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 10
            DateTime = getDateTime()
            DateTimeUTC = getDateTimeUTC()
            DateTimeUTC_ref = getDateTimeShifted(-DateTime(4) / 1440._RKG, DateTime)
            DateTimeUTC_ref(4) = 0_IKG
            assertion = logical(all(DateTimeUTC(1:6) == DateTimeUTC_ref(1:6)), LK)
            if (assertion) exit
        end do
        call test%assert(assertion, SK_"@test_getDateTimeUTC(): The current UTC time must be computed correctly.", int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 5000

            if (i == 1) then
                zone = -12_IKG * 60_IKG
            elseif (i == 2) then
                zone = +14_IKG * 60_IKG
            else
                zone = getUnifRand(-12_IKG * 60_IKG, +14_IKG * 60_IKG)
            end if

            DateTime = getDateTime()
            DateTime(4) = zone

            ub = 4
            DateTimeUTC = getDateTimeUTC(DateTime(1:ub))
            DateTimeUTC_ref = getDateTimeShifted(-zone / 1440._RKG, DateTime(1:ub))
            DateTimeUTC_ref(ub + 1:) = [DateTimeUTC(ub + 1:)]
            DateTimeUTC_ref(4) = 0_IKG
            call report(__LINE__)
            DateTimeUTC = getDateTimeUTC(DateTime(1), DateTime(2), DateTime(3), DateTime(4))

            ub = 5
            DateTimeUTC = getDateTimeUTC(DateTime(1:ub))
            DateTimeUTC_ref = getDateTimeShifted(-zone / 1440._RKG, DateTime(1:ub))
            DateTimeUTC_ref(ub + 1:) = [DateTimeUTC(ub + 1:)]
            DateTimeUTC_ref(4) = 0_IKG
            call report(__LINE__)
            DateTimeUTC = getDateTimeUTC(DateTime(1), DateTime(2), DateTime(3), DateTime(4), DateTime(5))

            ub = 6
            DateTimeUTC = getDateTimeUTC(DateTime(1:ub))
            DateTimeUTC_ref = getDateTimeShifted(-zone / 1440._RKG, DateTime(1:ub))
            DateTimeUTC_ref(ub + 1:) = [DateTimeUTC(ub + 1:)]
            DateTimeUTC_ref(4) = 0_IKG
            call report(__LINE__)
            DateTimeUTC = getDateTimeUTC(DateTime(1), DateTime(2), DateTime(3), DateTime(4), DateTime(5), DateTime(6))

            ub = 7
            DateTimeUTC = getDateTimeUTC(DateTime(1:ub))
            DateTimeUTC_ref = getDateTimeShifted(-zone / 1440._RKG, DateTime(1:ub))
            DateTimeUTC_ref(ub + 1:) = [DateTimeUTC(ub + 1:)]
            DateTimeUTC_ref(4) = 0_IKG
            call report(__LINE__)
            DateTimeUTC = getDateTimeUTC(DateTime(1), DateTime(2), DateTime(3), DateTime(4), DateTime(5), DateTime(7))

            ub = 8
            DateTimeUTC = getDateTimeUTC(DateTime(1:ub))
            DateTimeUTC_ref = getDateTimeShifted(-zone / 1440._RKG, DateTime(1:ub))
            DateTimeUTC_ref(ub + 1:) = [DateTimeUTC(ub + 1:)]
            DateTimeUTC_ref(4) = 0_IKG
            call report(__LINE__)
            DateTimeUTC = getDateTimeUTC(DateTime(1), DateTime(2), DateTime(3), DateTime(4), DateTime(5), DateTime(7), DateTime(8))

        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line)
            integer, intent(in) :: line
            assertion = assertion .and. logical(all(DateTimeUTC == DateTimeUTC_ref), LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "DateTimeUTC    ", DateTimeUTC
                write(test%disp%unit,"(*(g0,:,', '))") "DateTimeUTC_ref", DateTimeUTC_ref
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"@test_getDateTimeUTC(): The input time must be shifted correctly.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end  procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_getDateTimeNewZone

        use pm_distUnif, only: getUnifRand
        use pm_kind, only: RKG => RK, IKG => IK
        integer(IKG), allocatable :: DateTime(:), DateTimeNewZone(:), DateTimeNewZone_ref(:)
        integer(IKG) :: zone, newzone
        real(RKG) :: zoneshift
        integer :: i, ub
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 10
            DateTime = getDateTime()
            zone = DateTime(4)
            newzone = getUnifRand(-12_IKG * 60_IKG, +14_IKG * 60_IKG)
            DateTimeNewZone = getDateTimeNewZone(newzone)
            zoneshift = -(DateTime(4) - newzone) / 1440._RKG
            DateTimeNewZone_ref = getDateTimeShifted(zoneshift, DateTime)
            DateTimeNewZone_ref(6:8) = DateTimeNewZone(6:8)
            DateTimeNewZone_ref(4) = newzone
            call report(__LINE__, i /= 10)
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        newzone = -300_IK
        DateTime = [integer(IKG) :: 1999_IK, 3_IK, 1_IK, +660_IK, 8_IK, 21_IK, 35_IK, 847_IK]
        DateTimeNewZone_ref = [integer(IKG) :: +1999, +2, +28, -300, +16, +21, +35, +847]
        DateTimeNewZone = getDateTimeNewZone(newzone, DateTime)
        zone = DateTime(4)
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        newzone = -300_IK
        DateTime = [integer(IKG) :: 2000_IK, 3_IK, 1_IK, +660_IK, 8_IK, 21_IK, 35_IK, 847_IK]
        DateTimeNewZone_ref = [integer(IKG) :: +2000, +2, +29, -300, +16, +21, +35, +847]
        DateTimeNewZone = getDateTimeNewZone(newzone, DateTime)
        zone = DateTime(4)
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        newzone = +300_IK
        DateTime = [integer(IKG) :: 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK, 35_IK, 847_IK]
        DateTimeNewZone_ref = [integer(IKG) :: +2001, +1, +1, +300, +10, +21, +35, +847]
        DateTimeNewZone = getDateTimeNewZone(newzone, DateTime)
        zone = DateTime(4)
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        newzone = +300_IK
        DateTime = [integer(IKG) :: 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK, 35_IK]
        DateTimeNewZone_ref = [integer(IKG) :: +2001, +1, +1, +300, +10, +21, +35, +0]
        DateTimeNewZone = getDateTimeNewZone(newzone, DateTime)
        zone = DateTime(4)
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        newzone = +313_IK
        DateTime = [integer(IKG) :: 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK]
        DateTimeNewZone_ref = [integer(IKG) :: +2001, +1, +1, +313, +10, +34, +0, +0]
        DateTimeNewZone = getDateTimeNewZone(newzone, DateTime)
        zone = DateTime(4)
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        newzone = +300_IK
        DateTime = [integer(IKG) :: 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK, 21_IK]
        DateTimeNewZone_ref = [integer(IKG) :: +2001, +1, +1, +300, +10, +21, +0, +0]
        DateTimeNewZone = getDateTimeNewZone(newzone, DateTime)
        zone = DateTime(4)
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        newzone = +660_IK
        DateTime = [integer(IKG) :: 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK]
        DateTimeNewZone_ref = [integer(IKG) :: +2001, +1, +1, +660, +16, +0, +0, +0]
        DateTimeNewZone = getDateTimeNewZone(newzone, DateTime)
        zone = DateTime(4)
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        newzone = +0_IK
        DateTime = [integer(IKG) :: 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK]
        DateTimeNewZone_ref = [integer(IKG) :: +2001, +1, +1, +0, +5, +0, +0, +0]
        DateTimeNewZone = getDateTimeNewZone(newzone, DateTime)
        zone = DateTime(4)
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        newzone = -660_IK
        DateTime = [integer(IKG) :: 2000_IK, 12_IK, 31_IK, -660_IK, 18_IK]
        DateTimeNewZone_ref = [integer(IKG) :: +2000, +12, +31, -660, +18, +0, +0, +0]
        DateTimeNewZone = getDateTimeNewZone(newzone, DateTime)
        zone = DateTime(4)
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 5000

            if (i == 1) then
                zone = -12_IKG * 60_IKG
            elseif (i == 2) then
                zone = +14_IKG * 60_IKG
            else
                zone = getUnifRand(-12_IKG * 60_IKG, +14_IKG * 60_IKG)
            end if

            newzone = getUnifRand(-12_IKG * 60_IKG, +14_IKG * 60_IKG)
            zoneshift = -(zone - newzone) / 1440._RKG
            DateTime = getDateTime()
            DateTime(4) = zone

            ub = 4
            DateTimeNewZone = getDateTimeNewZone(newzone, DateTime(1:ub))
            DateTimeNewZone_ref = getDateTimeShifted(zoneshift, DateTime(1:ub))
            DateTimeNewZone_ref(ub + 1:) = [DateTimeNewZone(ub + 1:)]
            DateTimeNewZone_ref(4) = newzone
            call report(__LINE__)
            DateTimeNewZone = getDateTimeNewZone(newzone, DateTime(1), DateTime(2), DateTime(3), DateTime(4))

            ub = 5
            DateTimeNewZone = getDateTimeNewZone(newzone, DateTime(1:ub))
            DateTimeNewZone_ref = getDateTimeShifted(zoneshift, DateTime(1:ub))
            DateTimeNewZone_ref(ub + 1:) = [DateTimeNewZone(ub + 1:)]
            DateTimeNewZone_ref(4) = newzone
            call report(__LINE__)
            DateTimeNewZone = getDateTimeNewZone(newzone, DateTime(1), DateTime(2), DateTime(3), DateTime(4), DateTime(5))

            ub = 6
            DateTimeNewZone = getDateTimeNewZone(newzone, DateTime(1:ub))
            DateTimeNewZone_ref = getDateTimeShifted(zoneshift, DateTime(1:ub))
            DateTimeNewZone_ref(ub + 1:) = [DateTimeNewZone(ub + 1:)]
            DateTimeNewZone_ref(4) = newzone
            call report(__LINE__)
            DateTimeNewZone = getDateTimeNewZone(newzone, DateTime(1), DateTime(2), DateTime(3), DateTime(4), DateTime(5), DateTime(6))

            ub = 7
            DateTimeNewZone = getDateTimeNewZone(newzone, DateTime(1:ub))
            DateTimeNewZone_ref = getDateTimeShifted(zoneshift, DateTime(1:ub))
            DateTimeNewZone_ref(ub + 1:) = [DateTimeNewZone(ub + 1:)]
            DateTimeNewZone_ref(4) = newzone
            call report(__LINE__)
            DateTimeNewZone = getDateTimeNewZone(newzone, DateTime(1), DateTime(2), DateTime(3), DateTime(4), DateTime(5), DateTime(7))

            ub = 8
            DateTimeNewZone = getDateTimeNewZone(newzone, DateTime(1:ub))
            DateTimeNewZone_ref = getDateTimeShifted(zoneshift, DateTime(1:ub))
            DateTimeNewZone_ref(ub + 1:) = [DateTimeNewZone(ub + 1:)]
            DateTimeNewZone_ref(4) = newzone
            call report(__LINE__)
            DateTimeNewZone = getDateTimeNewZone(newzone, DateTime(1), DateTime(2), DateTime(3), DateTime(4), DateTime(5), DateTime(7), DateTime(8))

        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, continued)
            integer, intent(in) :: line
            logical, intent(in), optional :: continued
            assertion = assertion .and. logical(all(DateTimeNewZone == DateTimeNewZone_ref), LK)
            if (present(continued)) then
                assertion = assertion .or. logical(continued, LK)
            end if
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "zone               ", zone
                write(test%disp%unit,"(*(g0,:,', '))") "newzone            ", newzone
                write(test%disp%unit,"(*(g0,:,', '))") "DateTime           ", DateTime
                write(test%disp%unit,"(*(g0,:,', '))") "DateTimeNewZone    ", DateTimeNewZone
                write(test%disp%unit,"(*(g0,:,', '))") "DateTimeNewZone_ref", DateTimeNewZone_ref
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"@test_getDateTimeNewZone(): The time in the new zone must be computed correctly.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end  procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_getDateTime_1

        use pm_distUnif, only: getUnifRand
        use pm_kind, only: RKG => RK, IKG => IK
        integer(IKG), allocatable :: Values(:), values_ref(:)
        !integer(IKG) :: year, month, day, zone, hour, minute, second, millisecond
        real(RKG) :: julianDay
        integer :: i
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! check the condition `if (Values(5) == 24_IKG) then` in the code.
        values_ref = [integer(IKG) :: 1, 1, 1, -60, 23, 0, 0, 0]
        julianDay = getJulianDay(values_ref)
        values_ref = [integer(IKG) :: 1, 1, 2, 0, 0, 0, 0, 0]
        Values = getDateTimeNewZone(values_ref(4), getDateTime(julianDay))
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 500

            values_ref = getDateTime()
            julianDay = getJulianDay(values_ref)
            Values = getDateTimeNewZone(values_ref(4), getDateTime(julianDay))
            call report(__LINE__)

            Values = getDateTime(julianDay, values_ref(4))
            call report(__LINE__, values_ref(4))

            Values = getDateTime(julianDay, values_ref(4))
            call report(__LINE__, values_ref(4))

            values_ref(4) = getUnifRand(-12_IKG * 60_IKG, +14_IKG * 60_IKG)
            julianDay = getJulianDay(values_ref)
            Values = getDateTime(julianDay, values_ref(4))
            call report(__LINE__, values_ref(4))

            !year        = +values_ref(1)
            !month       = +values_ref(2)
            !day         = +values_ref(3)
            !zone        = +values_ref(4)
            !hour        = +values_ref(5)
            !minute      = +values_ref(6)
            !second      = +values_ref(7)
            !millisecond = +values_ref(8)

            Values = getDateTime(values_ref(1), values_ref(2), values_ref(3), values_ref(4), values_ref(5), values_ref(6), values_ref(7), values_ref(8))
            call report(__LINE__)

            Values = getDateTime(values_ref(1), values_ref(2), values_ref(3), values_ref(4), values_ref(5), values_ref(6), values_ref(7))
            values_ref(8) = 0_IKG
            call report(__LINE__)

            Values = getDateTime(values_ref(1), values_ref(2), values_ref(3), values_ref(4), values_ref(5), values_ref(6))
            values_ref(7) = 0_IKG
            call report(__LINE__)

            Values = getDateTime(values_ref(1), values_ref(2), values_ref(3), values_ref(4), values_ref(5))
            values_ref(6) = 0_IKG
            call report(__LINE__)

            Values = getDateTime(values_ref(1), values_ref(2), values_ref(3), values_ref(4))
            values_ref(5) = 0_IKG
            call report(__LINE__)

            Values = getDateTime(values_ref(1), values_ref(2), values_ref(3))
            values_ref(4) = 0_IKG
            call report(__LINE__)

            Values = getDateTime(values_ref(1), values_ref(2))
            values_ref(3) = 1_IKG
            call report(__LINE__)

            Values = getDateTime(values_ref(1))
            values_ref(2) = 1_IKG
            call report(__LINE__)

        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, zone)
            integer, intent(in) :: line
            integer, intent(in), optional :: zone
            assertion = assertion .and. logical(all(Values == values_ref), LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                if (present(zone)) then
                write(test%disp%unit,"(*(g0,:,', '))") "zone       ", zone
                end if
                write(test%disp%unit,"(*(g0,:,', '))") "Values     ", Values
                write(test%disp%unit,"(*(g0,:,', '))") "values_ref ", values_ref
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"@test_getDateTimeNewZone(): The time in the new zone must be computed correctly.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end  procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_getDateTime_2

        use pm_distUnif, only: getUnifRand
        use pm_kind, only: RKG => RK, IKG => IK, SKG => SK
        character(:,SKG), allocatable :: string, string_ref
        character(2,SKG), parameter :: SPECIFIER(*) =   [ "%a", "%A", "%b", "%B", "%c", "%C", "%d", "%D", "%e", "%F", "%f", "%g", "%G" &
                                                        , "%h", "%H", "%I", "%j", "%m", "%M", "%n", "%p", "%r", "%R", "%S", "%t", "%T" &
                                                        , "%u", "%U", "%V", "%w", "%W", "%x", "%X", "%y", "%Y", "%z", "%Z", "%%" &
                                                        ]
        integer(IKG), allocatable :: Values(:)
        integer :: i, j, k

        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_""
            string = getDateTime(format = SK_"")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): An empty format must yield an empty string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"Today is "//trim(WEEKDAY_NAME(getWeekDay())(1:3))//SKG_"."
            string = getDateTime(format = SK_"Today is %a.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %a format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"Today is "//trim(WEEKDAY_NAME(getWeekDay()))//SKG_"."
            string = getDateTime(format = SK_"Today is %A.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %A format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"This month is "//trim(MONTH_NAME(getMonth())(1:3))//SKG_"."
            string = getDateTime(format = SK_"This month is %b.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %b format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"This month is "//trim(MONTH_NAME(getMonth()))//SKG_"."
            string = getDateTime(format = SK_"This month is %B.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %B format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The date is "//trim(WEEKDAY_NAME(getWeekDay())(1:3))//SKG_" "//trim(MONTH_NAME(getMonth())(1:3))//SKG_" "//trim(getStr(getDay(), format = SK_"(I0.2)"))//SKG_" "// & ! LCOV_EXCL_LINE
            trim(getStr(getHour(), format = SK_"(I0.2)"))//SKG_":"//trim(getStr(getMinute(), format = SK_"(I0.2)"))//SKG_":"//trim(getStr(getSecond(), format = SK_"(I0.2)"))//SKG_" "//getStr(getYear())//SKG_"."
            string = getDateTime(format = SK_"The date is %c.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %c format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"This century is "//getStr(floor(getYear() / 100., IKG), signed = .true._LK)//SKG_"."
            string = getDateTime(format = SK_"This century is %C.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %C format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"Day of the month is "//getStr(getDay(), format = SK_"(I0.2)")//SKG_"."
            string = getDateTime(format = SK_"Day of the month is %d.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %d format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The short date is "//getDateTime(SKG_"%m/%d/%y")//SKG_"."
            string = getDateTime(format = SK_"The short date is %D.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %D format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"Day of the month is "//adjustr(getStr(getDay(), format = SK_"(I2)"))//SKG_"."
            string = getDateTime(format = SK_"Day of the month is %e.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %e format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The short date is "//getDateTime(SKG_"%Y-%m-%d")//SKG_"."
            string = getDateTime(format = SK_"The short date is %F.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %F format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            Values = getDateTime()
            string_ref = SKG_"The millisecond of the second is "//getStr(Values(8), format = SK_"(I3.3)")//SKG_"."
            string = getDateTime(format = SK_"The millisecond of the second is %f.", Values = Values)
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %f format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The week year is "//getStr(mod(abs(getWeekYear()), 100_IKG), format = SK_"(I0.2)")//SKG_"."
            string = getDateTime(format = SK_"The week year is %g.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %g format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The week year is "//getStr(getWeekYear())//SKG_"."
            string = getDateTime(format = SK_"The week year is %G.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %G format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"This month is "//trim(MONTH_NAME(getMonth())(1:3))//SKG_"."
            string = getDateTime(format = SK_"This month is %h.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %h format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The hour is "//trim(getStr(getHour(), format = SK_"(I0.2)"))//SKG_"."
            string = getDateTime(format = SK_"The hour is %H.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %H format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The hour is "//trim(getStr(getHour12(), format = SK_"(I0.2)"))//SKG_"."
            string = getDateTime(format = SK_"The hour is %I.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %I format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The day of year is "//trim(getStr(getOrdinalDay(), format = SK_"(I3.3)"))//SKG_"."
            string = getDateTime(format = SK_"The day of year is %j.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %j format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The month is "//trim(getStr(getMonth(), format = SK_"(I0.2)"))//SKG_"."
            string = getDateTime(format = SK_"The month is %m.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %m format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The minute is "//trim(getStr(getMinute(), format = SK_"(I0.2)"))//SKG_"."
            string = getDateTime(format = SK_"The minute is %M.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %M format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The newline character is "//achar(10, SKG)//SKG_"."
            string = getDateTime(format = SK_"The newline character is %n.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %n format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The designation is "//merge(SKG_"AM", SKG_"PM", isMorning())//SKG_"."
            string = getDateTime(format = SK_"The designation is %p.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %p format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The 12-hour clock time is "//getStr(getHour12(), format = SK_"(I0.2)")//SKG_":"//getStr(getMinute(), format = SK_"(I0.2)")//SKG_":"//getStr(getSecond(), format = SK_"(I0.2)")//SKG_" "//merge(SKG_"am", SKG_"pm", isMorning())//SKG_"."
            string = getDateTime(format = SK_"The 12-hour clock time is %r.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %r format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The 24-hour time is "//getDateTime(SKG_"%H:%M")//SKG_"."
            string = getDateTime(format = SK_"The 24-hour time is %R.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %R format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The second of time is "//getDateTime(getStr(getSecond(), format = SK_"(I0.2)"))//SKG_"."
            string = getDateTime(format = SK_"The second of time is %S.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %S format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The Horizontal-tab character is "//achar(9, SKG)//SKG_"."
            string = getDateTime(format = SK_"The Horizontal-tab character is %t.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %t format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The ISO 8601 time is "//getDateTime(SKG_"%H:%M:%S")//SKG_"."
            string = getDateTime(format = SK_"The ISO 8601 time is %T.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %T format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The ISO 8601 weekday is "//getStr(getWeekDayISO())//SKG_"."
            string = getDateTime(format = SK_"The ISO 8601 weekday is %u.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %u format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The ISO 8601 week number is "//getStr(getWeekNumber(), format = SK_"(I0.2)")//SKG_"."
            string = getDateTime(format = SK_"The ISO 8601 week number is %V.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %V format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The Weekday is "//getStr(getWeekDay())//SKG_"."
            string = getDateTime(format = SK_"The Weekday is %w.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %w format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = getDateTime(SKG_"The current date is %D.")
            string = getDateTime(format = SK_"The current date is %x.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %x format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = getDateTime(SKG_"The current time is %T.")
            string = getDateTime(format = SK_"The current time is %X.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %X format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The year is "//getStr(mod(abs(Values(1)), 100_IKG), format = SK_"(I0.2)")//SKG_"."
            string = getDateTime(format = SK_"The year is %y.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %y format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The year is "//getStr(Values(1), format = SK_"(I0.2)")//SKG_"."
            string = getDateTime(format = SK_"The year is %Y.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %Y format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The zone is "//getStr(Values(4), format = SK_"(sp,I0.4)")//SKG_"."
            string = getDateTime(format = SK_"The zone is %z.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %z format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The zone is "//getStr(Values(4), format = SK_"(sp,I0.4)")//SKG_"."
            string = getDateTime(format = SK_"The zone is %z.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %z format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, 10
            string_ref = SKG_"The percentage is %."
            string = getDateTime(format = SK_"The percentage is %%.")
            if (string == string_ref) exit
        end do
        call report(__LINE__, SK_"@test_getDateTime_2(): The %% format must yield a correctly formatted date and time string.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 500

            do j = 1, size(SPECIFIER)
                Values = getDateTime(julianDay = getUnifRand(-300000._RKG, +300000._RKG))
                string_ref = getDateTime_ref(SPECIFIER(j), Values)
                string = getDateTime(SPECIFIER(j), Values)
                call report(__LINE__)
            end do

            Values = getDateTime(julianDay = getUnifRand(-300000._RKG, +300000._RKG))
            j = getUnifRand(1, size(SPECIFIER))
            k = getUnifRand(j, size(SPECIFIER))
            string_ref = getDateTime_ref(getStr(SPECIFIER(j:k)), Values)
            string = getDateTime(getStr(SPECIFIER(j:k)), Values)
            call report(__LINE__)

            Values(1) = getUnifRand(300000_IKG, 300000_IKG)
            string_ref = getDateTime_ref(SKG_"%C", Values)
            string = getDateTime(SKG_"%C", Values)
            call report(__LINE__)

        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function getDateTime_ref(format, Values) result(string)

#define     RESIZE_STRING(LENSEG) \
            eposnew = epos + LENSEG; \
            if (eposnew > lenString) then; \
                if (allocated(tempstr)) deallocate(tempstr); \
                allocate(character(eposnew,SKG) :: tempstr); \
                tempstr(1:lenString) = string; \
                call move_alloc(tempstr, string); \
                lenString = eposnew; \
            end if;

            !>  \warning
            !>  The output of getStr() in this procedure is of kind \SK which is incompatible with any value of SKG /= SK.
            !>  For now, this is not an issue since both kinds point to the default character kind.
            !>  This will however become an issue once Fortran standard and compilers support non-default date and time characters.
            use pm_val2str, only: getStr
            integer(IK), intent(in), contiguous :: Values(:)
            character(*,SKG), intent(in) :: format
            character(:,SKG), allocatable :: string
            character(:,SKG), allocatable   :: tempstr, abbr
            character(9,SKG)                :: workspace9
            integer(IK)                     :: lenString, lenFormat, i, epos, eposnew, lenSeg
            integer(IKG)                    :: century
            allocate(character(127,SKG)     :: string)
            lenFormat = len(format, IK)
            lenString = 127_IK
            eposnew = 0_IK ! the last touched (end) position in the string
            i = 0_IK
            do
                i = i + 1_IK
                if (i > lenFormat) exit
                epos = eposnew
                if (format(i:i) == SKG_"%") then
                    if (i == lenFormat) exit
                    i = i + 1_IK
                    if (format(i:i) == SKG_"a") then ! Abbreviated weekday name *
                        lenSeg = 3_IK
                        RESIZE_STRING(lenSeg) ! fpp sets eposnew, and resizes string.
                        string(epos + 1 : eposnew) = WEEKDAY_NAME_ISO(getWeekDayISO(Values(1), Values(2), Values(3)))(1:lenSeg)
                    elseif (format(i:i) == SKG_"A") then ! Full weekday name *
                        workspace9 = WEEKDAY_NAME_ISO(getWeekDayISO(Values(1), Values(2), Values(3)))
                        lenSeg = len_trim(workspace9, IKG)
                        RESIZE_STRING(lenSeg) ! fpp sets eposnew, and resizes string.
                        string(epos + 1 : eposnew) = workspace9(1:lenSeg)
                    elseif (format(i:i) == SKG_"b" .or. format(i:i) == SKG_"h") then ! Abbreviated month name * .or. Abbreviated month name * (same as %b)
                        lenSeg = 3_IK
                        RESIZE_STRING(lenSeg) ! fpp sets eposnew, and resizes string.
                        string(epos + 1 : eposnew) = MONTH_NAME(Values(2))(1:lenSeg)
                    elseif (format(i:i) == SKG_"B") then ! Full month name *
                        lenSeg = len_trim(MONTH_NAME(Values(2)), IKG)
                        RESIZE_STRING(lenSeg) ! fpp sets eposnew, and resizes string.
                        string(epos + 1 : eposnew) = MONTH_NAME(Values(2))(1:lenSeg)
                    elseif (format(i:i) == SKG_"c") then ! Date and time representation *
                        if (Values(1) > 0_IKG) then
                            RESIZE_STRING(24_IK) ! fpp sets eposnew, and resizes string.
                        else
                            RESIZE_STRING(25_IK) ! fpp sets eposnew, and resizes string.
                        end if
                        string(epos + 1 : eposnew) =    WEEKDAY_NAME_ISO(getWeekDayISO(Values(1), Values(2), Values(3)))(1:3)//SKG_" "// & ! LCOV_EXCL_LINE
                                                        MONTH_NAME(Values(2))(1:3)//SKG_" "// & ! LCOV_EXCL_LINE
                                                        getStr(Values(3), length = 2_IK, format = "(1I0.2)")//SKG_" "// & ! LCOV_EXCL_LINE
                                                        getStr(Values(5), length = 2_IK, format = "(1I0.2)")//SKG_":"// & ! LCOV_EXCL_LINE
                                                        getStr(Values(6), length = 2_IK, format = "(1I0.2)")//SKG_":"// & ! LCOV_EXCL_LINE
                                                        getStr(Values(7), length = 2_IK, format = "(1I0.2)")//SKG_" "// & ! LCOV_EXCL_LINE
                                                        getStr(Values(1))
                    elseif (format(i:i) == SKG_"C") then ! Year divided by 100 and truncated to integer (00-99). \warning it works for years up to 9 digits.
                        century = floor(Values(1) / 100., IKG)
                        if (abs(century) < 100_IKG) then
                            !write(workspace9(1:3), "(sp,I0.2)") century
                            RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                            write(string(epos + 1 : eposnew), "(sp,I0.2)") century
                        else
                            workspace9 = getStr(century)
                            lenSeg = len_trim(workspace9, IKG)
                            RESIZE_STRING(lenSeg) ! fpp sets eposnew, and resizes string.
                            string(epos + 1 : eposnew) = workspace9(1:lenSeg)
                        end if
                    elseif (format(i:i) == SKG_"d") then ! Day of the month, zero-padded (01-31).
                        RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                        write(string(epos + 1 : eposnew), "(I0.2)") Values(3)
                    elseif (format(i:i) == SKG_"D") then ! Short MM/DD/YY date, equivalent to %m/%d/%y
                        RESIZE_STRING(8_IK) ! fpp sets eposnew, and resizes string.
                        write(string(epos + 1 : eposnew), "(I0.2,'/',I0.2,'/',I0.2)") Values(2:3), mod(abs(Values(1)), 100_IKG) ! last two digits of year.
                    elseif (format(i:i) == SKG_"e") then ! Day of the month, zero-padded (01-31).
                        RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                        write(string(epos + 1 : eposnew), "(I2)") Values(3)
                    elseif (format(i:i) == SKG_"f") then ! millisecond padded with leading zeros
                        RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                        write(string(epos + 1 : eposnew), "(I0.3)") Values(8)
                    elseif (format(i:i) == SKG_"F") then ! Short YYYY-MM-DD date, equivalent to %Y-%m-%d
                        if (Values(1) > 0_IKG) then
                            RESIZE_STRING(10_IK) ! fpp sets eposnew, and resizes string.
                        else
                            RESIZE_STRING(11_IK) ! fpp sets eposnew, and resizes string.
                        end if
                        write(string(epos + 1 : eposnew), "(I0.2,'-',I0.2,'-',I0.2)") Values(1:3)
                    elseif (format(i:i) == SKG_"g") then ! Week-based year, last two digits (00-99).
                        RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                        write(string(epos + 1 : eposnew), "(I0.2)") mod(abs(getWeekYear(Values(1:3))), 100_IKG) ! last two digits of the week year.
                    elseif (format(i:i) == SKG_"G") then ! Week-based year, full week year, possibly negative.
                        workspace9 = getStr(getWeekYear(Values(1:3))) ! WeekDate(1)) ! full week year, possibly negative.
                        lenSeg = len_trim(workspace9, IKG)
                        RESIZE_STRING(lenSeg) ! fpp sets eposnew, and resizes string.
                        string(epos + 1 : eposnew) = workspace9(1:lenSeg)
                    elseif (format(i:i) == SKG_"H") then ! Hour in 24h format (00-23)
                        RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                        write(string(epos + 1 : eposnew), "(I0.2)") Values(5)
                    elseif (format(i:i) == SKG_"I") then ! Hour in 12h format (01-12)
                        RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                        write(string(epos + 1 : eposnew), "(I0.2)") getHour12(Values(5))
                    elseif (format(i:i) == SKG_"j") then ! Day of the year (001-366)
                        RESIZE_STRING(3_IK) ! fpp sets eposnew, and resizes string.
                        write(string(epos + 1 : eposnew), "(I0.3)") getOrdinalDay(Values(1:3))
                    elseif (format(i:i) == SKG_"m") then ! Month as a decimal number (01-12)
                        RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                        write(string(epos + 1 : eposnew), "(I0.2)") Values(2)
                    elseif (format(i:i) == SKG_"M") then ! Minute (00-59)
                        RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                        write(string(epos + 1 : eposnew), "(I0.2)") Values(6)
                    elseif (format(i:i) == SKG_"n") then ! New-line character ('\n')
                        RESIZE_STRING(1_IK) ! fpp sets eposnew, and resizes string.
                        string(epos + 1 : eposnew) = achar(10, SKG) ! new_line(SKG_"a")
                    elseif (format(i:i) == SKG_"p") then ! New-line character ('\n')
                        RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                        if (Values(5) < 12_IKG) then
                            string(epos + 1 : eposnew) = SKG_"AM"
                        else
                            string(epos + 1 : eposnew) = SKG_"PM"
                        end if
                    elseif (format(i:i) == SKG_"r") then ! 12-hour clock time *
                        RESIZE_STRING(11_IK) ! fpp sets eposnew, and resizes string.
                        if (Values(5) < 12_IKG) then
                            write(string(epos + 1 : eposnew), "(I0.2,':',I0.2,':',I0.2,' am')") getHour12(Values(5)), Values(6:7)
                        else
                            write(string(epos + 1 : eposnew), "(I0.2,':',I0.2,':',I0.2,' pm')") getHour12(Values(5)), Values(6:7)
                        end if
                    elseif (format(i:i) == SKG_"R") then ! 24-hour HH:MM time, equivalent to %H:%M
                        RESIZE_STRING(5_IK) ! fpp sets eposnew, and resizes string.
                        write(string(epos + 1 : eposnew), "(I0.2,':',I0.2)") Values(5:6)
                    elseif (format(i:i) == SKG_"S") then ! Second (00-59)
                        RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                        write(string(epos + 1 : eposnew), "(I0.2)") Values(7)
                    elseif (format(i:i) == SKG_"t") then ! Horizontal-tab character ('\t')
                        RESIZE_STRING(1_IK) ! fpp sets eposnew, and resizes string.
                        string(epos + 1 : eposnew) = achar(9, SKG)
                    elseif (format(i:i) == SKG_"T") then ! ISO 8601 time format (HH:MM:SS), equivalent to %H:%M:%S
                        RESIZE_STRING(8_IK) ! fpp sets eposnew, and resizes string.
                        write(string(epos + 1 : eposnew), "(I0.2,':',I0.2,':',I0.2)") Values(5:7)
                    elseif (format(i:i) == SKG_"u") then ! ISO 8601 weekday as number with Monday as 1 (1-7)
                        RESIZE_STRING(1_IK) ! fpp sets eposnew, and resizes string.
                        write(string(epos + 1 : eposnew), "(I1)") getWeekDayISO(Values(1), Values(2), Values(3))
                    elseif (format(i:i) == SKG_"V") then ! ISO 8601 week number (01-53)
                        RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                        write(string(epos + 1 : eposnew), "(I0.2)") getWeekNumber(Values(1), Values(2), Values(3))
                    elseif (format(i:i) == SKG_"w") then ! Weekday as a decimal number with Sunday as 0 (0-6)
                        RESIZE_STRING(1_IK) ! fpp sets eposnew, and resizes string.
                        write(string(epos + 1 : eposnew), "(I1)") getWeekDay(Values(1), Values(2), Values(3))
                    elseif (format(i:i) == SKG_"x") then ! Date representation *
                        RESIZE_STRING(8_IK) ! fpp sets eposnew, and resizes string.
                        write(string(epos + 1 : eposnew), "(I0.2,'/',I0.2,'/',I0.2)") Values(2), Values(3), mod(abs(Values(1)), 100_IKG) ! last two digits of year.
                    elseif (format(i:i) == SKG_"X") then ! Time representation *
                        RESIZE_STRING(8_IK) ! fpp sets eposnew, and resizes string.
                        write(string(epos + 1 : eposnew), "(I0.2,':',I0.2,':',I0.2)") Values(5:7)
                    elseif (format(i:i) == SKG_"y") then ! Year, last two digits (00-99)
                        RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                        write(string(epos + 1 : eposnew), "(I0.2)") mod(abs(Values(1)), 100_IKG) ! last two digits of year.
                    elseif (format(i:i) == SKG_"Y") then ! Year, last two digits (00-99)
                        workspace9 = getStr(Values(1))
                        lenSeg = len_trim(workspace9, IKG)
                        RESIZE_STRING(lenSeg) ! fpp sets eposnew, and resizes string.
                        string(epos + 1 : eposnew) = workspace9(1:lenSeg)
                    elseif (format(i:i) == SKG_"z") then ! ISO 8601 offset from UTC in timezone in units of minutes
                        RESIZE_STRING(5_IK) ! fpp sets eposnew, and resizes string.
                        write(string(epos + 1 : eposnew), "(sp,I0.4)") Values(4)
                    elseif (format(i:i) == SKG_"Z") then ! Timezone name or abbreviation.
                        abbr = getZoneAbbr(Values(4))
                        RESIZE_STRING(len(abbr, IK))
                        string(epos + 1 : eposnew) = abbr
                    elseif (format(i:i) == SKG_"%") then ! add percentage.
                        RESIZE_STRING(1_IK) ! fpp sets eposnew, and resizes string.
                        string(epos + 1 : eposnew) = SKG_"%"
                    else ! Unrecognized format.
                        RESIZE_STRING(2_IK) ! fpp sets eposnew, and resizes string.
                        string(epos + 1 : eposnew) = format(i - 1 : i)
                    end if
                else ! normal characters.
                    RESIZE_STRING(1_IK) ! fpp sets eposnew, and resizes string.
                    string(epos + 1 : eposnew) = format(i : i)
                end if
            end do
            if (lenString > eposnew) then
                tempstr = string(1:eposnew)
                call move_alloc(tempstr, string)
            end if

#undef      RESIZE_STRING

        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, msg)
            use pm_option, only: getOption
            integer, intent(in) :: line
            character(*,SKG), intent(in), optional :: msg
            character(:,SKG), allocatable :: msg_def
            assertion = assertion .and. logical(string == string_ref, LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "string             ", string
                write(test%disp%unit,"(*(g0,:,', '))") "string_ref         ", string_ref
                write(test%disp%unit,"(*(g0,:,', '))") "len(string, IK)    ", len(string, IK)
                write(test%disp%unit,"(*(g0,:,', '))") "len(string_ref, IK)", len(string_ref, IK)
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            if (present(msg)) then
                msg_def = msg
            else
                msg_def = SKG_"@test_getDateTimeNewZone(): The output datetime string must have the correct format."
            end if
            call test%assert(assertion, msg_def, int(line, IK))
            assertion = assertion .and. logical(len(string, IK) == len(string_ref, IK), LK)
            call test%assert(assertion, SK_"@test_getDateTimeNewZone(): The output datetime string must have the correct length.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end  procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_getWeekDate

        use pm_distUnif, only: getUnifRand
        use pm_kind, only: SKG => SK, IKG => IK
        integer(IKG), allocatable :: Values(:), WeekDate(:), WeekDate_ref(:)
        integer :: i
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Values = int([1977_IK, 1_IK, 1_IK], IKG)
        WeekDate_ref = int([+1976, +53, +6], IKG)
        call runTests(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Values = int([1977_IK, 1_IK, 2_IK], IKG)
        WeekDate_ref = int([+1976, +53, +7], IKG)
        call runTests(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Values = int([1977_IK, 12_IK, 31_IK], IKG)
        WeekDate_ref = int([+1977, +52, +6], IKG)
        call runTests(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Values = int([1978_IK, 1_IK, 1_IK], IKG)
        WeekDate_ref = int([+1977, +52, +7], IKG)
        call runTests(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Values = int([1978_IK, 1_IK, 2_IK], IKG)
        WeekDate_ref = int([+1978, +1, +1], IKG)
        call runTests(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Values = int([1978_IK, 12_IK, 31_IK], IKG)
        WeekDate_ref = int([+1978, +52, +7], IKG)
        call runTests(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Values = int([1979_IK, 1_IK, 1_IK], IKG)
        WeekDate_ref = int([+1979, +1, +1], IKG)
        call runTests(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Values = int([1979_IK, 12_IK, 30_IK], IKG)
        WeekDate_ref = int([+1979, +52, +7], IKG)
        call runTests(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Values = int([1979_IK, 12_IK, 31_IK], IKG)
        WeekDate_ref = int([+1980, +1, +1], IKG)
        call runTests(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 10
            Values = getDateTime()
            WeekDate = getWeekDate()
            WeekDate_ref = getWeekDate(Values(1:3))
            if (all(WeekDate == WeekDate_ref)) exit
        end do
        call report(__LINE__, SK_"@test_getWeekDate(): Calling `getWeekDate()` without arguments must yield the current week date.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTests(line)
            integer, intent(in) :: line
            WeekDate = getWeekDate(Values)
            call report(line, SK_"@test_getWeekDate(): Calling `getWeekDate(Values(1:3))` must yield the current week date.")
            WeekDate = getWeekDate(Values(1), Values(2), Values(3))
            call report(line, SK_"@test_getWeekDate(): Calling `getWeekDate(Values(1), Values(2), Values(3))` must yield the current week date.")
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, msg)
            use pm_option, only: getOption
            integer, intent(in) :: line
            character(*,SKG), intent(in), optional :: msg
            character(:,SKG), allocatable :: msg_def
            assertion = assertion .and. logical(size(WeekDate, 1, IK) == size(WeekDate_ref, 1, IK), LK)
            call test%assert(assertion, SK_"@test_getWeekDate(): The output week date vector must have the correct length.", int(line, IK))
            assertion = assertion .and. logical(all(WeekDate == WeekDate_ref), LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "WeekDate               ", WeekDate
                write(test%disp%unit,"(*(g0,:,', '))") "WeekDate_ref           ", WeekDate_ref
                write(test%disp%unit,"(*(g0,:,', '))") "size(WeekDate, IK)     ", size(WeekDate, 1, IK)
                write(test%disp%unit,"(*(g0,:,', '))") "size(WeekDate_ref, IK) ", size(WeekDate_ref, 1, IK)
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            if (present(msg)) then
                msg_def = msg
            else
                msg_def = SKG_"@test_getWeekDate(): The output week date vector must be correct."
            end if
            call test%assert(assertion, msg_def, int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end  procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_getWeekYear

        use pm_distUnif, only: getUnifRand
        use pm_kind, only: SKG => SK, IKG => IK, RKG => RK
        integer(IKG) :: weekYear, weekYear_ref
        integer(IKG) :: WeekDate(3), Values(8)
        integer :: i
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        weekYear = getWeekYear()
        WeekDate = getWeekDate()
        weekYear_ref = WeekDate(1)
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 200

            Values = getDateTime(julianDay = getUnifRand(-300000._RKG, +300000._RKG))
            WeekDate = getWeekDate(Values(1:3))
            weekYear_ref = WeekDate(1)

            weekYear = getWeekYear(Values(1:3))
            call report(__LINE__)

            weekYear = getWeekYear(Values(1), Values(2), Values(3))
            call report(__LINE__)

        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, msg)
            use pm_option, only: getOption
            integer, intent(in) :: line
            character(*,SKG), intent(in), optional :: msg
            assertion = assertion .and. logical(weekYear == weekYear_ref, LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "weekYear               ", weekYear
                write(test%disp%unit,"(*(g0,:,', '))") "weekYear_ref           ", weekYear_ref
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, getOption(SK_"@test_getWeekYear(): The output week year must be correct.", msg), int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end  procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_isValidDateTime

        use pm_distUnif, only: getUnifRand
        use pm_kind, only: SKG => SK, IKG => IK, RKG => RK, LKG => LK
        logical(LKG) :: isValid, isValid_ref
        integer(IKG) :: Values(10)
        integer :: i, j
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        isValid_ref = .false._LKG
        isValid = isValidDateTime(Values(2:1))
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        isValid_ref = .false._LKG
        isValid = isValidDateTime(Values(1:10))
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        isValid_ref = .false._LKG
        isValid = isValidDateTime(Values(1:9))
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 200

            Values(1:8) = getDateTime(julianDay = getUnifRand(-300000._RKG, +300000._RKG))

            isValid = isValidDateTime(Values(1), Values(2), Values(3), Values(4), Values(5), Values(6), Values(7), Values(8))
            isValid_ref = isValidDateTime_ref(Values(1:8))
            call report(__LINE__)

            isValid = isValidDateTime(Values(1), Values(2), Values(3), Values(4), Values(5), Values(6), Values(7))
            isValid_ref = isValidDateTime_ref(Values(1:7))
            call report(__LINE__)

            isValid = isValidDateTime(Values(1), Values(2), Values(3), Values(4), Values(5), Values(6))
            isValid_ref = isValidDateTime_ref(Values(1:6))
            call report(__LINE__)

            isValid = isValidDateTime(Values(1), Values(2), Values(3), Values(4), Values(5))
            isValid_ref = isValidDateTime_ref(Values(1:5))
            call report(__LINE__)

            isValid = isValidDateTime(Values(1), Values(2), Values(3), Values(4))
            isValid_ref = isValidDateTime_ref(Values(1:4))
            call report(__LINE__)

            isValid = isValidDateTime(Values(1), Values(2), Values(3))
            isValid_ref = isValidDateTime_ref(Values(1:3))
            call report(__LINE__)

            isValid = isValidDateTime(Values(1), Values(2))
            isValid_ref = isValidDateTime_ref(Values(1:2))
            call report(__LINE__)

            isValid = isValidDateTime(Values(1))
            isValid_ref = isValidDateTime_ref(Values(1:1))
            call report(__LINE__)

            do j = 0, 10
                isValid = isValidDateTime(Values(1:j))
                isValid_ref = isValidDateTime_ref(Values(1:j))
                call report(__LINE__)
            end do

        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function isValidDateTime_ref(Values) result(isValid)
            integer(IKG), intent(in), contiguous :: Values(:)
            logical(LKG) :: isValid
            isValid = 0 < size(Values) .and. size(Values) < 9
            if (isValid .and. size(Values) > 1) then
                isValid = 0_IKG < Values(2) .and. Values(2) < 13_IKG
                if (isValid .and. size(Values) > 2) then
                    if (isLeapYear(Values(1))) then
                        isValid = 0_IKG < Values(3) .and. Values(3) <= DAYS_OF_MONTH_LEAP(Values(2))
                    else
                        isValid = 0_IKG < Values(3) .and. Values(3) <= DAYS_OF_MONTH(Values(2))
                    end if
                    if (isValid .and. size(Values) > 3) then
                        isValid = int(ZONE_MIN, IKG) <= Values(4) .and. Values(4) < int(ZONE_MAX, IKG)
                        if (isValid .and. size(Values) > 4) then
                            isValid = 0_IKG <= Values(5) .and. Values(5) < 24_IKG
                            if (isValid .and. size(Values) > 5) then
                                isValid = 0_IKG <= Values(6) .and. Values(6) < 60_IKG
                                if (isValid .and. size(Values) > 6) then
                                    isValid = 0_IKG <= Values(7) .and. Values(7) < 60_IKG
                                    if (isValid .and. size(Values) > 7) isValid = 0_IKG <= Values(8) .and. Values(8) < 1000_IKG
                                end if
                            end if
                        end if
                    end if
                end if
            end if
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, msg)
            use pm_option, only: getOption
            integer, intent(in) :: line
            character(*,SKG), intent(in), optional :: msg
            assertion = assertion .and. logical(isValid .eqv. isValid_ref, LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "Values         ", Values
                write(test%disp%unit,"(*(g0,:,', '))") "size(Values)   ", size(Values)
                write(test%disp%unit,"(*(g0,:,', '))") "isValid        ", isValid
                write(test%disp%unit,"(*(g0,:,', '))") "isValid_ref    ", isValid_ref
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, getOption(SK_"@test_isValidDateTime(): The procedure must correctly recognize valid and invalid date and time.", msg), int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end  procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_isLastDayInMonth

        use pm_distUnif, only: getUnifRand
        use pm_kind, only: SKG => SK, IKG => IK, RKG => RK, LKG => LK
        integer(IKG) :: Values(8), DateAfter(3)
        logical(LKG) :: isLast, isLast_ref
        integer :: i
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 10
            Values = getDateTime()
            isLast = isLastDayInMonth()
            isLast_ref = isLastDayInMonth(Values)
            if (isLast .eqv. isLast_ref) exit
        end do
        call report(__LINE__)


        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 10000

            Values(1:8) = getDateTime(julianDay = getUnifRand(-300000._RKG, +300000._RKG))
            DateAfter(1:3) = getDateAfter(Values(1:3))
            isLast_ref = logical(DateAfter(2) /= Values(2), LK)

            isLast = isLastDayInMonth(Values(1), Values(2), Values(3))
            call report(__LINE__)

            isLast = isLastDayInMonth(Values(1), Values(2), Values(3))
            call report(__LINE__)

        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, msg)
            use pm_option, only: getOption
            integer, intent(in) :: line
            character(*,SKG), intent(in), optional :: msg
            assertion = assertion .and. logical(isLast .eqv. isLast_ref, LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "Values         ", Values
                write(test%disp%unit,"(*(g0,:,', '))") "isLast_ref     ", isLast_ref
                write(test%disp%unit,"(*(g0,:,', '))") "isLast         ", isLast
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, getOption(SK_"@test_isLastDayInMonth(): The procedure must correctly recognize if input date and time is the last day of the month.", msg), int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end  procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_getDateAfter

        use pm_distUnif, only: getUnifRand
        use pm_kind, only: SKG => SK, IKG => IK, RKG => RK, LKG => LK
        integer(IKG) :: Values(8), DateBefore(3), DateAfter(3), DateAfter_ref(3)
        integer :: i
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 10000

            Values(1:8) = getDateTime(julianDay = getUnifRand(-300000._RKG, +300000._RKG))
            DateAfter_ref(1:3) = Values(1:3)
            DateBefore = getDateBefore(DateAfter_ref)

            DateAfter(1:3) = getDateAfter(DateBefore(1:3))
            call report(__LINE__)

            DateAfter(1:3) = getDateAfter(DateBefore(1), DateBefore(2), DateBefore(3))
            call report(__LINE__)

        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 10
            DateAfter = getDateAfter()
            DateAfter_ref = getDateAfter(getDateTime())
            if (all(DateAfter == DateAfter_ref)) exit
        end do
        call report(__LINE__)


        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, msg)
            use pm_option, only: getOption
            integer, intent(in) :: line
            character(*,SKG), intent(in), optional :: msg
            assertion = assertion .and. logical(all(DateAfter == DateAfter_ref), LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "DateAfter      ", DateAfter
                write(test%disp%unit,"(*(g0,:,', '))") "DateAfter_ref  ", DateAfter_ref
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, getOption(SK_"@test_getDateAfter(): The procedure must correctly compute the date after the specified or the current date.", msg), int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end  procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_getDateBefore

        use pm_distUnif, only: getUnifRand
        use pm_kind, only: SKG => SK, IKG => IK, RKG => RK, LKG => LK
        integer(IKG) :: Values(8), DateAfter(3), DateBefore(3), DateBefore_ref(3)
        integer :: i
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 10000

            Values(1:8) = getDateTime(julianDay = getUnifRand(-300000._RKG, +300000._RKG))
            DateBefore_ref(1:3) = Values(1:3)
            DateAfter = getDateAfter(DateBefore_ref)

            DateBefore(1:3) = getDateBefore(DateAfter(1:3))
            call report(__LINE__)

            DateBefore(1:3) = getDateBefore(DateAfter(1), DateAfter(2), DateAfter(3))
            call report(__LINE__)

        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 10
            DateBefore = getDateBefore()
            DateBefore_ref = getDateBefore(getDateTime())
            if (all(DateBefore == DateBefore_ref)) exit
        end do
        call report(__LINE__)


        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, msg)
            use pm_option, only: getOption
            integer, intent(in) :: line
            character(*,SKG), intent(in), optional :: msg
            assertion = assertion .and. logical(all(DateBefore == DateBefore_ref), LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "DateBefore     ", DateBefore
                write(test%disp%unit,"(*(g0,:,', '))") "DateBefore_ref ", DateBefore_ref
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, getOption(SK_"@test_getDateBefore(): The procedure must correctly compute the date after the specified or the current date.", msg), int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end  procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_getOrdinalDay

        use pm_distUnif, only: getUnifRand
        use pm_kind, only: SKG => SK, IKG => IK, RKG => RK, LKG => LK
        integer(IKG) :: Values(8), ordinalDay, ordinalDay_ref
        integer :: i
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 10000

            Values(1:8) = getDateTime(julianDay = getUnifRand(-300000._RKG, +300000._RKG))
            ordinalDay_ref = floor(getDateTimeDiff(Values, getDateTime(Values(1))), IKG) + 1_IKG

            ordinalDay = getordinalDay(Values)
            call report(__LINE__)

            ordinalDay = getordinalDay(Values(1), Values(2), Values(3))
            call report(__LINE__)

        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 10
            ordinalDay = getOrdinalDay()
            ordinalDay_ref = getOrdinalDay(getDateTime())
            if (ordinalDay == ordinalDay_ref) exit
        end do
        call report(__LINE__)


        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, msg)
            use pm_option, only: getOption
            integer, intent(in) :: line
            character(*,SKG), intent(in), optional :: msg
            assertion = assertion .and. logical(ordinalDay == ordinalDay_ref, LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "ordinalDay     ", ordinalDay
                write(test%disp%unit,"(*(g0,:,', '))") "ordinalDay_ref ", ordinalDay_ref
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, getOption(SK_"@test_getOrdinalDay(): The procedure must correctly compute the ordinal day for the specified or the current date.", msg), int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end  procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_getWeekNumber

        use pm_distUnif, only: getUnifRand
        use pm_kind, only: SKG => SK, IKG => IK, RKG => RK, LKG => LK
        integer(IKG) :: weekNumber, weekNumber_ref
        integer(IKG), allocatable :: Values(:)
        integer :: i
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        weekNumber_ref = +43_IKG
        Values = [integer(IKG) :: +2022, +10, +26]
        call runTests(__LINE__)

        weekNumber_ref = +18_IKG
        Values = [integer(IKG) :: +2022, +5, +8]
        call runTests(__LINE__)

        weekNumber_ref = +19_IKG
        Values = [integer(IKG) :: +2022, +5, +9]
        call runTests(__LINE__)

        weekNumber_ref = +19_IKG
        Values = [integer(IKG) :: +2022, +5, +12]
        call runTests(__LINE__)

        weekNumber_ref = +48_IKG
        Values = [integer(IKG) :: -4713_IK, 11_IK, 24_IK]
        call runTests(__LINE__)

        weekNumber_ref = +53_IKG
        Values = [integer(IKG) :: -1_IK, 1_IK, 1_IK]
        call runTests(__LINE__)

        weekNumber_ref = +52_IKG
        Values = [integer(IKG) :: -1_IK, 12_IK, 31_IK]
        call runTests(__LINE__)

        weekNumber_ref = +52_IKG
        Values = [integer(IKG) :: 0_IK, 12_IK, 31_IK]
        call runTests(__LINE__)

        weekNumber_ref = +1_IKG
        Values = [integer(IKG) :: 1_IK, 12_IK, 31_IK]
        call runTests(__LINE__)

        weekNumber_ref = +41_IKG
        Values = [integer(IKG) :: 1582_IK, 10_IK, 15_IK]
        call runTests(__LINE__)

        weekNumber_ref = +1_IKG
        Values = [integer(IKG) :: 1901_IK, 1_IK, 1_IK]
        call runTests(__LINE__)

        weekNumber_ref = +9_IKG
        Values = [integer(IKG) :: 1999_IK, 3_IK, 1_IK]
        call runTests(__LINE__)

        weekNumber_ref = +9_IKG
        Values = [integer(IKG) :: 2000_IK, 3_IK, 1_IK]
        call runTests(__LINE__)

        weekNumber_ref = +15_IKG
        Values = [integer(IKG) :: 1999_IK, 4_IK, 15_IK]
        call runTests(__LINE__)

        weekNumber_ref = +15_IKG
        Values = [integer(IKG) :: 2000_IK, 4_IK, 15_IK]
        call runTests(__LINE__)

        weekNumber_ref = +52_IKG
        Values = [integer(IKG) :: 9999_IK, 12_IK, 31_IK]
        call runTests(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 10
            weekNumber = getWeekNumber()
            weekNumber_ref = getWeekNumber(getDateTime())
            if (weekNumber == weekNumber_ref) exit
        end do
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTests(line, msg)
            integer, intent(in) :: line
            character(*,SKG), intent(in), optional :: msg
            weekNumber = getWeekNumber(Values)
            call report(line, msg)
            weekNumber = getWeekNumber(Values(1), Values(2), Values(3))
            call report(line, msg)
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, msg)
            use pm_option, only: getOption
            integer, intent(in) :: line
            character(*,SKG), intent(in), optional :: msg
            assertion = assertion .and. logical(weekNumber == weekNumber_ref, LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "Values         ", Values
                write(test%disp%unit,"(*(g0,:,', '))") "weekNumber     ", weekNumber
                write(test%disp%unit,"(*(g0,:,', '))") "weekNumber_ref ", weekNumber_ref
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, getOption(SK_"@test_getWeekNumber(): The procedure must correctly compute the week number for the specified or the current date.", msg), int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end  procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_getWeekDay

        use pm_distUnif, only: getUnifRand
        use pm_kind, only: SKG => SK, IKG => IK, RKG => RK, LKG => LK
        integer(IKG) :: weekDay, weekDay_ref
        integer(IKG), allocatable :: Values(:)
        integer :: i
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        weekDay_ref = +5_IKG
        Values = [integer(IKG) :: 1582_IK, 10_IK, 15_IK]
        call runTests(__LINE__)

        weekDay_ref = +1_IKG
        Values = [integer(IKG) :: 1_IK, 1_IK, 1_IK]
        call runTests(__LINE__)

        weekDay_ref = +6_IKG
        Values = [integer(IKG) :: 2000_IK, 1_IK, 1_IK]
        call runTests(__LINE__)

        weekDay_ref = +0_IKG
        Values = [integer(IKG) :: 2022_IK, 11_IK, 6_IK]
        call runTests(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 10
            weekDay = getWeekDay()
            weekDay_ref = getWeekDay(getDateTime())
            if (weekDay == weekDay_ref) exit
        end do
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTests(line, msg)
            integer, intent(in) :: line
            character(*,SKG), intent(in), optional :: msg
            weekDay = getWeekDay(Values)
            call report(line, msg)
            weekDay = getWeekDay(Values(1), Values(2), Values(3))
            call report(line, msg)
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, msg)
            use pm_option, only: getOption
            integer, intent(in) :: line
            character(*,SKG), intent(in), optional :: msg
            assertion = assertion .and. logical(weekDay == weekDay_ref, LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "Values         ", Values
                write(test%disp%unit,"(*(g0,:,', '))") "weekDay     ", weekDay
                write(test%disp%unit,"(*(g0,:,', '))") "weekDay_ref ", weekDay_ref
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, getOption(SK_"@test_getWeekDay(): The procedure must correctly compute the week day for the specified or the current date.", msg), int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end  procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_getWeekDayISO

        use pm_distUnif, only: getUnifRand
        use pm_kind, only: SKG => SK, IKG => IK, RKG => RK, LKG => LK
        integer(IKG) :: weekDay, weekDay_ref
        integer(IKG), allocatable :: Values(:)
        integer :: i
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        weekDay_ref = +5_IKG
        Values = [integer(IKG) :: 1582_IK, 10_IK, 15_IK]
        call runTests(__LINE__)

        weekDay_ref = +1_IKG
        Values = [integer(IKG) :: 1_IK, 1_IK, 1_IK]
        call runTests(__LINE__)

        weekDay_ref = +6_IKG
        Values = [integer(IKG) :: 2000_IK, 1_IK, 1_IK]
        call runTests(__LINE__)

        weekDay_ref = +7_IKG
        Values = [integer(IKG) :: 2022_IK, 11_IK, 6_IK]
        call runTests(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 10
            weekDay = getWeekDayISO()
            weekDay_ref = getWeekDayISO(getDateTime())
            if (weekDay == weekDay_ref) exit
        end do
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTests(line, msg)
            integer, intent(in) :: line
            character(*,SKG), intent(in), optional :: msg
            weekDay = getWeekDayISO(Values)
            call report(line, msg)
            weekDay = getWeekDayISO(Values(1), Values(2), Values(3))
            call report(line, msg)
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, msg)
            use pm_option, only: getOption
            integer, intent(in) :: line
            character(*,SKG), intent(in), optional :: msg
            assertion = assertion .and. logical(weekDay == weekDay_ref, LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "Values      ", Values
                write(test%disp%unit,"(*(g0,:,', '))") "weekDay     ", weekDay
                write(test%disp%unit,"(*(g0,:,', '))") "weekDay_ref ", weekDay_ref
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, getOption(SK_"@test_getWeekDayISO(): The procedure must correctly compute the ISO week day for the specified or the current date.", msg), int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end  procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_getCountDays

        use pm_distUnif, only: getUnifRand
        use pm_kind, only: SKG => SK, IKG => IK, RKG => RK, LKG => LK
        integer(IKG) :: countDays, countDays_ref
        integer(IKG), allocatable :: Values(:)
        integer :: i
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 1000

            Values = getDateTime(julianDay = getUnifRand(-300000._RKG, +300000._RKG))
            countDays_ref = merge(366_IKG, 365_IKG, isLeapYear(Values(1)))
            countDays = getCountDays(Values(1))
            call report(__LINE__)

            Values(1:8) = getDateTime(julianDay = getUnifRand(-300000._RKG, +300000._RKG))
            countDays_ref = merge(DAYS_OF_MONTH_LEAP(Values(2)), DAYS_OF_MONTH(Values(2)), isLeapYear(Values(1)))
            countDays = getCountDays(Values(1), Values(2))
            call report(__LINE__)

        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, msg)
            use pm_option, only: getOption
            integer, intent(in) :: line
            character(*,SKG), intent(in), optional :: msg
            assertion = assertion .and. logical(countDays == countDays_ref, LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "Values      ", Values
                write(test%disp%unit,"(*(g0,:,', '))") "countDays     ", countDays
                write(test%disp%unit,"(*(g0,:,', '))") "countDays_ref ", countDays_ref
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, getOption(SK_"@test_getCountDays(): The procedure must correctly compute the day count for the specified or the current date.", msg), int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end  procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_getCountWeeks

        use pm_arrayRange, only: getRange
        use pm_distUnif, only: getUnifRand
        use pm_kind, only: SKG => SK, IKG => IK, RKG => RK, LKG => LK
        integer(IKG), allocatable :: CountWeeks(:), CountWeeks_ref(:)
        integer :: i
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CountWeeks_ref = [integer(IKG) :: +52, +52, +52, +52, +53, +52, +52, +52, +52, +53, +52, +52, +52, +52, +52, +53, +52, +52, +52, +52, +53]
        CountWeeks = getCountWeeks(getRange(2000_IKG, 2020_IKG))
        call report(__LINE__)

        CountWeeks_ref = [integer(IKG) :: +4, +4, +5, +5, +5, +4, +4, +4, +5, +5, +4, +4, +4, +5, +5, +5, +4, +4, +4, +5, +5]
        CountWeeks = getCountWeeks(getRange(2000_IKG, 2020_IKG), 1_IKG)
        call report(__LINE__)

        CountWeeks_ref = [integer(IKG) :: +4, +4, +4, +4, +4, +4, +4, +4, +4, +4, +4, +4, +4, +4, +4, +4, +4, +4, +4, +4, +4]
        CountWeeks = getCountWeeks(getRange(2000_IKG, 2020_IKG), 2_IKG)
        call report(__LINE__)

        CountWeeks_ref = [integer(IKG) :: +5, +5, +4, +4, +4, +5, +5, +5, +4, +4, +4, +5, +5, +4, +4, +4, +5, +5, +5, +4, +4]
        CountWeeks = getCountWeeks(getRange(2000_IKG, 2020_IKG), 3_IKG)
        call report(__LINE__)

        CountWeeks_ref = [integer(IKG) :: +4, +4, +4, +4, +5, +4, +4, +4, +4, +5, +5, +4, +4, +4, +4, +5, +4, +4, +4, +4, +5]
        CountWeeks = getCountWeeks(getRange(2000_IKG, 2020_IKG), 4_IKG)
        call report(__LINE__)

        CountWeeks_ref = [integer(IKG) :: +4, +5, +5, +5, +4, +4, +4, +5, +5, +4, +4, +4, +5, +5, +5, +4, +4, +4, +5, +5, +4]
        CountWeeks = getCountWeeks(getRange(2000_IKG, 2020_IKG), 5_IKG)
        call report(__LINE__)

        CountWeeks_ref = [integer(IKG) :: +5, +4, +4, +4, +4, +5, +5, +4, +4, +4, +4, +5, +4, +4, +4, +4, +5, +5, +4, +4, +4]
        CountWeeks = getCountWeeks(getRange(2000_IKG, 2020_IKG), 6_IKG)
        call report(__LINE__)

        CountWeeks_ref = [integer(IKG) :: +4, +4, +4, +5, +5, +4, +4, +4, +5, +5, +5, +4, +4, +4, +5, +5, +4, +4, +4, +4, +5]
        CountWeeks = getCountWeeks(getRange(2000_IKG, 2020_IKG), 7_IKG)
        call report(__LINE__)

        CountWeeks_ref = [integer(IKG) :: +5, +5, +5, +4, +4, +4, +5, +5, +4, +4, +4, +4, +5, +5, +4, +4, +4, +5, +5, +5, +4]
        CountWeeks = getCountWeeks(getRange(2000_IKG, 2020_IKG), 8_IKG)
        call report(__LINE__)

        CountWeeks_ref = [integer(IKG) :: +4, +4, +4, +4, +5, +5, +4, +4, +4, +4, +5, +5, +4, +4, +4, +4, +5, +4, +4, +4, +4]
        CountWeeks = getCountWeeks(getRange(2000_IKG, 2020_IKG), 9_IKG)
        call report(__LINE__)

        CountWeeks_ref = [integer(IKG) :: +4, +4, +5, +5, +4, +4, +4, +4, +5, +5, +4, +4, +4, +5, +5, +5, +4, +4, +4, +5, +5]
        CountWeeks = getCountWeeks(getRange(2000_IKG, 2020_IKG), 10_IKG)
        call report(__LINE__)

        CountWeeks_ref = [integer(IKG) :: +5, +5, +4, +4, +4, +4, +5, +5, +4, +4, +4, +4, +5, +4, +4, +4, +4, +5, +5, +4, +4]
        CountWeeks = getCountWeeks(getRange(2000_IKG, 2020_IKG), 11_IKG)
        call report(__LINE__)

        CountWeeks_ref = [integer(IKG) :: +4, +4, +4, +4, +5, +5, +4, +4, +4, +5, +5, +5, +4, +4, +4, +5, +5, +4, +4, +4, +5]
        CountWeeks = getCountWeeks(getRange(2000_IKG, 2020_IKG), 12_IKG)
        call report(__LINE__)

        CountWeeks_ref = [integer(IKG) :: +4, +4, +5, +4, +4, +5, +4, +4, +5, +4, +4, +5]
        CountWeeks = getCountWeeks(2022_IKG, getRange(1_IKG, 12_IKG))
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, msg)
            use pm_option, only: getOption
            integer, intent(in) :: line
            character(*,SKG), intent(in), optional :: msg
            assertion = assertion .and. logical(size(CountWeeks) == size(CountWeeks_ref), LK)
            call test%assert(assertion, getOption(SK_"@test_getCountWeeks(): The sizes of `CountWeeks` and `CountWeeks_ref` must match.", msg), int(line, IK))
            do i = 1, size(CountWeeks)
                assertion = assertion .and. logical(CountWeeks(i) == CountWeeks_ref(i), LK)
                if (test%traceable .and. .not. assertion) then
                    ! LCOV_EXCL_START
                    write(test%disp%unit,"(*(g0,:,', '))")
                    write(test%disp%unit,"(*(g0,:,', '))") "CountWeeks     ", CountWeeks
                    write(test%disp%unit,"(*(g0,:,', '))") "CountWeeks_ref ", CountWeeks_ref
                    write(test%disp%unit,"(*(g0,:,', '))")
                    ! LCOV_EXCL_STOP
                end if
                call test%assert(assertion, getOption(SK_"@test_getCountWeeks(): The procedure must correctly compute the week counts for the specified or the current date.", msg), int(line, IK))
            end do
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end  procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_getCountLeapYears

        use pm_arrayRange, only: getRange
        use pm_distUnif, only: getUnifRand
        use pm_kind, only: SKG => SK, IKG => IK, RKG => RK, LKG => LK
        integer(IKG), allocatable :: CountLeapYears(:), CountLeapYears_ref(:)
        integer :: i
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CountLeapYears_ref = [integer(IKG) :: 0]
        CountLeapYears = getCountLeapYears(until = [1_IKG])
        call report(__LINE__)

        CountLeapYears_ref = [integer(IKG) :: +485]
        CountLeapYears = getCountLeapYears(until = [2000_IKG])
        call report(__LINE__)

        CountLeapYears_ref = [integer(IKG) :: +490]
        CountLeapYears = [getCountLeapYears(until = 2022_IKG)]
        call report(__LINE__)

        CountLeapYears_ref = [integer(IKG) :: +490]
        CountLeapYears = [getCountLeapYears(until = 2022_IKG, since = 1_IKG)]
        call report(__LINE__)

        CountLeapYears_ref = [integer(IKG) :: +0, +0, +0, +1]
        CountLeapYears = getCountLeapYears(until = [1_IKG, 2_IKG, 3_IKG, 4_IKG])
        call report(__LINE__)

        CountLeapYears_ref = [integer(IKG) :: +0, +1, +1, +1, +1, +2, +2, +2, +2, +3, +3]
        CountLeapYears = getCountLeapYears(until = getRange(-5_IKG, 5_IKG, 1_IKG), since = -5_IKG)
        call report(__LINE__)

        CountLeapYears_ref = [integer(IKG) :: -3, -2, -2, -2, -2, -1, -1, -1, -1, +0, +0, +0, +1, +1, +1, +1, +2]
        CountLeapYears = getCountLeapYears(until = getRange(-8_IKG, 8_IKG, 1_IKG))
        call report(__LINE__)

        CountLeapYears_ref = [integer(IKG) :: +6, +6, +5]
        CountLeapYears = getCountLeapYears(until = 2022_IKG, since = [1999_IKG, 2000_IKG, 2001_IKG])
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, msg)
            use pm_option, only: getOption
            integer, intent(in) :: line
            character(*,SKG), intent(in), optional :: msg
            assertion = assertion .and. logical(size(CountLeapYears) == size(CountLeapYears_ref), LK)
            call test%assert(assertion, SK_"@test_getCountLeapYears(): The sizes of `CountLeapYears` and `CountLeapYears_ref` must match.", int(line, IK))
            do i = 1, size(CountLeapYears)
                assertion = assertion .and. logical(CountLeapYears(i) == CountLeapYears_ref(i), LK)
                if (test%traceable .and. .not. assertion) then
                    ! LCOV_EXCL_START
                    write(test%disp%unit,"(*(g0,:,', '))")
                    write(test%disp%unit,"(*(g0,:,', '))") "CountLeapYears     ", CountLeapYears
                    write(test%disp%unit,"(*(g0,:,', '))") "CountLeapYears_ref ", CountLeapYears_ref
                    write(test%disp%unit,"(*(g0,:,', '))")
                    ! LCOV_EXCL_STOP
                end if
                call test%assert(assertion, getOption(SK_"@test_getCountLeapYears(): The procedure must correctly compute the leap year counts for the specified or the current date.", msg), int(line, IK))
            end do
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end  procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_isValidZone

        use pm_arrayRange, only: getRange
        use pm_distUnif, only: getUnifRand
        use pm_kind, only: SKG => SK, IKG => IK, RKG => RK, LKG => LK
        logical(LKG), allocatable :: IsValid(:), IsValid_ref(:)
        integer(IKG), allocatable :: Zone(:)
        integer :: i
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 200
            Zone = getUnifRand(-1000_IKG, 1000_IKG, getUnifRand(1_IK, 20_IK))
            IsValid_ref = int(ZONE_MIN, IKG) <= Zone .and. Zone <= int(ZONE_MAX, IKG)
            if (size(Zone) == 1 .and. getUnifRand()) then
                IsValid = [isValidZone(Zone(1))]
            else
                IsValid = isValidZone(Zone)
            end if
            call report(__LINE__)
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, msg)
            use pm_val2str, only: getStr
            use pm_option, only: getOption
            integer, intent(in) :: line
            character(*,SKG), intent(in), optional :: msg
            integer :: i
            assertion = assertion .and. logical(size(IsValid) == size(IsValid_ref), LK)
            call test%assert(assertion, SK_"@test_getCountLeapYears(): The sizes of `IsValid` and `IsValid_ref` must match. size(IsValid), size(IsValid) = "//getStr([size(IsValid), size(IsValid)]), int(line, IK))
            do i = 1, size(IsValid)
                assertion = assertion .and. logical(IsValid(i) .eqv. IsValid_ref(i), LK)
                if (test%traceable .and. .not. assertion) then
                    ! LCOV_EXCL_START
                    write(test%disp%unit,"(*(g0,:,', '))")
                    write(test%disp%unit,"(*(g0,:,', '))") "Zone           ", Zone
                    write(test%disp%unit,"(*(g0,:,', '))") "IsValid        ", IsValid
                    write(test%disp%unit,"(*(g0,:,', '))") "IsValid_ref    ", IsValid_ref
                    write(test%disp%unit,"(*(g0,:,', '))")
                    ! LCOV_EXCL_STOP
                end if
                call test%assert(assertion, getOption(SK_"@test_isValidZone(): The procedure must correctly recognize a valid zone.", msg), int(line, IK))
            end do
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end  procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_isMorning

        use pm_arrayRange, only: getRange
        use pm_distUnif, only: getUnifRand
        use pm_kind, only: SKG => SK, IKG => IK, RKG => RK, LKG => LK
        logical(LKG), allocatable :: Morning(:), Morning_ref(:)
        integer(IKG), allocatable :: Zone(:)
        real(RKG), allocatable :: JulianDay(:), JulianDayUTC(:)
        integer :: i, j
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 500

            do j = 1, 10
                Zone = getUnifRand(int(ZONE_MIN, IKG), int(ZONE_MAX, IKG), getUnifRand(1_IK, 10_IK))
                if (size(Zone) == 1 .and. getUnifRand()) then
                    Morning = [isMorning(Zone(1))]
                else
                    Morning = [isMorning(Zone)]
                end if
                Morning_ref = [getHour(Zone) < 12_IKG]
                if (all(Morning .eqv. Morning_ref)) exit
            end do
            call report(__LINE__, Zone = Zone)

            JulianDay = getUnifRand(-300000._RKG, +300000._RKG, getUnifRand(1_IK, 20_IK))
            Morning_ref = logical(JulianDay - real(floor(JulianDay, IK), RKG) >= 0.5_RKG, LK)
            if (size(JulianDay) == 1 .and. getUnifRand()) then
                Morning = [isMorning(JulianDay(1))]
            else
                Morning = isMorning(JulianDay)
            end if
            call report(__LINE__, JulianDay = JulianDay)

            Zone = getUnifRand(int(ZONE_MIN, IKG), int(ZONE_MAX, IKG), size(JulianDay, kind = IK))
            if (size(JulianDay) == 1 .and. getUnifRand()) then
                Morning = [isMorning(JulianDay(1), Zone(1))]
            else
                Morning = isMorning(JulianDay, Zone)
            end if
            JulianDayUTC = JulianDay + Zone / real(MINUTES_PER_DAY, RKG)
            Morning_ref = logical(JulianDayUTC - real(floor(JulianDayUTC, IK), RKG) >= 0.5_RKG, LK)
            call report(__LINE__, JulianDay = JulianDay, Zone = Zone)

        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, 10
            Morning = [isMorning()]
            Morning_ref = [getHour() < 12_IKG]
            if (all(Morning .eqv. Morning_ref)) exit
        end do
        call report(__LINE__)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, msg, JulianDay, Zone)
            use pm_val2str, only: getStr
            use pm_option, only: getOption
            integer, intent(in) :: line
            character(*,SKG), intent(in), optional :: msg
            integer(IKG), intent(in), contiguous, optional :: Zone(:)
            real(RKG), intent(in), contiguous, optional :: JulianDay(:)
            integer :: i
            assertion = assertion .and. logical(size(Morning) == size(Morning_ref), LK)
            call test%assert(assertion, SK_"@test_getCountLeapYears(): The sizes of `Morning` and `Morning_ref` must match. size(Morning), size(Morning) = "//getStr([size(Morning), size(Morning)]), int(line, IK))
            do i = 1, size(Morning)
                assertion = assertion .and. logical(Morning(i) .eqv. Morning_ref(i), LK)
                if (test%traceable .and. .not. assertion) then
                    ! LCOV_EXCL_START
                    write(test%disp%unit,"(*(g0,:,', '))")
                    if (present(JulianDay)) then
                    write(test%disp%unit,"(*(g0,:,', '))") "JulianDay      ", JulianDay
                    end if
                    if (present(Zone)) then
                    write(test%disp%unit,"(*(g0,:,', '))") "Zone           ", Zone
                    end if
                    write(test%disp%unit,"(*(g0,:,', '))") "Morning        ", Morning
                    write(test%disp%unit,"(*(g0,:,', '))") "Morning_ref    ", Morning_ref
                    write(test%disp%unit,"(*(g0,:,', '))")
                    ! LCOV_EXCL_STOP
                end if
                call test%assert(assertion, getOption(SK_"@test_isMorning(): The procedure must correctly recognize morning hours.", msg), int(line, IK))
            end do
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end  procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_getZone
        use pm_arrayRange, only: getRange
        use pm_distUnif, only: getUnifRand
        use pm_kind, only: SKG => SK, IKG => IK, RKG => RK, LKG => LK
        integer(IKG) :: zone, zone_ref, Values(8)
        integer :: i
        assertion = .true._LK
        do i = 1, 500
            zone = getZone()
            call date_and_time(Values = Values)
            zone_ref = ValueS(4)
            call report(__LINE__)
        end do
    contains
        subroutine report(line)
            use pm_val2str, only: getStr
            integer, intent(in) :: line
            assertion = assertion .and. logical(zone == zone_ref, LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "zone    ", zone
                write(test%disp%unit,"(*(g0,:,', '))") "zone_ref", zone_ref
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"@test_getZone(): The local zone must be inferred correctly.", int(line, IK))
        end subroutine
    end  procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_getMillisecond
        use pm_arrayRange, only: getRange
        use pm_distUnif, only: getUnifRand
        use pm_kind, only: SKG => SK, IKG => IK, RKG => RK, LKG => LK
        integer(IKG) :: millisecond
        integer :: i
        assertion = .true._LK
        do i = 1, 500
            millisecond = getMillisecond()
            call report(__LINE__)
        end do
    contains
        subroutine report(line)
            use pm_val2str, only: getStr
            integer, intent(in) :: line
            assertion = assertion .and. logical(0_IKG <= millisecond .and. millisecond < 1000_IKG, LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "millisecond", millisecond
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"@test_getMillisecond(): The local millisecond must be inferred correctly within the limits.", int(line, IK))
        end subroutine
    end  procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if 0
#define getDisMahalSq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getDisMahalSq_CK5
        use pm_kind, only: IK, CK => CK5
#include "test_pm_dateTime@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getDisMahalSq_CK4
        use pm_kind, only: IK, CK => CK4
#include "test_pm_dateTime@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getDisMahalSq_CK3
        use pm_kind, only: IK, CK => CK3
#include "test_pm_dateTime@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getDisMahalSq_CK2
        use pm_kind, only: IK, CK => CK2
#include "test_pm_dateTime@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getDisMahalSq_CK1
        use pm_kind, only: IK, CK => CK1
#include "test_pm_dateTime@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDisMahalSq_ENABLED
#endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE