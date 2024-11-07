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

!>  \brief
!>  This module contains classes and procedures for computing, manipulating, and styling dates and times.
!>
!>  \details
!>  This module strives to follow the conventions of \f$\ms{ISO 8601:2004}\f$.<br>
!>  -#  <b>\f$\ms{ISO 8601:2004}\f$</b>
!>      -#  The *International Organization for Standardization* (**ISO**) is the regulator of <i>\f$\ms{ISO 8601:2004}\f$</i>,
!>          which is the international standard for covering the exchange of data related to date or time.
!>      -#  Though it is used by the majority of the world, especially developed, not all adhere to this standard.
!>  -#  **Gregorian Calendar**
!>      -#  The Gregorian Calendar, named after Pope Gregory XIII, is based on the time it takes the moon to make
!>          one full revolution around the Earth (roughly one month) and the Sun to make a full revolution around the Earth (roughly one year).
!>      -#  Because celestial bodies of such close proximity can be widely encountered, it was common in ancient times to use them to tell
!>          time and/or date. This method was further refined into what we now know as the current Gregorian Calendar.
!>
!>  \details
!>  **A breif history of Julian Calendar and its evolution toward Gregorian calendar**<br>
!>
!>  -#  The Julian calendar is important to historians because it was used worldwide for over 16 centuries,
!>      and in various parts of the world for another three centuries after that.
!>  -#  It is also important to genealogists because it was used to record events in many countries as recently as the early 1900s.
!>  -#  There are several subtle but significant differences between the Julian Calendar and its successor that is currently used (the Gregorian Calendar.<br>
!>      For example, the birthday of George Washington is conventionally celebrated on February 22. However, as a result of the switch from the Julian calendar to the Gregorian,
!>      his birthday was in fact February 23 in the 19th century, February 24th in the 20th and 21st century, and will continue to advance in future centuries.
!>
!>  **The Calendar Requirements**
!>
!>  -#  The calendar has only one basic requirement -- that **the seasons do not migrate through the years**.<br>
!>      For example, We are used to going to the beach in July, and if after a few years it does not warm up sufficiently until November we might not be too happy.
!>  -#  The seasons are determined by the position of the Earth in its orbit around the sun.
!>  -#  Earth goes around the sun once every \f$\ms{365.2422}\f$ days.
!>  -#  The fractional part of the year \f$\ms{.2422}\f$ in units of days makes the design of precise calendars very difficult.
!>  -#  Since it is not desired to start a new year after a fractional number of days, the compromise is to pick some integral number of days that is close to \f$\ms{365.2422}\f$.
!>  -#  Picking an integral year length that is too low will shift the seasons to later and later times each year.<br>
!>      For example, if the year is taken to be 300 days, the earth will not yet have reached the correct point in its orbit for the winter next year and the winter season will not start until `65` days into the second year.
!>      This delay will propagate further into the following years and the third winter it will not start until `130` days into the third year.
!>      On the other hand, if the year is taken to be too many days, for example, `400`, the seasons will come earlier and earlier each year.
!>  -#  The closer the chosen length of the year is to the magical number \f$\ms{365.2422}\f$, such as `365` or `366`, the slower the drift in seasons will occur, but it will not not eliminate it.<br>
!>      If we wait long enough, even the best designed calendar will start to experience drift.
!>
!>  **A Historical Perspective**
!>
!>  -#  The seeds of the current Gregorian calendar lie in the calendars used in ancient Rome.
!>
!>  <b>The Calendar of Romulus, \f$\ms{753 BC}\f$</b>
!>
!>  -#  The first Roman calendar was introduced in approximately \f$\ms{753 BC}\f$ by Romulus, the first king of Rome.<br>
!>  -#  The calendar of Romulus had only ten months with each month having either `30` or `31` days as follows:
!>      -#  Martius -- 31 days
!>      -#  Aprilis -- 30 days
!>      -#  Maius -- 31 days
!>      -#  Iunius -- 30 days
!>      -#  Quintilis -- 31 days
!>      -#  Sextilis -- 30 days
!>      -#  September -- 30 days
!>      -#  October -- 31 days
!>      -#  November -- 30 days
!>      -#  December -- 30 days
!>
!>  -#  It is not clear where the name *Aprilis* came from, but the other three of the first four months were named after Roman gods.
!>  -#  Starting with the fifth month, the names reflect the order of the month in the calendar: *Quintilis* (meaning five in Latin), *Sextilis* (meaning six), ..., to *December* (meaning ten).
!>  -#  The total number of days in the ten months was `304` which is far less than the year length required to keep the seasons fixed.
!>  -#  People had to wait around until astronomers determined that it was time to start the next year.
!>  -#  This ancient calendar system left about `61` winter days unaccounted for that were not in any month.
!>
!>  <b>Calendar of Numa Pompilius, \f$\ms{713 BC}\f$</b>
!>
!>  -#  By \f$\ms{713 BC}\f$ the calendar was modified by *Numa Pompilius*, the second king of Rome.
!>  -#  Pompilius added two months **at the end** of the calendar, **Ianuarius** and **Februarius** to compensate for the unaccounted-for days.
!>  -#  Pompilius also introduced an *intercalary month* that occurred after *Februarius* in certain years. These years became known as **leap years**.
!>  -#  Pompilius also deleted one day from all the months that had `30` days so that they had `29` days instead.
!>      -#  Martius -- 31 days
!>      -#  Aprilis -- 29 days
!>      -#  Maius -- 31 days
!>      -#  Iunius -- 29 days
!>      -#  Quintilis -- 31 days
!>      -#  Sextilis -- 29 days
!>      -#  September -- 29 days
!>      -#  October -- 31 days
!>      -#  November -- 29 days
!>      -#  December -- 29 days
!>      -#  Ianuarius -- 29 days
!>      -#  Februarius -- 28 days (23 or 24 days in leap year)
!>      -#  Intercalarius -- 0 days (27 days in leap year)
!>  -#  This resulted in a total of `355` days in a common year and `377` days in a leap year.
!>  -#  This required to have a leap year just about every other year.
!>  -#  However, the Numa Calendar did not have a hard and fast rule to keep up with the leap years.<br>
!>      Instead it was left to the whim of the king, who would frequently choose the leap years for political gain rather than for sound astronomical reasons.
!>  -#  The lack of hard rules for the leap years made the Numa Calendar unstable, although it remained in widespread usage for the next `700` years.
!>  -#  By sometime around or before \f$\ms{450 BC}\f$ the starting point of the calendar was shifted from **Martius** to **Ianuarius**.<br>
!>      All other aspects of the calendar remained the same, so it was still effectively the Numa calendar.<br>
!>  -#  **This change of starting point resulted in month names that are now misnomers**, no longer corresponding to their position in the calendar.
!>  -#  The Numa Calendar in \f$\ms{450 BC}\f$ was as follows:
!>      -#  Ianuarius -- 29 days
!>      -#  Februarius -- 28 days (23 or 24 days in leap year)
!>      -#  Intercalarius -- 0 days (27 days in leap year)
!>      -#  Martius -- 31 days
!>      -#  Aprilis -- 29 days
!>      -#  Maius -- 31 days
!>      -#  Iunius -- 29 days
!>      -#  Quintilis -- 31 days
!>      -#  Sextilis -- 29 days
!>      -#  September -- 29 days
!>      -#  October -- 31 days
!>      -#  November -- 29 days
!>      -#  December -- 29 days
!>
!>  <b>The Calendar of Julius Caesar in \f$\ms{45 BC}\f$</b>
!>
!>  -#  Eventually the abuse of the leap years in the Numa Calendar became so egregious that the harvest festival was came before the summer planting season.
!>  -#  Therefore, in \f$\ms{45 BC}\f$ Julius Caesar reformed the calendar and introduced **the first stable calendar**.
!>  -#  Caesar incorporated fixed rules for determining which years were leap years.
!>  -#  Caesar also eliminated the **intercalary month** and replaced it with a single **intercalary day**.
!>  -#  Caesar also introduced a regular pattern of alternating `31` and `30` day counts in months.
!>  -#  Caesar also did a one-time insertion of three months in the year \f$\ms{46 BC}\f$ to give the seasons a chance to catch up.
!>  -#  The calendar of Julius Caesar is as follows:
!>      -#  Ianuarius -- 31 days
!>      -#  Februarius -- 29 days (30 days in leap year)
!>      -#  Martius -- 31 days
!>      -#  Aprilis -- 30 days
!>      -#  Maius -- 31 days
!>      -#  Iunius -- 30 days
!>      -#  Quintilis -- 31 days
!>      -#  Sextilis -- 30 days
!>      -#  September -- 31 days
!>      -#  October -- 30 days
!>      -#  November -- 31 days
!>      -#  December -- 30 days
!>  -#  The calendar of Julius Caesar required that every third year shall be a leap year.
!>  -#  It is believed that Julius intended for it to be every fourth year but the people who implemented it made a calculation error in the way they counted to four (the so-called **fence-post error**).
!>  -#  In spite of the leap-year error, every year of the calendar was now \f$\ms{365.3333}\f$ days on average, a number very close to the correct number of \f$\ms{365.2422}\f$.
!>
!>  <b>The Julian Calendar in \f$\ms{44 BC}\f$</b>
!>
!>  -#  Julius Caesar was killed on the *Ides of Martius* (March 15) in \f$\ms{44 BC}\f$, one year after his calendar went into effect.
!>  -#  The successor to Julius Caesar, Augustus Caesar, made some refinements to the calendar.
!>  -#  The changes introduced by Augustus created a stable calendar known as the **Julian calendar** which was used world-wide for the next `16` centuries.
!>  -#  Augustus changed the leap-year cycle to every four years instead of every three.
!>  -#  Augustus also renamed **Quintilis** to **Iulius**, to honor his predecessor, Julius Caesar.
!>  -#  Augustus also decided to give homage to himself by changing **Sextilis** to **Augustus**.
!>  -#  However, the month of Augustus had one fewer day than the month of Julius Caesar which made Augustus unhappy.
!>  -#  Therefore, Augustus also removed a day from *Februarius* and added it to *Augustus*, making it a `31`-day month, same as *Iulius*.
!>  -#  But since *September* had `31` days as well, that would make three consecutive months with `31` days.
!>  -#  Hence, Augustus also interchanged the number of days in *September* and *October*, as well as interchanging the number of days in *November* and *December*.
!>  -#  The revised calendar of Augustus (the Julian Calendar) is as follows:
!>      -#  Ianuarius -- 31 days
!>      -#  Februarius -- 28 days (29 days in leap year)
!>      -#  Martius -- 31 days
!>      -#  Aprilis -- 30 days
!>      -#  Maius -- 31 days
!>      -#  Iunius -- 30 days
!>      -#  Iulius -- 31 days
!>      -#  Augustus -- 31 days
!>      -#  September -- 30 days
!>      -#  October -- 31 days
!>      -#  November -- 30 days
!>      -#  December -- 31 days
!>  -#  **The above names and number of days of the months of the Julian Calendar remain the same up to the present day**.
!>  -#  Augustus also corrected the *fence-post error* by requiring three-year cycle for leap years.
!>
!>  <b>The Gregorian Calendar in \f$\ms{October 15, 1582 AD}\f$</b>
!>
!>  -#  The leap-year correction introduced by Augustus led to a Julian calendar that had an average of \f$\ms{365.25}\f$ days per year.
!>  -#  However, the corrected year length  was still slightly off the true number of \f$\ms{365.2422}\f$.
!>  -#  A difference of `0.0078` days per year yields `1` day every `128` years or about `3` days every `400` years.
!>  -#  By the 1500s, the Julian Calendar error amounted to about `10` days. The holidays were becoming noticeably misaligned with the seasons.
!>  -#  To get back in step, Pope Gregory XIII decreed in October 1582 that `10` days be stricken from the calendar.
!>  -#  Furthermore, Pope Gregory XIII required that *century years* (those ending in `00`) not be leap years unless they are evenly divisible by `400`.
!>  -#  The second correction by Pope Gregory XIII implies that there would be `3` fewer leap years every `400` years, which translates to `3` fewer days.
!>  -#  These corrections compensated for the accumulated errors of the Julian Calendar and made sure the corrections would remain in place for the foreseeable future.
!>  -#  Finally, Pope Gregory XIII decreed that the cutover date for the calendar should be \f$\ms{October 4, 1582}\f$.
!>  -#  The calendar with the above three fixes has become known as the **Gregorian Calendar** which is in widespread use in modern world.
!>  -#  Overall, Pope Gregory XIII striked **a total of `10` days** in the calendar starting \f$\ms{October 5, 1582}\f$ until and including \f$\ms{October 14, 1582}\f$.
!>  -#  On the day of the Gregorian Calendar cutover \f$\ms{October 4, 1582}\f$, the Catholic world (including Italy, Poland, Portugal, and Spain) was eager to follow the new calendar and all switched over to the Gregorian Calendar.
!>  -#  By the end of that year France, Holland, and part of Belgium made the switch.
!>  -#  The following year Austria, the rest of Belgium, and Catholic Germany fell in line.
!>  -#  And they were joined by Czechoslovakia and Catholic Switzerland in 1584, Hungary in 1587, and Transylvania in 1590.
!>  -#  The Protestant and Greek Orthodox countries did not immediately switch:
!>      -#  Germany switched piecemeal during the 1600s.
!>      -#  Denmark, Iceland, the rest of the Netherlands, Norway, and Protestant Switzerland switched in the year 1700.
!>      -#  Canada, Great Britain, Ireland, and the eastern US switched in 1752.
!>      -#  Japan switched in 1873 and Egypt in 1875.
!>      -#  Then between 1911 and 1923 Albania, Bulgaria, China, Estonia, Greece, Latvia, Lithuania, Romania, Russia, and Yugoslavia all switched over.
!>      -#  Finally Turkey switched in 1927.
!>  -#  The cutover was not simultaneous either within the United States mainland.
!>  -#  The calendar switch in the US depended on which country the specific territory was owned by:
!>      -#  Texas, Florida, California, Nevada, Arizona, and New Mexico all switched with Spain in 1582.
!>      -#  Mississippi switched with France in 1582.
!>      -#  The eastern seaboard switched with Great Britain in 1752.
!>      -#  Alaska switched in 1867 when it became part of the US.
!>  -#  In an attempt to have a gradual conversion, Sweden decided not to have leap years from 1700 to 1740.<br>
!>      This implied that there would be no jolt to the calendar.<br>
!>      But after skipping the leap year in 1700, they abandoned the plan.
!>      This put Sweden out of step with both the Julian and Gregorian calendars.<br>
!>      In 1712 Sweden reverted back to the Julian calendar by having 30 days in February that year to make up for the leap day that they missed in 1700.<br>
!>      Then in 1753 Sweden gave up and switched all-at-once to the Gregorian calendar.
!>  -#  Despite all switching, the Julian calendar is still used in parts of the world.<br>
!>      For example, it is used by Eastern Orthodox Church for calculating Easter and other feasts, by the Berber people in North Africa and on Mount Athos.<br>
!>      Also, Ethiopia uses the Alexandrian calendar which is based on the Julian calendar.
!>
!>  **The Gregorian Calendar Error**
!>
!>  -#  As mentioned above, the average Julian year is \f$\ms{365.25}\f$ days while the true length of a year is \f$\ms{365.2422}\f$ days.<br>
!>      The difference caused the Julian seasons to advance `1` day every `128` years.
!>  -#  The **average Gregorian year is \f$\ms{365.2425}\f$ days** which **causes the seasons to advance `1` day every `3333` years**.
!>  -#  The **Gregorian Calendar error** is about `1` day in `4000` years.<br>
!>  -#  The *Gregorian Calendar error* can be readily fixed by skipping the leap years in millennium years that are divisible by `4000`, requiring the month of February of year `4000` to have only `28` days.
!>  -#  However, there is currently no official fix introduced for the Gregorian Calendar.
!>  -#  The fix to the *Gregorian Calendar error* will make the average length of the modified Gregorian Calendar year \f$\ms{365.24225}\f$.
!>  -#  Therefore, the corrected Gregorian Calendar would have only a `1` day error in `20000` years (i.e., the seasons will shift by `1` day after `20000` years).
!>  -#  Any further fix to the modified Gregorian Calendar is likely effectively useless,
!>      since the error becomes comparable to the natural variations in the magic number \f$\ms{365.2422}\f$ (the number of days Earth takes to orbit the Sun).
!>
!>  **The zeroth year and the Gregorian Calendar**
!>
!>  A year zero does not exist in the **Anno Domini (AD)** calendar year system commonly used to number years in the Gregorian calendar
!>  (nor in its predecessor, the Julian calendar); in this system, the year \f$\ms{1 BC}\f$ is followed directly by year \f$\ms{AD 1}\f$.<br>
!>  However, <b>there is a year zero in both the astronomical year numbering system (where it coincides with the Julian year \f$\ms{1 BC}\f$),
!>  and the \f$\ms{ISO 8601:2004}\f$ international standard system</b>. This is the interchange standard for all calendar numbering systems.<br>
!>  The \f$\ms{ISO 8601:2004}\f$ convention for date and time explicitly requires year zero to coincide with the Gregorian year \f$\ms{1 BC}\f$.<br>
!>
!>  Gregorian year | ISO 8601 | Comments
!>  ---------------|----------|---------
!>  10000   BC     | −9999    | Beginning of the Holocene Era.
!>  9701    BC     | −9700    | End of the Pleistocene and beginning of the Holocene epoch.
!>  4714    BC     | −4713    | Epoch of the Julian day system: Julian day 0 starts at Greenwich noon on January 1, 4713 BC of the proleptic Julian calendar, which is November 24, 4714 BC in the proleptic Gregorian calendar.
!>  3761    BC     | −3760    | Beginning of the Anno Mundi calendar era in the Hebrew calendar.
!>  3102    BC     | −3101    | Beginning of the Kali Yuga in Hindu cosmology.
!>  2250    BC     | −2249    | Beginning of the Meghalayan age, the current and latest of the three stages in the Holocene era.
!>  45      BC     | −0044    | Introduction of the Julian calendar.
!>  1       BC     | +0000    | Year zero at ISO 8601.
!>  AD      1      | +0001    | Beginning of the Common Era and Anno Domini, from the estimate by Dionysius of the Incarnation of Jesus.
!>  AD      622    | +0622    | Migration of Muhammad from Mecca to Medina (Hegira), starting the Islamic calendar.
!>  AD      1582   | +1582    | Introduction of the Gregorian calendar.
!>  AD      1912   | +1912    | Epoch of the Juche[13] and Minguo calendars.
!>  AD      1950   | +1950    | Epoch of the Before Present dating scheme.
!>  AD      1960   | +1960    | UTC Epoch.
!>  AD      1970   | +1970    | Unix Epoch.
!>  AD      1993   | +1993    | Publication of the Holocene calendar.
!>  AD      2022   | +2022    | The future.
!>  AD      10000  | +10000   | The distant future.
!>
!>  The above information is partly based on the Wikipedia article [Year Zero](https://en.wikipedia.org/wiki/Year_zero)
!>  and the excellent historical review of the Julian Calendar by [Stephen P. Morse](https://stevemorse.org/).
!>
!>  \note
!>  If you ever need to generate uniformly-distributed random Gregorian calendar dates between two specified dates,<br>
!>  <ol>
!>      <li>    Convert the `lower` and `upper` limit dates to Julian Days via [getJulianDay()](@ref pm_dateTime::getJulianDay).<br>
!>      <li>    Generate a random Julian Day within the computed lower and upper Julian Day limits (`lowerJD`, `upperJD`) by calling [getUnifRand(lowerJD, upperJD)](@ref pm_distUnif::getUnifRand).<br>
!>              <ul>
!>                  <li>    Should the uniformly randomly generated Gregorian Calendar date be a UTC date and time (i.e., `zone = 0`),
!>                          then convert the uniformly randomly generated Julian Day (`randJD`) to the Gregorian Calendar date by calling [getDateTime(randJD)](@ref pm_dateTime::getDateTime).<br>
!>                  <li>    Should the uniformly randomly generated Gregorian Calendar date be local to a specific time zone `zone`,
!>                          then convert the uniformly randomly generated Julian Day (`randJD`) to the Gregorian Calendar date by calling [getDateTime(randJD, zone)](@ref pm_dateTime::getDateTime).<br>
!>              </ul>
!>  </ol>
!>  Here is an example:<br>
!>  \code{.F90}
!>
!>      use pm_kind, only: RKG => RK, IKG => IK
!>      integer(IKG) :: DateRand(8)
!>      DateRand = getDateTime(julianDay = getUnifRand(-300000._RKG, +300000._RKG))
!>      DateRand = getDateTime(julianDay = getUnifRand(getJulianDay(1_IK), getJulianDay())) ! uniform random date from the birth of Jesus until present.
!>      !
!>  \endcode
!>
!>  \see
!>  [M_time](https://github.com/urbanjost/M_time)<br>
!>  [datetime-fortran](https://github.com/wavebitscientific/datetime-fortran)<br>
!>
!>  \test
!>  [test_pm_dateTime](@ref test_pm_dateTime)
!>
!>  \todo
!>  A subroutine equivalent of the performance-critical functions with `allocatable`
!>  output (e.g., [getDateTime](@ref pm_dateTime::getDateTime)) should be added in the future.
!>
!>  \todo
!>  \phigh
!>  Most of the routines of this module are currently implemented for the default `integer` kind and the default `real64` real kind.<br>
!>  The choice of kinds is primarily dictated by the output of the Fortran intrinsic `date_and_time()`.<br>
!>  All routines should be extended to generic interfaces supporting multiple kinds in future.<br>
!>
!>  \final{pm_dateTime}
!>
!>  \author
!>  \AmirShahmoradi, January 30, 2021, 5:36 AM, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_dateTime

    use pm_kind, only: IK, LK, RK, SK

    implicit none

    !public
    !private :: queryDateTime

    character(*, SK), parameter :: MODULE_NAME = "@pm_dateTime"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    integer(IK) , parameter :: ZONE_MIN = -12_IK * 60_IK    !<   \public The scalar constant of default `integer` kind \IK representing the current minimum existing time zone value in the work in units of minutes.
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: ZONE_MIN
#endif

    integer(IK) , parameter :: ZONE_MAX = +14_IK * 60_IK    !<   \public The scalar constant of default `integer` kind \IK representing the current maximum existing time zone value in the work in units of minutes.
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: ZONE_MAX
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant vector of size `8` of type `integer` of default kind \IK containing origin of the Gregorian calendar
    !>  in the same format as returned by the `values` argument of the Fortran intrinsic `date_and_time()`.
    !>
    !>  \details
    !>  This date corresponds to the midnight of the first day of January of year `1` AD.
    !>
    !>  \see
    !>  [dateTimeInt_type](@ref pm_dateTime::dateTimeInt_type)<br>
    !>
    !>  \final{ORIGIN}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    integer(IK) , parameter :: ORIGIN(8) = [integer(IK) :: 1, 1, 1, 0, 0, 0, 0, 0]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: ORIGIN
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The `character` constant vector of shape `(0:6)` of length type parameter `9` containing
    !>  the names of the days of a week **assuming Sunday as the zeroth day of the week**.
    !>
    !>  \warning
    !>  Note that the order and index of this array is not ISO 8601 compliant where *the first day of the week is Monday*.<br>
    !>  To get the ISO 8601 compliant weekday names, see [WEEKDAY_NAME_ISO](@ref pm_dateTime::WEEKDAY_NAME_ISO)
    !>
    !>
    !>  \note
    !>  To generate a three-letters name of day abbreviations, simply slice the day name `(1:3)`.<br>
    !>  See examples below for usage.
    !>
    !>  \see
    !>  [getWeekDay](@ref pm_dateTime::getWeekDay)<br>
    !>  [getWeekDayISO](@ref pm_dateTime::getWeekDayISO)<br>
    !>  [WEEKDAY_NAME_ISO](@ref pm_dateTime::WEEKDAY_NAME_ISO)<br>
    !>
    !>  \example{WEEKDAY_NAME}
    !>  \include{lineno} example/pm_dateTime/WEEKDAY_NAME/main.F90
    !>  \compilef{WEEKDAY_NAME}
    !>  \output{WEEKDAY_NAME}
    !>  \include{lineno} example/pm_dateTime/WEEKDAY_NAME/main.out.F90
    !>
    !>  \final{WEEKDAY_NAME}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    character(9,SK), parameter :: WEEKDAY_NAME(0:6) =   [ "Sunday   " &
                                                        , "Monday   " &
                                                        , "Tuesday  " &
                                                        , "Wednesday" &
                                                        , "Thursday " &
                                                        , "Friday   " &
                                                        , "Saturday " &
                                                        ]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: WEEKDAY_NAME
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The `character` constant vector of size `7` of length type parameter `9` containing the names of the days of a week.
    !>
    !>  \details
    !>  Following ISO 8601 convention, **the first day of the week is Monday**.
    !>
    !>  \note
    !>  To generate a three-letters name of day abbreviations, simply slice the day name `(1:3)`.<br>
    !>  See examples below for usage.
    !>
    !>  \see
    !>  [getWeekDay](@ref pm_dateTime::getWeekDay)<br>
    !>  [getWeekDayISO](@ref pm_dateTime::getWeekDayISO)<br>
    !>  [WEEKDAY_NAME](@ref pm_dateTime::WEEKDAY_NAME)<br>
    !>
    !>  \example{WEEKDAY_NAME_ISO}
    !>  \include{lineno} example/pm_dateTime/WEEKDAY_NAME_ISO/main.F90
    !>  \compilef{WEEKDAY_NAME_ISO}
    !>  \output{WEEKDAY_NAME_ISO}
    !>  \include{lineno} example/pm_dateTime/WEEKDAY_NAME_ISO/main.out.F90
    !>
    !>  \final{WEEKDAY_NAME_ISO}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    character(9,SK), parameter :: WEEKDAY_NAME_ISO(1:7) =   [ "Monday   " &
                                                            , "Tuesday  " &
                                                            , "Wednesday" &
                                                            , "Thursday " &
                                                            , "Friday   " &
                                                            , "Saturday " &
                                                            , "Sunday   " &
                                                            ]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: WEEKDAY_NAME_ISO
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The `character` constant vector of size `12` of length type-parameter `9`,
    !>  containing full names of the months of the Gregorian calendar.
    !>
    !>  \see
    !>  [getMonth](@ref pm_dateTime::getMonth)<br>
    !>
    !>  \example{MONTH_NAME}
    !>  \include{lineno} example/pm_dateTime/MONTH_NAME/main.F90
    !>  \compilef{MONTH_NAME}
    !>  \output{MONTH_NAME}
    !>  \include{lineno} example/pm_dateTime/MONTH_NAME/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{MONTH_NAME}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    character(*, SK), parameter :: MONTH_NAME(12) = [ "January  " &
                                                    , "February " &
                                                    , "March    " &
                                                    , "April    " &
                                                    , "May      " &
                                                    , "June     " &
                                                    , "July     " &
                                                    , "August   " &
                                                    , "September" &
                                                    , "October  " &
                                                    , "November " &
                                                    , "December " &
                                                    ]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: MONTH_NAME
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The `integer` constant vector of size `12` containing number of days in each month of a common (non-leap) year of the Gregorian calendar.
    !>
    !>  \see
    !>  [getCountDays](@ref pm_dateTime::getCountDays)<br>
    !>
    !>  \example{DAYS_OF_MONTH}
    !>  \include{lineno} example/pm_dateTime/DAYS_OF_MONTH/main.F90
    !>  \compilef{DAYS_OF_MONTH}
    !>  \output{DAYS_OF_MONTH}
    !>  \include{lineno} example/pm_dateTime/DAYS_OF_MONTH/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{DAYS_OF_MONTH}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    integer(IK) , parameter :: DAYS_OF_MONTH(12) = [integer(IK) :: 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: DAYS_OF_MONTH
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The `integer` constant vector of size `12` containing number of days in each month of a leap year of the Gregorian calendar.
    !>
    !>  \details
    !>  Only the month of February becomes one day longer in the leap years.
    !>
    !>  \see
    !>  [isLeapYear](@ref pm_dateTime::isLeapYear)<br>
    !>  [getCountDays](@ref pm_dateTime::getCountDays)<br>
    !>
    !>  \example{DAYS_OF_MONTH_LEAP}
    !>  \include{lineno} example/pm_dateTime/DAYS_OF_MONTH_LEAP/main.F90
    !>  \compilef{DAYS_OF_MONTH_LEAP}
    !>  \output{DAYS_OF_MONTH_LEAP}
    !>  \include{lineno} example/pm_dateTime/DAYS_OF_MONTH_LEAP/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{DAYS_OF_MONTH_LEAP}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    integer(IK) , parameter :: DAYS_OF_MONTH_LEAP(12) = [integer(IK) :: 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: DAYS_OF_MONTH_LEAP
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant scalar of type `integer` of kind \IK containing the number of seconds per day.
    !>
    !>  \see
    !>  [MINUTES_PER_DAY](@ref pm_dateTime::MINUTES_PER_DAY)<br>
    !>
    !>  \final{SECONDS_PER_DAY}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    integer(IK) , parameter :: SECONDS_PER_DAY = 86400_IK
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: SECONDS_PER_DAY
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant scalar of type `integer` of kind \IK containing the number of minutes per day.
    !>
    !>  \see
    !>  [SECONDS_PER_DAY](@ref pm_dateTime::SECONDS_PER_DAY)<br>
    !>
    !>  \final{MINUTES_PER_DAY}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    integer(IK) , parameter :: MINUTES_PER_DAY = 1440_IK
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: MINUTES_PER_DAY
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant scalar of type `real` of default kind \RK containing the **average** number of days per month.
    !>
    !>  \details
    !>  This number takes into account the different length of the `February` month in leap years.
    !>
    !>  \see
    !>  [MEAN_DAYS_PER_MONTH](@ref pm_dateTime::MEAN_DAYS_PER_MONTH)<br>
    !>  [MEAN_WEEKS_PER_MONTH](@ref pm_dateTime::MEAN_WEEKS_PER_MONTH)<br>
    !>  [MEAN_HOURS_PER_MONTH](@ref pm_dateTime::MEAN_HOURS_PER_MONTH)<br>
    !>  [MEAN_MINUTES_PER_MONTH](@ref pm_dateTime::MEAN_MINUTES_PER_MONTH)<br>
    !>  [MEAN_SECONDS_PER_MONTH](@ref pm_dateTime::MEAN_SECONDS_PER_MONTH)<br>
    !>  [MEAN_DAYS_PER_YEAR](@ref pm_dateTime::MEAN_DAYS_PER_YEAR)<br>
    !>  [MEAN_WEEKS_PER_YEAR](@ref pm_dateTime::MEAN_WEEKS_PER_YEAR)<br>
    !>  [MEAN_HOURS_PER_YEAR](@ref pm_dateTime::MEAN_HOURS_PER_YEAR)<br>
    !>  [MEAN_MINUTES_PER_YEAR](@ref pm_dateTime::MEAN_MINUTES_PER_YEAR)<br>
    !>  [MEAN_SECONDS_PER_YEAR](@ref pm_dateTime::MEAN_SECONDS_PER_YEAR)<br>
    !>
    !>  \final{MEAN_DAYS_PER_MONTH}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    real(RK)    , parameter :: MEAN_DAYS_PER_MONTH = 0.25_RK * (0.25_RK * sum(DAYS_OF_MONTH) + sum(DAYS_OF_MONTH_LEAP) / 12._RK) ! = 30.4375_RK
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: MEAN_DAYS_PER_MONTH
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant scalar of type `real` of default kind \RK containing the *approximate (rounded)* **average** number of weeks per month.
    !>
    !>  \details
    !>  This number takes into account the different length of the `February` month in leap years.
    !>
    !>  \see
    !>  [MEAN_DAYS_PER_MONTH](@ref pm_dateTime::MEAN_DAYS_PER_MONTH)<br>
    !>  [MEAN_WEEKS_PER_MONTH](@ref pm_dateTime::MEAN_WEEKS_PER_MONTH)<br>
    !>  [MEAN_HOURS_PER_MONTH](@ref pm_dateTime::MEAN_HOURS_PER_MONTH)<br>
    !>  [MEAN_MINUTES_PER_MONTH](@ref pm_dateTime::MEAN_MINUTES_PER_MONTH)<br>
    !>  [MEAN_SECONDS_PER_MONTH](@ref pm_dateTime::MEAN_SECONDS_PER_MONTH)<br>
    !>  [MEAN_DAYS_PER_YEAR](@ref pm_dateTime::MEAN_DAYS_PER_YEAR)<br>
    !>  [MEAN_WEEKS_PER_YEAR](@ref pm_dateTime::MEAN_WEEKS_PER_YEAR)<br>
    !>  [MEAN_HOURS_PER_YEAR](@ref pm_dateTime::MEAN_HOURS_PER_YEAR)<br>
    !>  [MEAN_MINUTES_PER_YEAR](@ref pm_dateTime::MEAN_MINUTES_PER_YEAR)<br>
    !>  [MEAN_SECONDS_PER_YEAR](@ref pm_dateTime::MEAN_SECONDS_PER_YEAR)<br>
    !>
    !>  \final{MEAN_WEEKS_PER_MONTH}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    real(RK)    , parameter :: MEAN_WEEKS_PER_MONTH = MEAN_DAYS_PER_MONTH / 7._RK ! = 4.34821428571428571428571428571429_RK
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: MEAN_WEEKS_PER_MONTH
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant scalar of type `real` of default kind \RK containing the **average** number of hours per month.
    !>
    !>  \details
    !>  This number takes into account the different length of the `February` month in leap years.
    !>
    !>  \see
    !>  [MEAN_DAYS_PER_MONTH](@ref pm_dateTime::MEAN_DAYS_PER_MONTH)<br>
    !>  [MEAN_WEEKS_PER_MONTH](@ref pm_dateTime::MEAN_WEEKS_PER_MONTH)<br>
    !>  [MEAN_HOURS_PER_MONTH](@ref pm_dateTime::MEAN_HOURS_PER_MONTH)<br>
    !>  [MEAN_MINUTES_PER_MONTH](@ref pm_dateTime::MEAN_MINUTES_PER_MONTH)<br>
    !>  [MEAN_SECONDS_PER_MONTH](@ref pm_dateTime::MEAN_SECONDS_PER_MONTH)<br>
    !>  [MEAN_DAYS_PER_YEAR](@ref pm_dateTime::MEAN_DAYS_PER_YEAR)<br>
    !>  [MEAN_WEEKS_PER_YEAR](@ref pm_dateTime::MEAN_WEEKS_PER_YEAR)<br>
    !>  [MEAN_HOURS_PER_YEAR](@ref pm_dateTime::MEAN_HOURS_PER_YEAR)<br>
    !>  [MEAN_MINUTES_PER_YEAR](@ref pm_dateTime::MEAN_MINUTES_PER_YEAR)<br>
    !>  [MEAN_SECONDS_PER_YEAR](@ref pm_dateTime::MEAN_SECONDS_PER_YEAR)<br>
    !>
    !>  \final{MEAN_HOURS_PER_MONTH}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    real(RK)    , parameter :: MEAN_HOURS_PER_MONTH = 24_IK * MEAN_DAYS_PER_MONTH ! = 730.5_RK
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: MEAN_HOURS_PER_MONTH
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant scalar of type `real` of default kind \RK containing the **average** number of minutes per month.
    !>
    !>  \details
    !>  This number takes into account the different length of the `February` month in leap years.
    !>
    !>  \see
    !>  [MEAN_DAYS_PER_MONTH](@ref pm_dateTime::MEAN_DAYS_PER_MONTH)<br>
    !>  [MEAN_WEEKS_PER_MONTH](@ref pm_dateTime::MEAN_WEEKS_PER_MONTH)<br>
    !>  [MEAN_HOURS_PER_MONTH](@ref pm_dateTime::MEAN_HOURS_PER_MONTH)<br>
    !>  [MEAN_MINUTES_PER_MONTH](@ref pm_dateTime::MEAN_MINUTES_PER_MONTH)<br>
    !>  [MEAN_SECONDS_PER_MONTH](@ref pm_dateTime::MEAN_SECONDS_PER_MONTH)<br>
    !>  [MEAN_DAYS_PER_YEAR](@ref pm_dateTime::MEAN_DAYS_PER_YEAR)<br>
    !>  [MEAN_WEEKS_PER_YEAR](@ref pm_dateTime::MEAN_WEEKS_PER_YEAR)<br>
    !>  [MEAN_HOURS_PER_YEAR](@ref pm_dateTime::MEAN_HOURS_PER_YEAR)<br>
    !>  [MEAN_MINUTES_PER_YEAR](@ref pm_dateTime::MEAN_MINUTES_PER_YEAR)<br>
    !>  [MEAN_SECONDS_PER_YEAR](@ref pm_dateTime::MEAN_SECONDS_PER_YEAR)<br>
    !>
    !>  \final{MEAN_MINUTES_PER_MONTH}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    real(RK)    , parameter :: MEAN_MINUTES_PER_MONTH = 60_IK * MEAN_HOURS_PER_MONTH ! = 43830._RK
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: MEAN_MINUTES_PER_MONTH
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant scalar of type `real` of default kind \RK containing the **average** number of seconds per month.
    !>
    !>  \details
    !>  This number takes into account the different length of the `February` month in leap years.
    !>
    !>  \see
    !>  [MEAN_DAYS_PER_MONTH](@ref pm_dateTime::MEAN_DAYS_PER_MONTH)<br>
    !>  [MEAN_WEEKS_PER_MONTH](@ref pm_dateTime::MEAN_WEEKS_PER_MONTH)<br>
    !>  [MEAN_HOURS_PER_MONTH](@ref pm_dateTime::MEAN_HOURS_PER_MONTH)<br>
    !>  [MEAN_MINUTES_PER_MONTH](@ref pm_dateTime::MEAN_MINUTES_PER_MONTH)<br>
    !>  [MEAN_SECONDS_PER_MONTH](@ref pm_dateTime::MEAN_SECONDS_PER_MONTH)<br>
    !>  [MEAN_DAYS_PER_YEAR](@ref pm_dateTime::MEAN_DAYS_PER_YEAR)<br>
    !>  [MEAN_WEEKS_PER_YEAR](@ref pm_dateTime::MEAN_WEEKS_PER_YEAR)<br>
    !>  [MEAN_HOURS_PER_YEAR](@ref pm_dateTime::MEAN_HOURS_PER_YEAR)<br>
    !>  [MEAN_MINUTES_PER_YEAR](@ref pm_dateTime::MEAN_MINUTES_PER_YEAR)<br>
    !>  [MEAN_SECONDS_PER_YEAR](@ref pm_dateTime::MEAN_SECONDS_PER_YEAR)<br>
    !>
    !>  \final{MEAN_SECONDS_PER_MONTH}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    real(RK)    , parameter :: MEAN_SECONDS_PER_MONTH = 60_IK * MEAN_MINUTES_PER_MONTH ! = 2629800._RK
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: MEAN_SECONDS_PER_MONTH
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant scalar of type `real` of default kind \RK containing the **average** number of days per year.
    !>
    !>  \details
    !>  This number takes into account the leap years.
    !>
    !>  \see
    !>  [MEAN_DAYS_PER_MONTH](@ref pm_dateTime::MEAN_DAYS_PER_MONTH)<br>
    !>  [MEAN_WEEKS_PER_MONTH](@ref pm_dateTime::MEAN_WEEKS_PER_MONTH)<br>
    !>  [MEAN_HOURS_PER_MONTH](@ref pm_dateTime::MEAN_HOURS_PER_MONTH)<br>
    !>  [MEAN_MINUTES_PER_MONTH](@ref pm_dateTime::MEAN_MINUTES_PER_MONTH)<br>
    !>  [MEAN_SECONDS_PER_MONTH](@ref pm_dateTime::MEAN_SECONDS_PER_MONTH)<br>
    !>  [MEAN_DAYS_PER_YEAR](@ref pm_dateTime::MEAN_DAYS_PER_YEAR)<br>
    !>  [MEAN_WEEKS_PER_YEAR](@ref pm_dateTime::MEAN_WEEKS_PER_YEAR)<br>
    !>  [MEAN_HOURS_PER_YEAR](@ref pm_dateTime::MEAN_HOURS_PER_YEAR)<br>
    !>  [MEAN_MINUTES_PER_YEAR](@ref pm_dateTime::MEAN_MINUTES_PER_YEAR)<br>
    !>  [MEAN_SECONDS_PER_YEAR](@ref pm_dateTime::MEAN_SECONDS_PER_YEAR)<br>
    !>
    !>  \final{MEAN_DAYS_PER_YEAR}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    real(RK)    , parameter :: MEAN_DAYS_PER_YEAR = 0.25_RK * (3 * 365 + 366) ! = 365.25_RK
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: MEAN_DAYS_PER_YEAR
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant scalar of type `real` of default kind \RK containing the *approximate (rounded)* **average** number of weeks per year.
    !>
    !>  \details
    !>  This number takes into account the different length of the `February` month in leap years.
    !>
    !>  \see
    !>  [MEAN_DAYS_PER_MONTH](@ref pm_dateTime::MEAN_DAYS_PER_MONTH)<br>
    !>  [MEAN_WEEKS_PER_MONTH](@ref pm_dateTime::MEAN_WEEKS_PER_MONTH)<br>
    !>  [MEAN_HOURS_PER_MONTH](@ref pm_dateTime::MEAN_HOURS_PER_MONTH)<br>
    !>  [MEAN_MINUTES_PER_MONTH](@ref pm_dateTime::MEAN_MINUTES_PER_MONTH)<br>
    !>  [MEAN_SECONDS_PER_MONTH](@ref pm_dateTime::MEAN_SECONDS_PER_MONTH)<br>
    !>  [MEAN_DAYS_PER_YEAR](@ref pm_dateTime::MEAN_DAYS_PER_YEAR)<br>
    !>  [MEAN_WEEKS_PER_YEAR](@ref pm_dateTime::MEAN_WEEKS_PER_YEAR)<br>
    !>  [MEAN_HOURS_PER_YEAR](@ref pm_dateTime::MEAN_HOURS_PER_YEAR)<br>
    !>  [MEAN_MINUTES_PER_YEAR](@ref pm_dateTime::MEAN_MINUTES_PER_YEAR)<br>
    !>  [MEAN_SECONDS_PER_YEAR](@ref pm_dateTime::MEAN_SECONDS_PER_YEAR)<br>
    !>
    !>  \final{MEAN_WEEKS_PER_YEAR}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    real(RK)    , parameter :: MEAN_WEEKS_PER_YEAR = MEAN_DAYS_PER_YEAR / 7._RK ! = 52.1785714285714285714285714285714_RK
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: MEAN_WEEKS_PER_YEAR
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant scalar of type `real` of default kind \RK containing the **average** number of hours per year.
    !>
    !>  \details
    !>  This number takes into account the different length of the `February` month in leap years.
    !>
    !>  \see
    !>  [MEAN_DAYS_PER_MONTH](@ref pm_dateTime::MEAN_DAYS_PER_MONTH)<br>
    !>  [MEAN_WEEKS_PER_MONTH](@ref pm_dateTime::MEAN_WEEKS_PER_MONTH)<br>
    !>  [MEAN_HOURS_PER_MONTH](@ref pm_dateTime::MEAN_HOURS_PER_MONTH)<br>
    !>  [MEAN_MINUTES_PER_MONTH](@ref pm_dateTime::MEAN_MINUTES_PER_MONTH)<br>
    !>  [MEAN_SECONDS_PER_MONTH](@ref pm_dateTime::MEAN_SECONDS_PER_MONTH)<br>
    !>  [MEAN_DAYS_PER_YEAR](@ref pm_dateTime::MEAN_DAYS_PER_YEAR)<br>
    !>  [MEAN_WEEKS_PER_YEAR](@ref pm_dateTime::MEAN_WEEKS_PER_YEAR)<br>
    !>  [MEAN_HOURS_PER_YEAR](@ref pm_dateTime::MEAN_HOURS_PER_YEAR)<br>
    !>  [MEAN_MINUTES_PER_YEAR](@ref pm_dateTime::MEAN_MINUTES_PER_YEAR)<br>
    !>  [MEAN_SECONDS_PER_YEAR](@ref pm_dateTime::MEAN_SECONDS_PER_YEAR)<br>
    !>
    !>  \final{MEAN_HOURS_PER_YEAR}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    real(RK)    , parameter :: MEAN_HOURS_PER_YEAR = 24_IK * MEAN_DAYS_PER_YEAR ! = 8766._RK
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: MEAN_HOURS_PER_YEAR
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant scalar of type `real` of default kind \RK containing the **average** number of minutes per year.
    !>
    !>  \details
    !>  This number takes into account the different length of the `February` month in leap years.
    !>
    !>  \see
    !>  [MEAN_DAYS_PER_MONTH](@ref pm_dateTime::MEAN_DAYS_PER_MONTH)<br>
    !>  [MEAN_WEEKS_PER_MONTH](@ref pm_dateTime::MEAN_WEEKS_PER_MONTH)<br>
    !>  [MEAN_HOURS_PER_MONTH](@ref pm_dateTime::MEAN_HOURS_PER_MONTH)<br>
    !>  [MEAN_MINUTES_PER_MONTH](@ref pm_dateTime::MEAN_MINUTES_PER_MONTH)<br>
    !>  [MEAN_SECONDS_PER_MONTH](@ref pm_dateTime::MEAN_SECONDS_PER_MONTH)<br>
    !>  [MEAN_DAYS_PER_YEAR](@ref pm_dateTime::MEAN_DAYS_PER_YEAR)<br>
    !>  [MEAN_WEEKS_PER_YEAR](@ref pm_dateTime::MEAN_WEEKS_PER_YEAR)<br>
    !>  [MEAN_HOURS_PER_YEAR](@ref pm_dateTime::MEAN_HOURS_PER_YEAR)<br>
    !>  [MEAN_MINUTES_PER_YEAR](@ref pm_dateTime::MEAN_MINUTES_PER_YEAR)<br>
    !>  [MEAN_SECONDS_PER_YEAR](@ref pm_dateTime::MEAN_SECONDS_PER_YEAR)<br>
    !>
    !>  \final{MEAN_MINUTES_PER_YEAR}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    real(RK)    , parameter :: MEAN_MINUTES_PER_YEAR = 60_IK * MEAN_HOURS_PER_YEAR ! = 525960._RK
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: MEAN_MINUTES_PER_YEAR
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant scalar of type `real` of default kind \RK containing the **average** number of seconds per year.
    !>
    !>  \details
    !>  This number takes into account the different length of the `February` month in leap years.
    !>
    !>  \see
    !>  [MEAN_DAYS_PER_MONTH](@ref pm_dateTime::MEAN_DAYS_PER_MONTH)<br>
    !>  [MEAN_WEEKS_PER_MONTH](@ref pm_dateTime::MEAN_WEEKS_PER_MONTH)<br>
    !>  [MEAN_HOURS_PER_MONTH](@ref pm_dateTime::MEAN_HOURS_PER_MONTH)<br>
    !>  [MEAN_MINUTES_PER_MONTH](@ref pm_dateTime::MEAN_MINUTES_PER_MONTH)<br>
    !>  [MEAN_SECONDS_PER_MONTH](@ref pm_dateTime::MEAN_SECONDS_PER_MONTH)<br>
    !>  [MEAN_DAYS_PER_YEAR](@ref pm_dateTime::MEAN_DAYS_PER_YEAR)<br>
    !>  [MEAN_WEEKS_PER_YEAR](@ref pm_dateTime::MEAN_WEEKS_PER_YEAR)<br>
    !>  [MEAN_HOURS_PER_YEAR](@ref pm_dateTime::MEAN_HOURS_PER_YEAR)<br>
    !>  [MEAN_MINUTES_PER_YEAR](@ref pm_dateTime::MEAN_MINUTES_PER_YEAR)<br>
    !>  [MEAN_SECONDS_PER_YEAR](@ref pm_dateTime::MEAN_SECONDS_PER_YEAR)<br>
    !>
    !>  \final{MEAN_SECONDS_PER_YEAR}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    real(RK)    , parameter :: MEAN_SECONDS_PER_YEAR = 60_IK * MEAN_MINUTES_PER_YEAR ! = 31557600._RK
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: MEAN_SECONDS_PER_YEAR
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type containing the components of a numeric date and time in the Gregorian calendar.
    !>
    !>  \details
    !>  This class is a simple container for the values returned by the Fortran intrinsic `date_and_time()`.
    !>
    !>  \interface{dateTimeInt_type}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: dateTimeInt_type
    !>      type(dateTimeInt_type) :: dateTimeInt
    !>
    !>      dateTimeInt = dateTimeInt_type( year = year &
    !>                                    , month = month &
    !>                                    , day = day &
    !>                                    , zone = zone &
    !>                                    , hour = hour &
    !>                                    , minute = minute &
    !>                                    , second = second &
    !>                                    , millisecond = millisecond &
    !>                                    )
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [dateTimeInt_typer](@ref pm_dateTime::dateTimeInt_typer) (class constructor)<br>
    !>  [dateTimeStr_type](@ref pm_dateTime::dateTimeStr_type)<br>
    !>
    !>  \example{dateTimeInt_type}
    !>  \include{lineno} example/pm_dateTime/dateTimeInt_type/main.F90
    !>  \compilef{dateTimeInt_type}
    !>  \output{dateTimeInt_type}
    !>  \include{lineno} example/pm_dateTime/dateTimeInt_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{dateTimeInt_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    type                    :: dateTimeInt_type
       !integer(IK), private:: values(8)    = 0_IK  !   \private    The vector `integer` of default kind \IK of length `8` containing `[year, month, day, zone, hour, minute, millisecond]`.
        integer(IK)         :: year         = 1_IK  !<  \public     The scalar `integer` of default kind \IK containing the year of the Gregorian calendar.
        integer(IK)         :: month        = 1_IK  !<  \public     The scalar `integer` of default kind \IK containing the month of the year.
        integer(IK)         :: day          = 1_IK  !<  \public     The scalar `integer` of default kind \IK containing the day of the month.
        integer(IK)         :: zone         = 0_IK  !<  \public     The scalar `integer` of default kind \IK containing the time **difference in minutes**
                                                    !!              with respect to the [UTC](https://www.timeanddate.com/worldclock/timezone/utc).
        integer(IK)         :: hour         = 0_IK  !<  \public     The scalar `integer` of default kind \IK containing the hour of the day.
        integer(IK)         :: minute       = 0_IK  !<  \public     The scalar `integer` of default kind \IK containing the minute of the hour.
        integer(IK)         :: second       = 0_IK  !<  \public     The scalar `integer` of default kind \IK containing the second of the minute.
        integer(IK)         :: millisecond  = 0_IK  !<  \public     The scalar `integer` of default kind \IK containing the milliseconds of the second.
    contains
        procedure, pass     :: getValues => getDateTimeIntValues
    end type

    interface dateTimeInt_type
        module procedure :: dateTimeInt_typer
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type containing the string components that contain numeric date and time in the Gregorian calendar.
    !>
    !>  \details
    !>  This class is a simple container for the values returned by the Fortran intrinsic `date_and_time()`.
    !>
    !>  \interface{dateTimeStr_type}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: dateTimeStr_type
    !>      type(dateTimeStr_type) :: dateTimeStr
    !>
    !>      dateTimeStr = dateTimeStr_type( year = year &
    !>                                    , month = month &
    !>                                    , day = day &
    !>                                    , zone = zone &
    !>                                    , hour = hour &
    !>                                    , minute = minute &
    !>                                    , second = second &
    !>                                    , millisecond = millisecond &
    !>                                    )
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [dateTimeStr_typer](@ref pm_dateTime::dateTimeStr_typer) (class constructor)<br>
    !>  [dateTimeInt_type](@ref pm_dateTime::dateTimeInt_type)<br>
    !>
    !>  \example{dateTimeStr_type}
    !>  \include{lineno} example/pm_dateTime/dateTimeStr_type/main.F90
    !>  \compilef{dateTimeStr_type}
    !>  \output{dateTimeStr_type}
    !>  \include{lineno} example/pm_dateTime/dateTimeStr_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{dateTimeStr_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    type                            :: dateTimeStr_type
       !character(8, SK), private   :: date         = "00000000"    !<  \private    The scalar `character` of default kind \SK of length `8`    containing the Gregorian calendar date in the form `yyyymmdd`.
       !character(10,SK), private   :: time         = "0000000000"  !<  \private    The scalar `character` of default kind \SK of length `10`   containing the Gregorian calendar date in the form `hhmmss.sss`.
       !character(4, SK)            :: century      = "00"          !<  \public     The scalar `character` of default kind \SK of length `2`    containing the century of the Gregorian calendar in the form `cc`.
        character(4, SK)            :: year         = "0001"        !<  \public     The scalar `character` of default kind \SK of length `4`    containing the year of the Gregorian calendar in the form `yyyy`.
        character(2, SK)            :: month        = "01"          !<  \public     The scalar `character` of default kind \SK of length `2`    containing the month of the year in the form `mm`.
        character(2, SK)            :: day          = "01"          !<  \public     The scalar `character` of default kind \SK of length `2`    containing the day of the month in the form `dd`.
        character(5, SK)            :: zone         = "+0000"       !<  \public     The scalar `character` of default kind \SK of length `5`    containing the time difference between
                                                                    !!              local time and UTC (also known as Greenwich mean Time) in the form `Shhmm`,
                                                                    !!              corresponding to *sign*, *hours*, and *minutes*. For example, `-0500` (New York).
        character(2, SK)            :: hour         = "00"          !<  \public     The scalar `character` of default kind \SK of length `2`    containing the hour of the day in the form `hh`.
        character(2, SK)            :: minute       = "00"          !<  \public     The scalar `character` of default kind \SK of length `2`    containing the minute of the hour in the form `mm`.
        character(2, SK)            :: second       = "00"          !<  \public     The scalar `character` of default kind \SK of length `2`    containing the second of the minute in the form `ss`.
        character(3, SK)            :: millisecond  = "000"         !<  \public     The scalar `character` of default kind \SK of length `3`    containing the milliseconds of the second in the form `sss`.
   !contains
   !    procedure, pass             :: getDate => getDateTimeStrDate
   !    procedure, pass             :: time => getDateTimeStrTime
    end type

    interface dateTimeStr_type
        module procedure :: dateTimeStr_typer
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! >  \brief
    ! >  This is the derived type for containing a given Gregorian calendar date and time in both numeric and string formats.
    ! >
    ! >  \details
    ! >  This class is a simple container for the values returned by the Fortran intrinsic `date_and_time()`.
    ! >
    ! >  \interface{DateTime_type}
    ! >  \code{.F90}
    ! >
    ! >      use pm_dateTime, only: DateTime_type
    ! >      type(DateTime_type) :: DateTime
    ! >
    ! >      DateTime = DateTime_type( Int = dateTimeInt_type() &
    ! >                              , Str = dateTimeStr_type() &
    ! >                              )
    ! >
    ! >  \endcode
    ! >
    ! >  \see
    ! >  [constructDateTime](@ref pm_dateTime::constructDateTime) (class constructor)<br>
    ! >  [dateTimeInt_type](@ref pm_dateTime::dateTimeInt_type)<br>
    ! >  [dateTimeStr_type](@ref pm_dateTime::dateTimeStr_type)<br>
    ! >
    ! >  \example{DateTime_type}
    ! >  \include{lineno} example/pm_dateTime/DateTime_type/main.F90
    ! >  \compile{DateTime_type}
    ! >  \output{DateTime_type}
    ! >  \include{lineno} example/pm_dateTime/DateTime_type/main.out.F90
    ! >
    ! >  \test
    ! >  [test_pm_dateTime](@ref test_pm_dateTime)
    ! >
    ! >  \final{DateTime_type}
    ! >
    ! >  \author
    ! >  Amir Shahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    !type                        :: DateTime_type
    !    type(dateTimeInt_type)  :: Int  !<  \public The scalar object of type [dateTimeInt_type](@ref pm_dateTime::dateTimeInt_type) containing the Gregorian calendar date and time in integer format.
    !    type(dateTimeStr_type)  :: Str  !<  \public The scalar object of type [dateTimeStr_type](@ref pm_dateTime::dateTimeStr_type) containing the Gregorian calendar date and time in string format.
    !end type DateTime_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating object parameters containing a list of time zones and their representative abbreviations.<br>
    !>
    !>  \details
    !>  The time zone abbreviation is difficult to infer and generally requires communication with the operating system, which may or may not have it.<br>
    !>  The current implementation of this generic interface relies on a predefined internal list of time zone abbreviations to convert the specified `zone` **in units of minutes** to an abbreviation.<br>
    !>  Where there are multiple time zone abbreviations available for a single time zone, only a single representative abbreviation is kept in the list.<br>
    !>  The following is the current list of internal time zone abbreviations used in the ParaMonte library.<br>
    !>
    !>  Timezone                    |   Abbreviation        |   Timezone name
    !>  ----------------------------|-----------------------|----------------------------------------
    !>  \f$-60 \times 12     \f$    |   \f$\ms{IDLW }\f$    |   International Day Line West time zone
    !>  \f$-60 \times 11     \f$    |   \f$\ms{SST  }\f$    |   Samoa Standard Time
    !>  \f$-60 \times 10     \f$    |   \f$\ms{HST  }\f$    |   Hawaii–Aleutian Standard Time
    !>  \f$-60 \times 9 - 30 \f$    |   \f$\ms{MIT  }\f$    |   Marquesas Islands Time
    !>  \f$-60 \times 9      \f$    |   \f$\ms{AKST }\f$    |   Alaska Standard Time
    !>  \f$-60 \times 8      \f$    |   \f$\ms{PST  }\f$    |   Pacific Standard Time (North America)
    !>  \f$-60 \times 7      \f$    |   \f$\ms{MST  }\f$    |   Mountain Standard Time (North America)
    !>  \f$-60 \times 6      \f$    |   \f$\ms{CST  }\f$    |   Central Standard Time (North America)
    !>  \f$-60 \times 5      \f$    |   \f$\ms{EST  }\f$    |   Eastern Standard Time (North America)
    !>  \f$-60 \times 3 - 30 \f$    |   \f$\ms{NST  }\f$    |   Newfoundland Standard Time
    !>  \f$-60 \times 3      \f$    |   \f$\ms{UYT  }\f$    |   Uruguay Standard Time
    !>  \f$-60 \times 2 - 30 \f$    |   \f$\ms{NDT  }\f$    |   Newfoundland Daylight Time
    !>  \f$-60 \times 2      \f$    |   \f$\ms{UYST }\f$    |   Uruguay Summer Time
    !>  \f$-60 \times 1      \f$    |   \f$\ms{EGT  }\f$    |   Eastern Greenland Time
    !>  \f$+60 \times 0      \f$    |   \f$\ms{UTC  }\f$    |   Coordinated Universal Time
    !>  \f$+60 \times 1      \f$    |   \f$\ms{CET  }\f$    |   Central European Time
    !>  \f$+60 \times 2      \f$    |   \f$\ms{EET  }\f$    |   Eastern European Time
    !>  \f$+60 \times 3      \f$    |   \f$\ms{AST  }\f$    |   Arabia Standard Time
    !>  \f$+60 \times 3 + 30 \f$    |   \f$\ms{IRST }\f$    |   Iran Standard Time
    !>  \f$+60 \times 4      \f$    |   \f$\ms{GET  }\f$    |   Georgia Standard Time
    !>  \f$+60 \times 4 + 30 \f$    |   \f$\ms{AFT  }\f$    |   Afghanistan Time
    !>  \f$+60 \times 5      \f$    |   \f$\ms{PKT  }\f$    |   Pakistan Standard Time
    !>  \f$+60 \times 5 + 30 \f$    |   \f$\ms{IST  }\f$    |   Indian Standard Time
    !>  \f$+60 \times 5 + 45 \f$    |   \f$\ms{NPT  }\f$    |   Nepal Time
    !>  \f$+60 \times 6      \f$    |   \f$\ms{BST  }\f$    |   Bangladesh Standard Time
    !>  \f$+60 \times 6 + 30 \f$    |   \f$\ms{MMT  }\f$    |   Myanmar Standard Time
    !>  \f$+60 \times 7      \f$    |   \f$\ms{THA  }\f$    |   Thailand Standard Time
    !>  \f$+60 \times 8      \f$    |   \f$\ms{SST  }\f$    |   Singapore Standard Time
    !>  \f$+60 \times 8 + 45 \f$    |   \f$\ms{CWST }\f$    |   Central Western Standard Time (Australia)
    !>  \f$+60 \times 9      \f$    |   \f$\ms{JST  }\f$    |   Japan Standard Time
    !>  \f$+60 \times 9 + 30 \f$    |   \f$\ms{ACST }\f$    |   Australian Central Standard Time
    !>  \f$+60 \times 10     \f$    |   \f$\ms{AEST }\f$    |   Australian Eastern Standard Time
    !>  \f$+60 \times 10 + 30\f$    |   \f$\ms{LHST }\f$    |   Lord Howe Standard Time
    !>  \f$+60 \times 11     \f$    |   \f$\ms{PONT }\f$    |   Pohnpei Standard Time
    !>  \f$+60 \times 12     \f$    |   \f$\ms{NZST }\f$    |   New Zealand Standard Time
    !>  \f$+60 \times 12 + 45\f$    |   \f$\ms{CHAST}\f$    |   Chatham Standard Time
    !>  \f$+60 \times 13     \f$    |   \f$\ms{TOT  }\f$    |   Tonga Time
    !>  \f$+60 \times 13 + 45\f$    |   \f$\ms{CHADT}\f$    |   Chatham Daylight Time
    !>  \f$+60 \times 14     \f$    |   \f$\ms{LINT }\f$    |   Line Islands Time
    !>
    !>  \see
    !>  [timeZone](@ref pm_dateTime::timeZone)<br>
    !>  [getZoneAbbr](@ref pm_dateTime::getZoneAbbr)<br>
    !>
    !>  \final{timeZone_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, January 30, 2021, 5:36 AM, Dallas, TX
    type :: timeZone_type
        integer(IK)     :: Zone(39) =   [ -60 * 12          &
                                        , -60 * 11          &
                                        , -60 * 10          &
                                        , -60 * 9 - 30      &
                                        , -60 * 9           &
                                        , -60 * 8           &
                                        , -60 * 7           &
                                        , -60 * 6           &
                                        , -60 * 5           &
                                        , -60 * 3 - 30      &
                                        , -60 * 3           &
                                        , -60 * 2 - 30      &
                                        , -60 * 2           &
                                        , -60 * 1           &
                                        , +60 * 0           &
                                        , +60 * 1           &
                                        , +60 * 2           &
                                        , +60 * 3           &
                                        , +60 * 3 + 30      &
                                        , +60 * 4           &
                                        , +60 * 4 + 30      &
                                        , +60 * 5           &
                                        , +60 * 5 + 30      &
                                        , +60 * 5 + 45      &
                                        , +60 * 6           &
                                        , +60 * 6 + 30      &
                                        , +60 * 7           &
                                        , +60 * 8           &
                                        , +60 * 8 + 45      &
                                        , +60 * 9           &
                                        , +60 * 9 + 30      &
                                        , +60 * 10          &
                                        , +60 * 10 + 30     &
                                        , +60 * 11          &
                                        , +60 * 12          &
                                        , +60 * 12 + 45     &
                                        , +60 * 13          &
                                        , +60 * 13 + 45     &
                                        , +60 * 14          &
                                        ]
        character(5,SK) :: Abbr(39) =   [ "IDLW " &
                                        , "SST  " &
                                        , "HST  " &
                                        , "MIT  " &
                                        , "AKST " &
                                        , "PST  " &
                                        , "MST  " &
                                        , "CST  " &
                                        , "EST  " &
                                        , "NST  " &
                                        , "UYT  " &
                                        , "NDT  " &
                                        , "UYST " &
                                        , "EGT  " &
                                        , "UTC  " &
                                        , "CET  " &
                                        , "EET  " &
                                        , "AST  " &
                                        , "IRST " &
                                        , "GET  " &
                                        , "AFT  " &
                                        , "PKT  " &
                                        , "IST  " &
                                        , "NPT  " &
                                        , "BST  " &
                                        , "MMT  " &
                                        , "THA  " &
                                        , "SST  " &
                                        , "CWST " &
                                        , "JST  " &
                                        , "ACST " &
                                        , "AEST " &
                                        , "LHST " &
                                        , "PONT " &
                                        , "NZST " &
                                        , "CHAST" &
                                        , "TOT  " &
                                        , "CHADT" &
                                        , "LINT " &
                                        ]
    end type

    !>  \brief
    !>  This is an object parameter of type [timeZone_type](@ref pm_dateTime::timeZone_type) containing a list of time zones and their representative abbreviations.<br>
    !>
    !>  \details
    !>  This list contains only a representative set of zone abbreviations.<br>
    !>  See the documentation of [timeZone_type](@ref pm_dateTime::timeZone_type) for more information.<br>
    !>  See the documentation of [getZoneAbbr](@ref pm_dateTime::getZoneAbbr) for relevant usage and examples.<br>
    !>  The primary use case of this object is in the implementation of [getZoneAbbr()](@ref pm_dateTime::getZoneAbbr).<br>
    !>
    !>  \see
    !>  [getZoneAbbr](@ref pm_dateTime::getZoneAbbr)<br>
    !>  [timeZone_type](@ref pm_dateTime::timeZone_type)<br>
    !>
    !>  \final{timeZone}
    !>
    !>  \author
    !>  \AmirShahmoradi, January 30, 2021, 5:36 AM, Dallas, TX
    type(timeZone_type), parameter :: timeZone = timeZone_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: timeZone
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the time zone abbreviation corresponding to the current local time or the specified `zone` in **minutes**.
    !>
    !>  \details
    !>  The time zone abbreviation is difficult to infer and generally requires communication with the operating system, which may or may not have it.<br>
    !>  The current implementation of this generic interface relies on a predefined internal list of time zone abbreviations to convert the specified `zone` **in units of minutes** to an abbreviation.<br>
    !>  Where there are multiple time zone abbreviations available for a single time zone, only a single representative abbreviation is kept in the list.<br>
    !>  See the documentation of [timeZone_type](@ref pm_dateTime::timeZone_type) for the current internal list of zones and abbreviations used by the ParaMonte library.<br>
    !>
    !>  \param[in]  zone    :   The input scalar of type `integer` of default kind \IK, containing the <b>local time zone</b> of the Gregorian calendar <b>in minutes</b>.<br>
    !>                          (**optional**, default = [getZone()](@ref pm_dateTime::getZone))
    !>  \return
    !>  `abbr`              :   The output scalar of type `integer` of default kind \IK containing the local
    !>                          time difference **in minutes** with respect to the Coordinated Universal Time (UTC).
    !>
    !>  \interface{getZoneAbbr}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getZoneAbbr
    !>      character(:, SK), allocatable :: abbr
    !>      integer(IK) :: zone
    !>
    !>      abbr = getZoneAbbr()
    !>      abbr = getZoneAbbr(zone)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  When the input argument `zone` is out of the supported range (\f$-12\times60\f$ minutes UTC and \f$+14\times60\f$ minutes UTC), the output string will be empty.<br>
    !>
    !>  \warning
    !>  This generic interface does not currently take into account the daylight savings calendar.<br>
    !>  This can lead to shifts in the output zones by up to one hour.<br>
    !>
    !>  \see
    !>  [getMillisecond](@ref pm_dateTime::getMillisecond)<br>
    !>  [getSecond](@ref pm_dateTime::getSecond)<br>
    !>  [getMinute](@ref pm_dateTime::getMinute)<br>
    !>  [getHour12](@ref pm_dateTime::getHour12)<br>
    !>  [getHour](@ref pm_dateTime::getHour)<br>
    !>  [getZone](@ref pm_dateTime::getZone)<br>
    !>  [getZoneAbbr](@ref pm_dateTime::getZoneAbbr)<br>
    !>  [getDay](@ref pm_dateTime::getDay)<br>
    !>  [getMonth](@ref pm_dateTime::getMonth)<br>
    !>  [getYear](@ref pm_dateTime::getYear)<br>
    !>  [isLeapYear](@ref pm_dateTime::isLeapYear)<br>
    !>
    !>  \example{getZoneAbbr}
    !>  \include{lineno} example/pm_dateTime/getZoneAbbr/main.F90
    !>  \compilef{getZoneAbbr}
    !>  \output{getZoneAbbr}
    !>  \include{lineno} example/pm_dateTime/getZoneAbbr/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \todo
    !>  \phigh
    !>  Currently, the zone abbreviation is derived from a predefined list.<br>
    !>  The zone abbreviation must be inferred directly from the operating system.<br>
    !>  On Windows, this could be also be done via the Powershell command `(Get-Date).IsDaylightSavingTime()` to test whether the daylight savings is activated.<br>
    !>  On Unix, The command `date +"\%Z"` outputs the zone abbreviation.<br>
    !>
    !>  \final{getZoneAbbr}
    !>
    !>  \author
    !>  \AmirShahmoradi, January 30, 2021, 5:36 AM, Dallas, TX
    interface getZoneAbbr

    module function getZoneAbbrC() result(abbr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getZoneAbbrC
#endif
        use pm_kind, only: SKG => SK
        character(:,SKG), allocatable   :: abbr
    end function

    PURE module function getZoneAbbrZ(zone) result(abbr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getZoneAbbrZ
#endif
        use pm_kind, only: IKG => IK, SKG => SK
        integer(IKG)    , intent(in)    :: zone
        character(:,SKG), allocatable   :: abbr
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the current **12-hour-clock** local hour of the current day of the Gregorian calendar.
    !>
    !>  \param[in]  hour    :   The input scalar or array of arbitrary shape of type `integer` of default kind \IK
    !>                          containing the **24-hour-clock** to be converted to a **12-hour-clock** hour.<br>
    !>                          (**optional**, default = [getHour()](@ref pm_dateTime::getHour))
    !>
    !>  \return
    !>  `hour`              :   The output scalar or array of the same rank, shape and size as the input `hour`,
    !>                          of the same type and kind as `hour`, containing the input `hour` converted to **12-hour-clock** hour.
    !>
    !>  \interface{getHour12}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getHour12
    !>      integer(IK) :: hour, hour12
    !>
    !>      hour12 = getHour12()
    !>      hour12 = getHour12(hour)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>  The procedures under this generic interface are always `impure` when the input argument `hour` is missing.<br>
    !>
    !>  \remark
    !>  The procedures under this generic interface are `elemental` when the input argument `hour` is present.<br>
    !>
    !>  \see
    !>  [getMillisecond](@ref pm_dateTime::getMillisecond)<br>
    !>  [getSecond](@ref pm_dateTime::getSecond)<br>
    !>  [getMinute](@ref pm_dateTime::getMinute)<br>
    !>  [getHour12](@ref pm_dateTime::getHour12)<br>
    !>  [getHour](@ref pm_dateTime::getHour)<br>
    !>  [getZone](@ref pm_dateTime::getZone)<br>
    !>  [getZoneAbbr](@ref pm_dateTime::getZoneAbbr)<br>
    !>  [getDay](@ref pm_dateTime::getDay)<br>
    !>  [getMonth](@ref pm_dateTime::getMonth)<br>
    !>  [getYear](@ref pm_dateTime::getYear)<br>
    !>  [isLeapYear](@ref pm_dateTime::isLeapYear)<br>
    !>
    !>  \example{getHour12}
    !>  \include{lineno} example/pm_dateTime/getHour12/main.F90
    !>  \compilef{getHour12}
    !>  \output{getHour12}
    !>  \include{lineno} example/pm_dateTime/getHour12/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getHour12}
    !>
    !>  \author
    !>  \AmirShahmoradi, January 30, 2021, 5:36 AM, Dallas, TX
    interface getHour12

    module function getHour12C() result(hour12)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHour12C
#endif
        use pm_kind, only: IKG => IK
        integer(IKG) :: hour12
    end function

    PURE elemental module function getHour12H(hour) result(hour12)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHour12H
#endif
        use pm_kind, only: IKG => IK
        integer(IKG), intent(in) :: hour
        integer(IKG)             :: hour12
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Julian Date (Julian Day Number JDN + the fractional part of the day) from
    !>  the input `[year, month, day, zone, hour, minute, second, millisecond]` of the Gregorian calendar date.
    !>
    !>  \details
    !>  The algrithm of this generic interface is valid for any Gregorian date, even proleptic Gregorian dates including those for negative years.
    !>  The Julian Day Number, Julian Day, or JD of a particular instant of time is the number of days and fractions of a day since `12` hours Universal Time (Greenwich mean noon)
    !>  on January 1 of the year `-4712`, where the year is given in the Julian proleptic calendar. The idea of using this reference date was originally proposed by Joseph Scalizer
    !>  in 1582 to count years but it was modified by 19th century astronomers to count days.<br>
    !>  **Julian days are Julian Day Numbers and are not to be confused with Julian dates.**<br>
    !>  The Julian Day Number (JDN) is expressed as an integer and it represents the number of whole days since the reference instant to **noon of that day**. For example,
    !>  <br>
    !>  Gregorian Date      | Time of Day                                   | JD
    !>  --------------------|-----------------------------------------------|----------
    !>  November 24, -4713  | start of the day (just after midnight)        | -0.5
    !>  November 24, -4713  | noon (start of JD in Gregorian Calendar)      | 0.0
    !>  November 25, -4713  | start of the day (just after midnight)        | +0.5
    !>  January 1, -1       | start of day                                  | 1720694.5
    !>  October 15, 1582    | start of day (first day of Gregorian reform)  | 2299160.5
    !>  January 1, 1901     | start of day (start of the 20th century)      | 2415385.5
    !>  January 1, 1970     | start of day (Unix reference date)            | 2440587.5
    !>  December 31, 1979   | noon                                          | 2444239.0
    !>  January 1, 1980     | start of the day (just after midnight)        | 2444239.5
    !>  January 1, 1980     | noon (Microsoft DOS reference date)           | 2444240.0
    !>  January 1, 1980     | midnight commencing January 2                 | 2444240.5
    !>  <br>
    !>  A **Julian date** is a date in the Julian calendar, similar to a **Gregorian date** in the Gregorian calendar.<br>
    !>  Note that the Julian Day corresponding to a Julian Date does not have the same value as the JD corresponding to the same date in the Gregorian Calendar.
    !>
    !>  \param[in]  year        :   The input scalar of type `integer` of default kind \IK, containing the year of the Gregorian calendar.<br>
    !>                              (**optional**, default =  the current value if all input arguments are missing. It can be present only if `values` is missing.)
    !>  \param[in]  month       :   The input scalar of type `integer` of default kind \IK, containing the month of the year of the Gregorian calendar.<br>
    !>                              (**optional**, default = `1` (or the current value if all input arguments are missing). It can be present only if `year` is present.)
    !>  \param[in]  day         :   The input scalar of type `integer` of default kind \IK, containing the day of the month of the year of the Gregorian calendar.<br>
    !>                              (**optional**, default = `1`. (or the current value if all input arguments are missing). It can be present only if `month` is present.)
    !>  \param[in]  zone        :   The input scalar of type `integer` of default kind \IK, containing the <b>time zone</b> of the Gregorian calendar <b>in minutes</b>.<br>
    !>                              (**optional**, default = `0` (UTC) (or the current value if all input arguments are missing). It can be present only if `day` is present.)
    !>  \param[in]  hour        :   The input scalar of type `integer` of default kind \IK, containing the hour of the day of the Gregorian calendar.<br>
    !>                              (**optional**, default = `0` (or the current value if all input arguments are missing). It can be present only if `zone` is present.)
    !>  \param[in]  minute      :   The input scalar of type `integer` of default kind \IK, containing the minute of the hour of the day of the Gregorian calendar.<br>
    !>                              (**optional**, default = `0` (or the current value if all input arguments are missing). It can be present only if `hour` is present.)
    !>  \param[in]  second      :   The input scalar of type `integer` of default kind \IK, containing the second of the minute of the hour of the day of the Gregorian calendar.<br>
    !>                              (**optional**, default = `0` (or the current value if all input arguments are missing). It can be present only if `minute` is present.)
    !>  \param[in]  millisecond :   The input scalar of type `integer` of default kind \IK, containing the millisecond of the second of the minute of the hour of the day of the Gregorian calendar.<br>
    !>                              (**optional**, default = `0` (or the current value if all input arguments are missing). It can be present only if `second` is present.)
    !>  \param[in]  values      :   The input `contiguous` vector of maximum size `8` of type `integer` of default kind \IK, containing the values
    !>                              `[year, month, day, zone, hour, minute, seconds, milliseconds]` of the Gregorian calendar or a subset of the octuple starting with `year`.<br>
    !>                              The order of the elements of the vector follows that of the `values` returned by the Fortran intrinsic `date_and_time()`.<br>
    !>                              (**optional**, default = the current date and time. It can be present **only if** all other arguments are missing.)
    !>
    !>  \return
    !>  `julianDay`            :   The output scalar or array of the same shape as the input array-like arguments (except `values`), of type `real` of default kind \RK,
    !>                              representing the <b>Julian Date equivalent (in units of days, possibly fractional)</b> of the specified Gregorian Calendar date.
    !>
    !>  \interface{getJulianDay}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getJulianDay
    !>      real(RK) :: julianDay
    !>
    !>      julianDay = getJulianDay() ! Current Julian Date (JD)
    !>      julianDay = getJulianDay(year)
    !>      julianDay = getJulianDay(year, month)
    !>      julianDay = getJulianDay(year, month, day)
    !>      julianDay = getJulianDay(year, month, day, zone)
    !>      julianDay = getJulianDay(year, month, day, zone, hour)
    !>      julianDay = getJulianDay(year, month, day, zone, hour, minute)
    !>      julianDay = getJulianDay(year, month, day, zone, hour, minute, second)
    !>      julianDay = getJulianDay(year, month, day, zone, hour, minute, second, millisecond)
    !>      julianDay = getJulianDay(values(:)) ! values = [year, month, day, zone, hour, minute, second, millisecond] or a subset starting with `year`.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < size(values) < 9` must hold.<br>
    !>  The input values for `[year, month, day, zone, hour, minute, second, millisecond]` must be valid and consistent with each other.<br>
    !>  For example, the input `month` must be a number between `1` and `12` and `day` must be between `1` and `31`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  The procedures under this generic interface are always non-elemental and `impure` when no input argument is present.
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  The procedures under this generic interface are non-elemental when the input argument `values(:)` is present.
    !>
    !>  \note
    !>  This generic interface is particularly useful for efficient computation of the number of days between two Gregorian dates.
    !>
    !>  \see
    !>  [getDateTimeDiff](@ref pm_dateTime::getDateTimeDiff)<br>
    !>  [JPL Gregorian to Julian Day Number Converter](https://ssd.jpl.nasa.gov/tools/jdc/#/cd)<br>
    !>  [A One-Line Algorithm for Julian Date](https://ui.adsabs.harvard.edu/abs/1983IAPPP..13...16F/abstract)<br>
    !>  [Hatcher, 1984, Simple Formulae for Julian Day Numbers and Calendar Dates](https://ui.adsabs.harvard.edu/abs/1984QJRAS..25...53H/abstract)<br>
    !>  [Baum, 2017, Date Algorithms](https://www.researchgate.net/publication/316558298_Date_Algorithms)<br>
    !>
    !>  \example{getJulianDay}
    !>  \include{lineno} example/pm_dateTime/getJulianDay/main.F90
    !>  \compilef{getJulianDay}
    !>  \output{getJulianDay}
    !>  \include{lineno} example/pm_dateTime/getJulianDay/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getJulianDay}
    !>  If you use this generic interface, you must also cite the paper by [Peter Baum, 2017, Date Algorithms](https://www.researchgate.net/publication/316558298_Date_Algorithms).<br>
    !>
    !>  \author
    !>  \AmirShahmoradi, May 9, 2022, 5:20 AM, Albuquerque, New Mexico
    interface getJulianDay

    impure module function getJulianDayC() result(julianDay)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getJulianDayC
#endif
        real(RK)                                    :: julianDay
    end function

    PURE module function getJulianDayV(values) result(julianDay)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getJulianDayV
#endif
        integer(IK)     , intent(in), contiguous    :: values(:)
        real(RK)                                    :: julianDay
    end function

    PURE elemental module function getJulianDayY(year) result(julianDay)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getJulianDayY
#endif
        integer(IK)     , intent(in)                :: year
        real(RK)                                    :: julianDay
    end function

    PURE elemental module function getJulianDayYM(year, month) result(julianDay)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getJulianDayYM
#endif
        integer(IK)     , intent(in)                :: year, month
        real(RK)                                    :: julianDay
    end function

    PURE elemental module function getJulianDayYMD(year, month, day) result(julianDay)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getJulianDayYMD
#endif
        integer(IK)     , intent(in)                :: year, month, day
        real(RK)                                    :: julianDay
    end function

    PURE elemental module function getJulianDayYMDZ(year, month, day, zone) result(julianDay)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getJulianDayYMDZ
#endif
        integer(IK)     , intent(in)                :: year, month, day, zone
        real(RK)                                    :: julianDay
    end function

    PURE elemental module function getJulianDayYMDZH(year, month, day, zone, hour) result(julianDay)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getJulianDayYMDZH
#endif
        integer(IK)     , intent(in)                :: year, month, day, zone, hour
        real(RK)                                    :: julianDay
    end function

    PURE elemental module function getJulianDayYMDZHM(year, month, day, zone, hour, minute) result(julianDay)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getJulianDayYMDZHM
#endif
        integer(IK)     , intent(in)                :: year, month, day, zone, hour, minute
        real(RK)                                    :: julianDay
    end function

    PURE elemental module function getJulianDayYMDZHMS(year, month, day, zone, hour, minute, second) result(julianDay)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getJulianDayYMDZHMS
#endif
        integer(IK)     , intent(in)                :: year, month, day, zone, hour, minute, second
        real(RK)                                    :: julianDay
    end function

    PURE elemental module function getJulianDayYMDZHMSM(year, month, day, zone, hour, minute, second, millisecond) result(julianDay)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getJulianDayYMDZHMSM
#endif
        integer(IK)     , intent(in)                :: year, month, day, zone, hour, minute, second, millisecond
        real(RK)                                    :: julianDay
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the current or the requested Gregorian date and time shifted by the specified `amount` in units of days.
    !>
    !>  \details
    !>  This algorithm carefully takes into account the possibility of leap years.<br>
    !>  If the input `amount` is in units other than days, such as hours, minutes or seconds, simply use the relevant module constants
    !>  (e.g., [MINUTES_PER_DAY](@ref pm_dateTime::MINUTES_PER_DAY) or [SECONDS_PER_DAY](@ref pm_dateTime::SECONDS_PER_DAY))
    !>  to convert the non-day measures of time to days.
    !>
    !>  \param[in]  amount      :   The input scalar of type `real` of default kind \RK, containing the <b>number of days</b> (possibly fractional) with which the input date and time must be shifted.<br>
    !>                              It can practically be any positive or negative value. Note that a day is `86400` seconds ([SECONDS_PER_DAY](@ref pm_dateTime::SECONDS_PER_DAY)).<br>
    !>  \param[in]  year        :   The input scalar of type `integer` of default kind \IK, containing the year of the Gregorian calendar.<br>
    !>                              (**optional**, default =  the current value if all input arguments except `amount` are missing. It can be present only if `values` is missing.)
    !>  \param[in]  month       :   The input scalar of type `integer` of default kind \IK, containing the month of the year of the Gregorian calendar.<br>
    !>                              (**optional**, default = `1` (or the current value if all input arguments are missing). It can be present only if `year` is also present.)
    !>  \param[in]  day         :   The input scalar of type `integer` of default kind \IK, containing the day of the month of the year of the Gregorian calendar.<br>
    !>                              (**optional**, default = `1`. (or the current value if all input arguments are missing). It can be present only if `month` is also present.)
    !>  \param[in]  zone        :   The input scalar of type `integer` of default kind \IK, containing the <b>local time zone</b> of the Gregorian calendar <b>in minutes</b>.<br>
    !>                              The input argument `zone` has no effects on the output shifted date and time. It is only used to construct the output date and time vector.<br>
    !>                              (**optional**, default = `0` (UTC) (or the current value if all input arguments are missing). It can be present only if `day` is also present.)
    !>  \param[in]  hour        :   The input scalar of type `integer` of default kind \IK, containing the hour of the day of the Gregorian calendar.<br>
    !>                              (**optional**, default = `0` (or the current value if all input arguments are missing). It can be present only if `zone` is present.)
    !>  \param[in]  minute      :   The input scalar of type `integer` of default kind \IK, containing the minute of the hour of the day of the Gregorian calendar.<br>
    !>                              (**optional**, default = `0` (or the current value if all input arguments are missing). It can be present only if `hour` is present.)
    !>  \param[in]  second      :   The input scalar of type `integer` of default kind \IK, containing the second of the minute of the hour of the day of the Gregorian calendar.<br>
    !>                              (**optional**, default = `0` (or the current value if all input arguments are missing). It can be present only if `minute` is present.)
    !>  \param[in]  millisecond :   The input scalar of type `integer` of default kind \IK, containing the millisecond of the second of the minute of the hour of the day of the Gregorian calendar.<br>
    !>                              (**optional**, default = `0` (or the current value if all input arguments are missing). It can be present only if `second` is present.)
    !>  \param[in]  values      :   The input `contiguous` vector of minimum size `1` and maximum size `8` of type `integer` of default kind \IK, containing the values
    !>                              `[year, month, day, zone, hour, minute, seconds, milliseconds]` of the Gregorian calendar or a subset of the octuple starting with `year`.<br>
    !>                              The order of the elements of the vector follows that of the `values` returned by the Fortran intrinsic `date_and_time()`.<br>
    !>                              (**optional**, default = the current date and time. It can be present **only if** all other arguments except `amount` are missing.)
    !>
    !>  \return
    !>  `dateTimeShifted(1:8)`  :   The output vector of size `8` of type `integer` of default kind \IK containing
    !>                              the requested local date and time of the Gregorian Calendar corresponding to the input
    !>                              date and time shifted by the specified `amount` in the order `[year, month, day, zone, hour, minute, seconds, milliseconds]`.
    !>                              By definition, `dateTimeShifted(4)` (corresponding to the local time zone) is the same as `values(4)` or the `zone` input argument.
    !>
    !>  \interface{getDateTimeShifted}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: SK, IK
    !>      use pm_dateTime, only: getDateTimeShifted
    !>      integer(IK) :: values(8)
    !>      integer(IK) :: dateTimeShifted(8)
    !>
    !>      dateTimeShifted(1:8) = getDateTimeShifted(amount) ! return the current local date and time shifted by the specified `amount`.
    !>      dateTimeShifted(1:8) = getDateTimeShifted(amount, year, month, day, zone) ! return the specified date shifted by the specified `amount`.
    !>      dateTimeShifted(1:8) = getDateTimeShifted(amount, year, month, day, zone, hour)
    !>      dateTimeShifted(1:8) = getDateTimeShifted(amount, year, month, day, zone, hour, minute)
    !>      dateTimeShifted(1:8) = getDateTimeShifted(amount, year, month, day, zone, hour, minute, second)
    !>      dateTimeShifted(1:8) = getDateTimeShifted(amount, year, month, day, zone, hour, minute, second, millisecond)
    !>      dateTimeShifted(1:8) = getDateTimeShifted(amount, values(1:8)) ! values(1:8) = [year, month, day, zone, hour, minute, second, millisecond]
    !>      dateTimeShifted(1:8) = getDateTimeShifted(amount, values(1:7)) ! values(1:7) = [year, month, day, zone, hour, minute, second]
    !>      dateTimeShifted(1:8) = getDateTimeShifted(amount, values(1:6)) ! values(1:6) = [year, month, day, zone, hour, minute]
    !>      dateTimeShifted(1:8) = getDateTimeShifted(amount, values(1:5)) ! values(1:5) = [year, month, day, zone, hour]
    !>      dateTimeShifted(1:8) = getDateTimeShifted(amount, values(1:4)) ! values(1:4) = [year, month, day, zone]
    !>      dateTimeShifted(1:8) = getDateTimeShifted(amount, values(1:3)) ! values(1:4) = [year, month, day]
    !>      dateTimeShifted(1:8) = getDateTimeShifted(amount, values(1:2)) ! values(1:4) = [year, month]
    !>      dateTimeShifted(1:8) = getDateTimeShifted(amount, values(1:1)) ! values(1:1) = [year]
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The input `month` must be a number between `1` and `12`.<br>
    !>  The input `day` must be a number between `1` and `31` and be consistent with the specified month and (leap) year.<br>
    !>  The input `zone` must be a valid time zone in units of minutes.<br>
    !>  The input `hour` must be a number between `0` and `23`.<br>
    !>  The input `minute` must be a number between `0` and `59`.<br>
    !>  The input `second` must be a number between `0` and `59`.<br>
    !>  The input `millisecond` must be a number between `0` and `999`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  The procedures under this generic interface are always `impure` when all input arguments are missing.
    !>
    !>  \see
    !>  [getJulianDay](@ref pm_dateTime::getJulianDay)<br>
    !>  [getDateTimeDiff](@ref pm_dateTime::getDateTimeDiff)<br>
    !>
    !>  \example{getDateTimeShifted}
    !>  \include{lineno} example/pm_dateTime/getDateTimeShifted/main.F90
    !>  \compilef{getDateTimeShifted}
    !>  \output{getDateTimeShifted}
    !>  \include{lineno} example/pm_dateTime/getDateTimeShifted/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getDateTimeShifted}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 3:59 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface getDateTimeShifted

    impure module function getDateTimeShiftedC(amount) result(dateTimeShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeShiftedC
#endif
        real(RK)        , intent(in)                :: amount
        integer(IK)                                 :: dateTimeShifted(8)
    end function

    PURE module function getDateTimeShiftedV(amount, values) result(dateTimeShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeShiftedV
#endif
        real(RK)        , intent(in)                :: amount
        integer(IK)     , intent(in), contiguous    :: values(:)
        integer(IK)                                 :: dateTimeShifted(8)
    end function

    PURE module function getDateTimeShiftedY(amount, year) result(dateTimeShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeShiftedY
#endif
        real(RK)        , intent(in)                :: amount
        integer(IK)     , intent(in)                :: year
        integer(IK)                                 :: dateTimeShifted(8)
    end function

    PURE module function getDateTimeShiftedYM(amount, year, month) result(dateTimeShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeShiftedYM
#endif
        real(RK)        , intent(in)                :: amount
        integer(IK)     , intent(in)                :: year, month
        integer(IK)                                 :: dateTimeShifted(8)
    end function

    PURE module function getDateTimeShiftedYMD(amount, year, month, day) result(dateTimeShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeShiftedYMD
#endif
        real(RK)        , intent(in)                :: amount
        integer(IK)     , intent(in)                :: year, month, day
        integer(IK)                                 :: dateTimeShifted(8)
    end function

    PURE module function getDateTimeShiftedYMDZ(amount, year, month, day, zone) result(dateTimeShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeShiftedYMDZ
#endif
        real(RK)        , intent(in)                :: amount
        integer(IK)     , intent(in)                :: year, month, day, zone
        integer(IK)                                 :: dateTimeShifted(8)
    end function

    PURE module function getDateTimeShiftedYMDZH(amount, year, month, day, zone, hour) result(dateTimeShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeShiftedYMDZH
#endif
        real(RK)        , intent(in)                :: amount
        integer(IK)     , intent(in)                :: year, month, day, zone, hour
        integer(IK)                                 :: dateTimeShifted(8)
    end function

    PURE module function getDateTimeShiftedYMDZHM(amount, year, month, day, zone, hour, minute) result(dateTimeShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeShiftedYMDZHM
#endif
        real(RK)        , intent(in)                :: amount
        integer(IK)     , intent(in)                :: year, month, day, zone, hour, minute
        integer(IK)                                 :: dateTimeShifted(8)
    end function

    PURE module function getDateTimeShiftedYMDZHMS(amount, year, month, day, zone, hour, minute, second) result(dateTimeShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeShiftedYMDZHMS
#endif
        real(RK)        , intent(in)                :: amount
        integer(IK)     , intent(in)                :: year, month, day, zone, hour, minute, second
        integer(IK)                                 :: dateTimeShifted(8)
    end function

    PURE module function getDateTimeShiftedYMDZHMSM(amount, year, month, day, zone, hour, minute, second, millisecond) result(dateTimeShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeShiftedYMDZHMSM
#endif
        real(RK)        , intent(in)                :: amount
        integer(IK)     , intent(in)                :: year, month, day, zone, hour, minute, second, millisecond
        integer(IK)                                 :: dateTimeShifted(8)
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the calendarical difference in days between the two input Gregorian dates octuples (Values1(:) - Values2(:)).
    !>
    !>  \details
    !>  This algorithm carefully takes into account the possibility of leap years and different time zones.<br>
    !>  If `Values1` is smaller than `Values2` the difference between the two dates is negative.
    !>
    !>  \param[in]  Values1 :   The input `contiguous` vector of size `8` of type `integer` of default kind \IK, containing the octuple
    !>                          `[year, month, day, zone, hour, minute, seconds, milliseconds]` of the Gregorian calendar from which `Values2` should be subtracted.<br>
    !>                          The order of the elements of the vector follows that of the `values` returned by the Fortran intrinsic `date_and_time()`.<br>
    !>  \param[in]  Values2 :   The input `contiguous` vector of size `8` of type `integer` of default kind \IK, containing the octuple
    !>                          `[year, month, day, zone, hour, minute, seconds, milliseconds]` of the Gregorian calendar that should be subtracted from `Values1`.<br>
    !>                          The order of the elements of the vector follows that of the `values` returned by the Fortran intrinsic `date_and_time()`.<br>
    !>
    !>  \return
    !>  `dateTimeDiff`      :   The output scalar of type `real` of default kind \RK containing the result (<b>in units of days</b>, possible fractional)
    !>                          of the subtraction of `Values2` date octuple from `Values1` date octuple (`Values1 - Values2`).
    !>
    !>  \interface{getDateTimeDiff}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getDateTimeDiff
    !>      use pm_kind, only: IK, RK
    !>      integer(IK) :: Values1(8)
    !>      integer(IK) :: Values2(8)
    !>      real(RK) :: dateTimeDiff
    !>
    !>      dateTimeDiff = getDateTimeDiff(Values1(1:8), Values2(1:8)) ! return Values1 - Values2 in units of days.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The sizes of the two input vector arguments must be the same and equal to `8`.<br>
    !>  The input `month` must be a number between `1` and `12`.<br>
    !>  The input `day` must be a number between `1` and `31` and be consistent with the specified month and (leap) year.<br>
    !>  The input `zone` must be a valid time zone in units of minutes.<br>
    !>  The input `hour` must be a number between `0` and `23`.<br>
    !>  The input `minute` must be a number between `0` and `59`.<br>
    !>  The input `second` must be a number between `0` and `59`.<br>
    !>  The input `millisecond` must be a number between `0` and `999`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  The procedures under this generic interface are always `impure` when all input arguments are missing.<br>
    !>
    !>  \see
    !>  [getJulianDay](@ref pm_dateTime::getJulianDay)<br>
    !>  [getDateTimeShifted](@ref pm_dateTime::getDateTimeShifted)<br>
    !>
    !>  \example{getDateTimeDiff}
    !>  \include{lineno} example/pm_dateTime/getDateTimeDiff/main.F90
    !>  \compilef{getDateTimeDiff}
    !>  \output{getDateTimeDiff}
    !>  \include{lineno} example/pm_dateTime/getDateTimeDiff/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getDateTimeDiff}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 3:59 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface getDateTimeDiff

    PURE module function getDateTimeDiffValues(Values1, Values2) result(dateTimeDiff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeDiffValues
#endif
        integer(IK)     , intent(in), contiguous    :: Values1(:), Values2(:)
        real(RK)                                    :: dateTimeDiff
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the current or the requested date and time converted to the
    !>  corresponding Coordinated Universal Time (UTC) as an integer-valued array of size `8`.
    !>
    !>  \details
    !>  The first four input arguments or the first four elements of the input vector `values(:)`
    !>  are mandatory and must be provided since they determine the date and local time zone of interest.<br>
    !>  The rest of the elements containing the time of the day of the Gregorian Calendar are optional and taken to be zero if not provided.<br>
    !>
    !>  \param[in]  year        :   The input scalar of type `integer` of default kind \IK, containing the year of the Gregorian calendar.<br>
    !>                              (**optional**, default =  the current value if all input arguments are missing. It can be present if only if `month`, `day`, and `zone` are also present and `values` is missing.)
    !>  \param[in]  month       :   The input scalar of type `integer` of default kind \IK, containing the month of the year of the Gregorian calendar.<br>
    !>                              (**optional**, default = `1` (or the current value if all input arguments are missing). It **must** be present if `year` is also present.)
    !>  \param[in]  day         :   The input scalar of type `integer` of default kind \IK, containing the day of the month of the year of the Gregorian calendar.<br>
    !>                              (**optional**, default = `1`. (or the current value if all input arguments are missing). It **must**  be present if `month` is also present.)
    !>  \param[in]  zone        :   The input scalar of type `integer` of default kind \IK, containing the <b>local time zone</b> of the Gregorian calendar <b>in minutes</b>.<br>
    !>                              (**optional**, default = `0` (UTC) (or the current value if all input arguments are missing). It **must** be present if `day` is also present.)
    !>  \param[in]  hour        :   The input scalar of type `integer` of default kind \IK, containing the hour of the day of the Gregorian calendar.<br>
    !>                              (**optional**, default = `0` (or the current value if all input arguments are missing). It can be present only if `zone` is present.)
    !>  \param[in]  minute      :   The input scalar of type `integer` of default kind \IK, containing the minute of the hour of the day of the Gregorian calendar.<br>
    !>                              (**optional**, default = `0` (or the current value if all input arguments are missing). It can be present only if `hour` is present.)
    !>  \param[in]  second      :   The input scalar of type `integer` of default kind \IK, containing the second of the minute of the hour of the day of the Gregorian calendar.<br>
    !>                              (**optional**, default = `0` (or the current value if all input arguments are missing). It can be present only if `minute` is present.)
    !>  \param[in]  millisecond :   The input scalar of type `integer` of default kind \IK, containing the millisecond of the second of the minute of the hour of the day of the Gregorian calendar.<br>
    !>                              (**optional**, default = `0` (or the current value if all input arguments are missing). It can be present only if `second` is present.)
    !>  \param[in]  values      :   The input `contiguous` vector of minimum size `4` and maximum size `8` of type `integer` of default kind \IK, containing the values
    !>                              `[year, month, day, zone, hour, minute, seconds, milliseconds]` of the Gregorian calendar or a subset of the octuple starting with `year`.<br>
    !>                              The order of the elements of the vector follows that of the `values` returned by the Fortran intrinsic `date_and_time()`.<br>
    !>                              (**optional**, default = the current date and time. It can be present **only if** all other arguments are missing.)
    !>
    !>  \return
    !>  `DateTimeUTC(1:8)`      :   The output vector of size `8` of type `integer` of default kind \IK containing
    !>                              the requested UTC date and time of the Gregorian Calendar corresponding to the input local date and time
    !>                              in the order `[year, month, day, zone, hour, minute, seconds, milliseconds]`.
    !>                              By definition, `DateTimeUTC(4)` (corresponding to the UTC time zone) is always zero.
    !>
    !>  \interface{getDateTimeUTC}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: SK, IK
    !>      use pm_dateTime, only: getDateTimeUTC
    !>      integer(IK) :: DateTimeUTC(8)
    !>      integer(IK) :: DateTime(8)
    !>
    !>      DateTimeUTC(1:8) = getDateTimeUTC() ! return the current UTC date and time.
    !>      DateTimeUTC(1:8) = getDateTimeUTC(year, month, day, zone)
    !>      DateTimeUTC(1:8) = getDateTimeUTC(year, month, day, zone, hour)
    !>      DateTimeUTC(1:8) = getDateTimeUTC(year, month, day, zone, hour, minute)
    !>      DateTimeUTC(1:8) = getDateTimeUTC(year, month, day, zone, hour, minute, second)
    !>      DateTimeUTC(1:8) = getDateTimeUTC(year, month, day, zone, hour, minute, second, millisecond)
    !>      DateTimeUTC(1:8) = getDateTimeUTC(DateTime(1:8)) ! DateTime(1:8) = [year, month, day, zone, hour, minute, second, millisecond]
    !>      DateTimeUTC(1:8) = getDateTimeUTC(DateTime(1:7)) ! DateTime(1:7) = [year, month, day, zone, hour, minute, second]
    !>      DateTimeUTC(1:8) = getDateTimeUTC(DateTime(1:6)) ! DateTime(1:6) = [year, month, day, zone, hour, minute]
    !>      DateTimeUTC(1:8) = getDateTimeUTC(DateTime(1:5)) ! DateTime(1:5) = [year, month, day, zone, hour]
    !>      DateTimeUTC(1:8) = getDateTimeUTC(DateTime(1:4)) ! DateTime(1:4) = [year, month, day, zone]
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  An input argument `year = 0` corresponds to the historic 1 BC notation of the Gregorian calendar.<br>
    !>  This is in accordance with the convention in astronomical year numbering and the international standard date system, **ISO 8601**.<br>
    !>  The input `month` must be a number between `1` and `12`.<br>
    !>  The input `day` must be a number between `1` and `31` and be consistent with the specified month and (leap) year.<br>
    !>  The input `zone` must be a valid time zone in units of minutes.<br>
    !>  The input `hour` must be a number between `0` and `23`.<br>
    !>  The input `minute` must be a number between `0` and `59`.<br>
    !>  The input `second` must be a number between `0` and `59`.<br>
    !>  The input `millisecond` must be a number between `0` and `999`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  The procedures under this generic interface are always `impure` when all input arguments are missing.
    !>
    !>  \see
    !>  [getJulianDay](@ref pm_dateTime::getJulianDay)<br>
    !>  [getDateTimeNewZone](@ref pm_dateTime::getDateTimeNewZone)<br>
    !>
    !>  \example{getDateTimeUTC}
    !>  \include{lineno} example/pm_dateTime/getDateTimeUTC/main.F90
    !>  \compilef{getDateTimeUTC}
    !>  \output{getDateTimeUTC}
    !>  \include{lineno} example/pm_dateTime/getDateTimeUTC/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getDateTimeUTC}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 3:59 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface getDateTimeUTC

    impure module function getDateTimeUTCC() result(DateTimeUTC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeUTCC
#endif
        integer(IK)                                 :: DateTimeUTC(8)
    end function

    PURE module function getDateTimeUTCV(values) result(DateTimeUTC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeUTCV
#endif
        integer(IK)     , intent(in), contiguous    :: values(:)
        integer(IK)                                 :: DateTimeUTC(8)
    end function

!    PURE module function getDateTimeUTCY(year) result(DateTimeUTC)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeUTCY
!#endif
!        integer(IK)     , intent(in)                :: year
!        integer(IK)                                 :: DateTimeUTC(8)
!    end function
!
!    PURE module function getDateTimeUTCYM(year, month) result(DateTimeUTC)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeUTCYM
!#endif
!        integer(IK)     , intent(in)                :: year, month
!        integer(IK)                                 :: DateTimeUTC(8)
!    end function
!
!    PURE module function getDateTimeUTCYMD(year, month, day) result(DateTimeUTC)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeUTCYMD
!#endif
!        integer(IK)     , intent(in)                :: year, month, day
!        integer(IK)                                 :: DateTimeUTC(8)
!    end function

    PURE module function getDateTimeUTCYMDZ(year, month, day, zone) result(DateTimeUTC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeUTCYMDZ
#endif
        integer(IK)     , intent(in)                :: year, month, day, zone
        integer(IK)                                 :: DateTimeUTC(8)
    end function

    PURE module function getDateTimeUTCYMDZH(year, month, day, zone, hour) result(DateTimeUTC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeUTCYMDZH
#endif
        integer(IK)     , intent(in)                :: year, month, day, zone, hour
        integer(IK)                                 :: DateTimeUTC(8)
    end function

    PURE module function getDateTimeUTCYMDZHM(year, month, day, zone, hour, minute) result(DateTimeUTC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeUTCYMDZHM
#endif
        integer(IK)     , intent(in)                :: year, month, day, zone, hour, minute
        integer(IK)                                 :: DateTimeUTC(8)
    end function

    PURE module function getDateTimeUTCYMDZHMS(year, month, day, zone, hour, minute, second) result(DateTimeUTC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeUTCYMDZHMS
#endif
        integer(IK)     , intent(in)                :: year, month, day, zone, hour, minute, second
        integer(IK)                                 :: DateTimeUTC(8)
    end function

    PURE module function getDateTimeUTCYMDZHMSM(year, month, day, zone, hour, minute, second, millisecond) result(DateTimeUTC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeUTCYMDZHMSM
#endif
        integer(IK)     , intent(in)                :: year, month, day, zone, hour, minute, second, millisecond
        integer(IK)                                 :: DateTimeUTC(8)
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the current or the requested date and time in the requested time zone `newzone` (in units of minutes) as an integer-valued array of size `8`.<br>
    !>
    !>  \details
    !>  If there is no input argument besides `newzone`, then the output is the current date and time in the specified zone `newzone`.<br>
    !>  When the input argument `values(:)` is present, the first four input arguments or the first four elements of the input vector `values(:)`
    !>  are mandatory and must be provided since they determine the old date and time zone of interest to be converted to the new time zone.<br>
    !>  The rest of the elements containing the time of the day of the Gregorian Calendar are optional and taken to be zero if not provided.<br>
    !>
    !>  \param[in]  newzone     :   The input scalar of type `integer` of default kind \IK, containing the new <b>time zone</b> of the Gregorian calendar <b>in minutes</b>.<br>
    !>  \param[in]  year        :   The input scalar of type `integer` of default kind \IK, containing the year of the Gregorian calendar.<br>
    !>                              (**optional**, default =  the current value if all input arguments are missing. It can be present only if `values` is missing.)
    !>  \param[in]  month       :   The input scalar of type `integer` of default kind \IK, containing the month of the year of the Gregorian calendar.<br>
    !>                              (**optional**, default = `1` (or the current value if all input arguments are missing). It **must** be present if `year` is also present.)
    !>  \param[in]  day         :   The input scalar of type `integer` of default kind \IK, containing the day of the month of the year of the Gregorian calendar.<br>
    !>                              (**optional**, default = `1`. (or the current value if all input arguments are missing). It **must**  be present if `month` is also present.)
    !>  \param[in]  zone        :   The input scalar of type `integer` of default kind \IK, containing the <b>local time zone</b> of the Gregorian calendar <b>in minutes</b>
    !>                              to be used to convert the UTC date and time to the local date and time.<br>
    !>                              (**optional**, default = `0` (UTC) (or the current value if all input arguments are missing). It **must** be present if `day` is also present.)
    !>  \param[in]  hour        :   The input scalar of type `integer` of default kind \IK, containing the hour of the day of the Gregorian calendar.<br>
    !>                              (**optional**, default = `0` (or the current value if all input arguments are missing). It can be present only if `zone` is present.)
    !>  \param[in]  minute      :   The input scalar of type `integer` of default kind \IK, containing the minute of the hour of the day of the Gregorian calendar.<br>
    !>                              (**optional**, default = `0` (or the current value if all input arguments are missing). It can be present only if `hour` is present.)
    !>  \param[in]  second      :   The input scalar of type `integer` of default kind \IK, containing the second of the minute of the hour of the day of the Gregorian calendar.<br>
    !>                              (**optional**, default = `0` (or the current value if all input arguments are missing). It can be present only if `minute` is present.)
    !>  \param[in]  millisecond :   The input scalar of type `integer` of default kind \IK, containing the millisecond of the second of the minute of the hour of the day of the Gregorian calendar.<br>
    !>                              (**optional**, default = `0` (or the current value if all input arguments are missing). It can be present only if `second` is present.)
    !>  \param[in]  values      :   The input `contiguous` vector of minimum size `4` and maximum size `8` of type `integer` of default kind \IK, containing the values
    !>                              `[year, month, day, zone, hour, minute, seconds, milliseconds]` of the Gregorian calendar or a subset of the octuple starting with `year`.<br>
    !>                              The order of the elements of the vector follows that of the `values` returned by the Fortran intrinsic `date_and_time()`.<br>
    !>                              (**optional**, default = the current date and time. It can be present **only if** all other arguments except `newzone` are missing.)
    !>
    !>  \return
    !>  `DateTimeNewZone(1:8)`  :   The output vector of size `8` of type `integer` of default kind \IK containing
    !>                              the requested date and time of the Gregorian Calendar in the new time zone `newzone` corresponding
    !>                              to the input date and time in the order `[year, month, day, zone, hour, minute, seconds, milliseconds]`.
    !>                              By definition, `DateTimeNewZone(4)` (corresponding to the zone in minutes) has the same values as the input argument `newzone`.
    !>
    !>  \interface{getDateTimeNewZone}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: SK, IK
    !>      use pm_dateTime, only: getDateTimeNewZone
    !>      integer(IK) :: DateTimeNewZone(8)
    !>      integer(IK) :: DateTime(8)
    !>
    !>      DateTimeNewZone(1:8) = getDateTimeNewZone(newzone) ! return the current date and time in the specified zone `newzone`.
    !>      DateTimeNewZone(1:8) = getDateTimeNewZone(newzone, year, month, day, zone)
    !>      DateTimeNewZone(1:8) = getDateTimeNewZone(newzone, year, month, day, zone, hour)
    !>      DateTimeNewZone(1:8) = getDateTimeNewZone(newzone, year, month, day, zone, hour, minute)
    !>      DateTimeNewZone(1:8) = getDateTimeNewZone(newzone, year, month, day, zone, hour, minute, second)
    !>      DateTimeNewZone(1:8) = getDateTimeNewZone(newzone, year, month, day, zone, hour, minute, second, millisecond)
    !>      DateTimeNewZone(1:8) = getDateTimeNewZone(newzone, DateTime(1:8)) ! DateTime(1:8) = [year, month, day, zone, hour, minute, second, millisecond]
    !>      DateTimeNewZone(1:8) = getDateTimeNewZone(newzone, DateTime(1:7)) ! DateTime(1:7) = [year, month, day, zone, hour, minute, second]
    !>      DateTimeNewZone(1:8) = getDateTimeNewZone(newzone, DateTime(1:6)) ! DateTime(1:6) = [year, month, day, zone, hour, minute]
    !>      DateTimeNewZone(1:8) = getDateTimeNewZone(newzone, DateTime(1:5)) ! DateTime(1:5) = [year, month, day, zone, hour]
    !>      DateTimeNewZone(1:8) = getDateTimeNewZone(newzone, DateTime(1:4)) ! DateTime(1:4) = [year, month, day, zone]
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  An input argument `year = 0` corresponds to the historic 1 BC notation of the Gregorian calendar.<br>
    !>  This is in accordance with the convention in astronomical year numbering and the international standard date system, **ISO 8601**.<br>
    !>  The input `month` must be a number between `1` and `12`.<br>
    !>  The input `day` must be a number between `1` and `31` and be consistent with the specified month and (leap) year.<br>
    !>  The input `zone` must be a valid time zone in units of minutes.<br>
    !>  The input `hour` must be a number between `0` and `23`.<br>
    !>  The input `minute` must be a number between `0` and `59`.<br>
    !>  The input `second` must be a number between `0` and `59`.<br>
    !>  The input `millisecond` must be a number between `0` and `999`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  The procedures under this generic interface are always `impure` when all `optional` input arguments are missing.<br>
    !>
    !>  \see
    !>  [getDateTime](@ref pm_dateTime::getDateTime)<br>
    !>  [getJulianDay](@ref pm_dateTime::getJulianDay)<br>
    !>  [getDateTimeUTC](@ref pm_dateTime::getDateTimeUTC)<br>
    !>
    !>  \example{getDateTimeNewZone}
    !>  \include{lineno} example/pm_dateTime/getDateTimeNewZone/main.F90
    !>  \compilef{getDateTimeNewZone}
    !>  \output{getDateTimeNewZone}
    !>  \include{lineno} example/pm_dateTime/getDateTimeNewZone/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getDateTimeNewZone}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 3:59 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface getDateTimeNewZone

    impure module function getDateTimeNewZoneC(newzone) result(DateTimeNewZone)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeNewZoneC
#endif
        integer(IK)     , intent(in)                :: newzone
        integer(IK)                                 :: DateTimeNewZone(8)
    end function

    PURE module function getDateTimeNewZoneV(newzone, values) result(DateTimeNewZone)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeNewZoneV
#endif
        integer(IK)     , intent(in)                :: newzone
        integer(IK)     , intent(in), contiguous    :: values(:)
        integer(IK)                                 :: DateTimeNewZone(8)
    end function

!    PURE module function getDateTimeNewZoneY(year) result(DateTimeNewZone)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeNewZoneY
!#endif
!        integer(IK)     , intent(in)                :: newzone
!        integer(IK)     , intent(in)                :: year
!        integer(IK)                                 :: DateTimeNewZone(8)
!    end function
!
!    PURE module function getDateTimeNewZoneYM(year, month) result(DateTimeNewZone)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeNewZoneYM
!#endif
!        integer(IK)     , intent(in)                :: newzone
!        integer(IK)     , intent(in)                :: year, month
!        integer(IK)                                 :: DateTimeNewZone(8)
!    end function
!
!    PURE module function getDateTimeNewZoneYMD(year, month, day) result(DateTimeNewZone)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeNewZoneYMD
!#endif
!        integer(IK)     , intent(in)                :: newzone
!        integer(IK)     , intent(in)                :: year, month, day
!        integer(IK)                                 :: DateTimeNewZone(8)
!    end function

    PURE module function getDateTimeNewZoneYMDZ(newzone, year, month, day, zone) result(DateTimeNewZone)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeNewZoneYMDZ
#endif
        integer(IK)     , intent(in)                :: newzone
        integer(IK)     , intent(in)                :: year, month, day, zone
        integer(IK)                                 :: DateTimeNewZone(8)
    end function

    PURE module function getDateTimeNewZoneYMDZH(newzone, year, month, day, zone, hour) result(DateTimeNewZone)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeNewZoneYMDZH
#endif
        integer(IK)     , intent(in)                :: newzone
        integer(IK)     , intent(in)                :: year, month, day, zone, hour
        integer(IK)                                 :: DateTimeNewZone(8)
    end function

    PURE module function getDateTimeNewZoneYMDZHM(newzone, year, month, day, zone, hour, minute) result(DateTimeNewZone)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeNewZoneYMDZHM
#endif
        integer(IK)     , intent(in)                :: newzone
        integer(IK)     , intent(in)                :: year, month, day, zone, hour, minute
        integer(IK)                                 :: DateTimeNewZone(8)
    end function

    PURE module function getDateTimeNewZoneYMDZHMS(newzone, year, month, day, zone, hour, minute, second) result(DateTimeNewZone)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeNewZoneYMDZHMS
#endif
        integer(IK)     , intent(in)                :: newzone
        integer(IK)     , intent(in)                :: year, month, day, zone, hour, minute, second
        integer(IK)                                 :: DateTimeNewZone(8)
    end function

    PURE module function getDateTimeNewZoneYMDZHMSM(newzone, year, month, day, zone, hour, minute, second, millisecond) result(DateTimeNewZone)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeNewZoneYMDZHMSM
#endif
        integer(IK)     , intent(in)                :: newzone
        integer(IK)     , intent(in)                :: year, month, day, zone, hour, minute, second, millisecond
        integer(IK)                                 :: DateTimeNewZone(8)
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the current or the requested date and time as an integer-valued array of size `8` or a formatted string.<br>
    !>
    !>  \details
    !>  This generic interface performs and exceeds the functionalities of the well-known `strftime()` and `strptime()` in the C family of programming languages.<br>
    !>
    !>  -#  **When no input argument is provided**, the procedures under this generic interface return the
    !>      current date and time as an `integer` vector of `values` of size `8`, similar to the `values` argument of the Fortran intrinsic procedure `date_and_time()`.<br>
    !>  -#  **When the `format` argument is present**, the procedures under this generic interface return the current date and time as a string with the specified `format`.<br>
    !>  -#  When any components of the date and time are specified as input argument,<br>
    !>      -#  if `format` is missing, the procedures under this generic interface return the requested date and time as an `integer` vector of `values` of size `8` whose elements are set to the
    !>          corresponding specified arguments and all other elements corresponding to the missing arguments are set to zero.<br>
    !>      -#  if `format` is present, the procedures under this generic interface return the requested date and time as a formated string.<br>
    !>
    !>  The following specifier are recognized in the input `format` argument:<br>
    !>  (Except where mentioned, the convensions are identical to those of the `format` argument in the `strftime()` function of the C/C++ programming languages)<br>
    !>  (The specifiers marked with an asterisk (*) are locale-dependent in C/C++ but not currently in this library.)<br>
    !>
    !>  Specifier       |   Example                 |   Description
    !>  ----------------|---------------------------|-------------------------------------------------------------------------------------------------------------------
    !>  \f$\ms{%a}\f$   |   Thu                     |   Abbreviated weekday name *
    !>  \f$\ms{%A}\f$   |   Thursday                |   Full weekday name *
    !>  \f$\ms{%b}\f$   |   Aug                     |   Abbreviated month name *
    !>  \f$\ms{%B}\f$   |   August                  |   Full month name *
    !>  \f$\ms{%c}\f$   |   Thu Aug 23 14:55:02 2001|   Date and time representation *
    !>  \f$\ms{%C}\f$   |   20                      |   Year divided by 100 and floored to integer (-2,147,483,647 to +2,147,483,647) (extends the C/C++ `strftime` behavior)
    !>  \f$\ms{%d}\f$   |   23                      |   Day of the month, zero-padded (01-31)
    !>  \f$\ms{%D}\f$   |   08/23/01                |   Short `MM/DD/YY` date, equivalent to \f$\ms{%m/%d/%y}\f$
    !>  \f$\ms{%e}\f$   |   23                      |   Day of the month, space-padded ( 1-31)
    !>  \f$\ms{%F}\f$   |   2001-08-23              |   Short `YYYY-MM-DD` date, equivalent to \f$\ms{%Y-%m-%d}\f$
    !>  \f$\ms{%f}\f$   |   037                     |   The millisecond padded with leading zeros (**unique to this library**, follows the convension of Python for microsecond).
    !>  \f$\ms{%g}\f$   |   01                      |   Week-based year, last two digits (00-99)
    !>  \f$\ms{%G}\f$   |   2001                    |   Week-based year
    !>  \f$\ms{%h}\f$   |   Aug                     |   Abbreviated month name * (same as \f$\ms{%b}\f$)
    !>  \f$\ms{%H}\f$   |   14                      |   Hour in 24h format (00-23)
    !>  \f$\ms{%I}\f$   |   02                      |   Hour in 12h format (01-12)
    !>  \f$\ms{%j}\f$   |   235                     |   Day of the year (001-366)
    !>  \f$\ms{%m}\f$   |   08                      |   Month as a decimal number (01-12)
    !>  \f$\ms{%M}\f$   |   55                      |   Minute (00-59)
    !>  \f$\ms{%n}\f$   |   \f$\ms{\n}\f$           |   New-line character
    !>  \f$\ms{%p}\f$   |   PM                      |   `AM` or `PM` designation
    !>  \f$\ms{%r}\f$   |   02:55:02 pm             |   12-hour clock time *
    !>  \f$\ms{%R}\f$   |   14:55                   |   24-hour `HH:MM` time, equivalent to \f$\ms{%H:%M}\f$
    !>  \f$\ms{%S}\f$   |   02                      |   Second (00-59)
    !>  \f$\ms{%t}\f$   |   \f$\ms{\t}\f$           |   Horizontal-tab character
    !>  \f$\ms{%T}\f$   |   14:55:02                |   ISO 8601 time format (`HH:MM:SS`), equivalent to \f$\ms{%H:%M:%S}\f$
    !>  \f$\ms{%u}\f$   |   4                       |   ISO 8601 weekday as number with Monday as 1 (1-7)
    !>  \f$\ms{%U}\f$   |   33                      |   Week number with the first Sunday as the first day of week one (00-53) (**currently not implemented**)
    !>  \f$\ms{%V}\f$   |   34                      |   ISO 8601 week number (01-53)
    !>  \f$\ms{%w}\f$   |   4                       |   Weekday as a decimal number with Sunday as 0 (0-6)
    !>  \f$\ms{%W}\f$   |   34                      |   Week number with the first Monday as the first day of week one (00-53) (**currently not implemented**)
    !>  \f$\ms{%x}\f$   |   08/23/01                |   Date representation *
    !>  \f$\ms{%X}\f$   |   14:55:02                |   Time representation *
    !>  \f$\ms{%y}\f$   |   01                      |   Year, last two digits (00-99)
    !>  \f$\ms{%Y}\f$   |   2001                    |   Year
    !>  \f$\ms{%z}\f$   |   +0100                   |   ISO 8601 offset from UTC in time zone in units of minutes (as returned by the Fortran intrinsic `date_and_time()`).
    !>  \f$\ms{%Z}\f$   |   CDT                     |   Time zone abbreviation (**daylight saving not implemented**; different from the C/C++ `strftime` behavior which uses system time zone name)
    !>  \f$\ms{%%}\f$   |   \%                      |   A `%` sign
    !>
    !>  \param[in]  format      :   The input scalar or array of the same shape as other array-like arguments of type `character` of default kind \SK,
    !>                              containing the output format of the specified date and time (See the table above for a list of possible specifiers).
    !>                              (**optional**. If missing, the output date and time will be an `integer` vector of size `8`.)
    !>  \param[in]  julianDay   :   The input scalar of type `real` of default kind \RK, containing the Julian Day to be converted to the corresponding Gregorian date and time.<br>
    !>                              (**optional**, default = it can be present only if all arguments are missing or if only `format` and/or `zone` are present.)
    !>  \param[in]  year        :   The input scalar of type `integer` of default kind \IK, containing the year of the Gregorian calendar.<br>
    !>                              (**optional**, default =  the current value if all arguments are missing or only `format` is present.)
    !>  \param[in]  month       :   The input scalar of type `integer` of default kind \IK, containing the month of the year of the Gregorian calendar.<br>
    !>                              (**optional**, default = `1` or the current value if all arguments are missing or only `format` is present. It can be present only if `year` is present.)
    !>  \param[in]  day         :   The input scalar of type `integer` of default kind \IK, containing the day of the month of the year of the Gregorian calendar.<br>
    !>                              (**optional**, default = `1` or the current value if all arguments are missing or only `format` is present. It can be present only if `month` is present.)
    !>  \param[in]  zone        :   The input scalar of type `integer` of default kind \IK, containing the <b>time zone</b> of the Gregorian calendar <b>in minutes</b>.<br>
    !>                              (**optional**, default = `0` (UTC) or the current value if all arguments are missing or only `format` is present. It can be present only if either `day` or `julianDay` (but noth both) is present.)
    !>  \param[in]  hour        :   The input scalar of type `integer` of default kind \IK, containing the hour of the day of the Gregorian calendar.<br>
    !>                              (**optional**, default = `0` or the current value if all arguments are missing or only `format` is present. It can be present only if `zone` is present.)
    !>  \param[in]  minute      :   The input scalar of type `integer` of default kind \IK, containing the minute of the hour of the day of the Gregorian calendar.<br>
    !>                              (**optional**, default = `0` or the current value if all arguments are missing or only `format` is present. It can be present only if `hour` is present.)
    !>  \param[in]  second      :   The input scalar of type `integer` of default kind \IK, containing the second of the minute of the hour of the day of the Gregorian calendar.<br>
    !>                              (**optional**, default = `0` or the current value if all arguments are missing or only `format` is present. It can be present only if `minute` is present.)
    !>  \param[in]  millisecond :   The input scalar of type `integer` of default kind \IK, containing the millisecond of the second of the minute of the hour of the day of the Gregorian calendar.<br>
    !>                              (**optional**, default = `0` or the current value if all arguments are missing or only `format` is present. It can be present only if `second` is present.)
    !>  \param[in]  values      :   The input `contiguous` vector of maximum size `8` of type `integer` of default kind \IK, containing
    !>                              the values `[year, month, day, zone, hour, minute, seconds, milliseconds]` of the Gregorian calendar.<br>
    !>                              The order of the elements of the vector follows that of the `values` returned by the Fortran intrinsic `date_and_time()`.<br>
    !>                              (**optional**, default = the current date and time. It can be present **only if** `format` is also present and all other arguments are missing.)
    !>
    !>  \return
    !>  `string / values(1:8)`  :   If the input argument `format` is present, the output is a scalar formatted string containing the requested date and time.<br>
    !>                              If the input argument `format` is missing, the output is a vector of size `8` of type `integer` of default kind \IK containing
    !>                              the requested date and time of the Gregorian Calendar.
    !>
    !>  \interface{getDateTime}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getDateTime
    !>      use pm_kind, only: SK, IK
    !>      character(:, SK), allocatable :: string
    !>      integer(IK) :: values(8)
    !>
    !>      values(1:8) = getDateTime() ! return the current date and time.
    !>      values(1:8) = getDateTime(year)
    !>      values(1:8) = getDateTime(year, month)
    !>      values(1:8) = getDateTime(year, month, day)
    !>      values(1:8) = getDateTime(year, month, day, zone)
    !>      values(1:8) = getDateTime(year, month, day, zone, hour)
    !>      values(1:8) = getDateTime(year, month, day, zone, hour, minute)
    !>      values(1:8) = getDateTime(year, month, day, zone, hour, minute, second)
    !>      values(1:8) = getDateTime(year, month, day, zone, hour, minute, second, millisecond)
    !>      values(1:8) = getDateTime(julianDay, zone) ! convert the Julian Day to Gregorian date and time in the requested zone (in minutes) with respect to UTC.
    !>      values(1:8) = getDateTime(julianDay) ! convert the Julian Day to the corresponding UTC Gregorian date and time.
    !>      string = getDateTime(format, values(:)) ! return the current `values` date and time as a formatted string.
    !>      string = getDateTime(format) ! impure : return the current local date and time with the requested format.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < size(values)` must hold for the corresponding input arguments.<br>
    !>  The condition `9 > size(values)` must hold for the corresponding input arguments.<br>
    !>  The specified number of elements of `values`, if less than `8`, must be consistent with the specifiers in input `format`.<br>
    !>  For example, when \f$\ms{%Z}\f$ is specified together with `values`, then the size of `values` must be at least `4`.<br>
    !>  An input argument `year = 0` corresponds to the historic 1 BC notation of the Gregorian calendar.<br>
    !>  This is in accordance with the convention in astronomical year numbering and the international standard date system, **ISO 8601**.<br>
    !>  The input `month` must be a number between `1` and `12`.<br>
    !>  The input `day` must be a number between `1` and `31` and be consistent with the specified month and (leap) year.<br>
    !>  The input `zone` must be a valid time zone in units of minutes.<br>
    !>  The input `hour` must be a number between `0` and `23`.<br>
    !>  The input `minute` must be a number between `0` and `59`.<br>
    !>  The input `second` must be a number between `0` and `59`.<br>
    !>  The input `millisecond` must be a number between `0` and `999`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  The procedures under this generic interface are always `impure` when all input arguments are missing.
    !>
    !>  \see
    !>  [getMillisecond](@ref pm_dateTime::getMillisecond)<br>
    !>  [getSecond](@ref pm_dateTime::getSecond)<br>
    !>  [getMinute](@ref pm_dateTime::getMinute)<br>
    !>  [getHour12](@ref pm_dateTime::getHour12)<br>
    !>  [getHour](@ref pm_dateTime::getHour)<br>
    !>  [getZone](@ref pm_dateTime::getZone)<br>
    !>  [getZoneAbbr](@ref pm_dateTime::getZoneAbbr)<br>
    !>  [getDay](@ref pm_dateTime::getDay)<br>
    !>  [getMonth](@ref pm_dateTime::getMonth)<br>
    !>  [getYear](@ref pm_dateTime::getYear)<br>
    !>  [isLeapYear](@ref pm_dateTime::isLeapYear)<br>
    !>  [getJulianDay](@ref pm_dateTime::getJulianDay)<br>
    !>  [getDateTimeNewZone](@ref pm_dateTime::getDateTimeNewZone)<br>
    !>
    !>  \example{getDateTime}
    !>  \include{lineno} example/pm_dateTime/getDateTime/main.F90
    !>  \compilef{getDateTime}
    !>  \output{getDateTime}
    !>  \include{lineno} example/pm_dateTime/getDateTime/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \todo
    !>  \plow
    !>  Locale-dependent date and time formats can be added to this generic interface in the future, depending on the demand.<br>
    !>
    !>  \todo
    !>  \pmed
    !>  This interface can be extended to add `format` input argument to all possible calling interfaces.<br>
    !>
    !>  \todo
    !>  \phigh
    !>  The format specifiers \f$\ms{%U}\f$ and \f$\ms{%W}\f$ must be implemented.<br>
    !>
    !>  \final{getDateTime}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface getDateTime

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module function getDateTimeValuesJ(julianDay) result(values)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeValuesJ
#endif
        real(RK)        , intent(in)                :: julianDay
        integer(IK)                                 :: values(8)
    end function

    PURE module function getDateTimeValuesJZ(julianDay, zone) result(values)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeValuesJZ
#endif
        real(RK)        , intent(in)                :: julianDay
        integer(IK)     , intent(in)                :: zone
        integer(IK)                                 :: values(8)
    end function

    impure module function getDateTimeValuesC() result(values)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeValuesC
#endif
        integer(IK)                                 :: values(8)
    end function

    PURE module function getDateTimeValuesY(year) result(values)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeValuesY
#endif
        integer(IK)     , intent(in)                :: year
        integer(IK)                                 :: values(8)
    end function

    PURE module function getDateTimeValuesYM(year, month) result(values)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeValuesYM
#endif
        integer(IK)     , intent(in)                :: year, month
        integer(IK)                                 :: values(8)
    end function

    PURE module function getDateTimeValuesYMD(year, month, day) result(values)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeValuesYMD
#endif
        integer(IK)     , intent(in)                :: year, month, day
        integer(IK)                                 :: values(8)
    end function

    PURE module function getDateTimeValuesYMDZ(year, month, day, zone) result(values)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeValuesYMDZ
#endif
        integer(IK)     , intent(in)                :: year, month, day, zone
        integer(IK)                                 :: values(8)
    end function

    PURE module function getDateTimeValuesYMDZH(year, month, day, zone, hour) result(values)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeValuesYMDZH
#endif
        integer(IK)     , intent(in)                :: year, month, day, zone, hour
        integer(IK)                                 :: values(8)
    end function

    PURE module function getDateTimeValuesYMDZHM(year, month, day, zone, hour, minute) result(values)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeValuesYMDZHM
#endif
        integer(IK)     , intent(in)                :: year, month, day, zone, hour, minute
        integer(IK)                                 :: values(8)
    end function

    PURE module function getDateTimeValuesYMDZHMS(year, month, day, zone, hour, minute, second) result(values)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeValuesYMDZHMS
#endif
        integer(IK)     , intent(in)                :: year, month, day, zone, hour, minute, second
        integer(IK)                                 :: values(8)
    end function

    PURE module function getDateTimeValuesYMDZHMSM(year, month, day, zone, hour, minute, second, millisecond) result(values)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeValuesYMDZHMSM
#endif
        integer(IK)     , intent(in)                :: year, month, day, zone, hour, minute, second, millisecond
        integer(IK)                                 :: values(8)
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function getDateTimeStringC(format) result(string)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeStringC
#endif
        character(*, SK), intent(in)                :: format
        character(:, SK), allocatable               :: string
    end function

    PURE module function getDateTimeStringV(format, values) result(string)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeStringV
#endif
        character(*, SK), intent(in)                :: format
        integer(IK)     , intent(in), contiguous    :: values(:)
        character(:, SK), allocatable               :: string
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **ISO 8601 Week Date** triple `[week year, week number, week day]` corresponding to the input
    !>  Gregorian date and time `values(1:3)` or `(year, month, day)` triple.
    !>
    !>  \details
    !>  The **ISO week date system** is a *leap week calendar* system that is part of the \f$\ms{ISO 8601}\f$
    !>  date and time standard issued by the International Organization for Standardization (**ISO**) since 1988
    !>  It is (mainly) used in government and business for fiscal years, as well as in timekeeping.<br>
    !>  The system specifies a week year atop the Gregorian calendar by defining a notation for **ordinal weeks** of the year.<br>
    !>  The following table contains examples of the Gregorian Calendar date and their corresponding Week Dates.
    !>
    !>  Gregorian Date  |   ISO Date    |   ISO Week Date
    !>  ----------------|---------------|----------------
    !>  Sat 1  Jan 1977 |   1977-01-01  |   1976-W53-6
    !>  Sun 2  Jan 1977 |   1977-01-02  |   1976-W53-7
    !>  Sat 31 Dec 1977 |   1977-12-31  |   1977-W52-6
    !>  Sun 1  Jan 1978 |   1978-01-01  |   1977-W52-7
    !>  Mon 2  Jan 1978 |   1978-01-02  |   1978-W01-1
    !>  Sun 31 Dec 1978 |   1978-12-31  |   1978-W52-7
    !>  Mon 1  Jan 1979 |   1979-01-01  |   1979-W01-1
    !>  Sun 30 Dec 1979 |   1979-12-30  |   1979-W52-7
    !>  Mon 31 Dec 1979 |   1979-12-31  |   1980-W01-1
    !>  Tue 1  Jan 1980 |   1980-01-01  |   1980-W01-2
    !>  Sun 28 Dec 1980 |   1980-12-28  |   1980-W52-7
    !>  Mon 29 Dec 1980 |   1980-12-29  |   1981-W01-1
    !>  Tue 30 Dec 1980 |   1980-12-30  |   1981-W01-2
    !>  Wed 31 Dec 1980 |   1980-12-31  |   1981-W01-3
    !>  Thu 1  Jan 1981 |   1981-01-01  |   1981-W01-4
    !>  Thu 31 Dec 1981 |   1981-12-31  |   1981-W53-4
    !>  Fri 1  Jan 1982 |   1982-01-01  |   1981-W53-5
    !>  Sat 2  Jan 1982 |   1982-01-02  |   1981-W53-6
    !>  Sun 3  Jan 1982 |   1982-01-03  |   1981-W53-7
    !>
    !>  When all input arguments are missing, the procedures under this generic interface return the current ISO Week Date.
    !>
    !>  \param[in]  values  :   The input `contiguous` array of shape `(:)`, of size `3` or larger, of type `integer` of default kind \IK, containing the `[year, month, day]` triple of the Gregorian calendar.<br>
    !>                          For the current local date, this triple can be obtained from the Fortran intrinsic `date_and_time()` or [getDateTime()](@ref pm_dateTime::getDateTime).<br>
    !>                          Only the first three elements (`values(1:3)`) are used to compute the output.<br>
    !>                          The ability to pass longer vectors as input is to allow the output `values(1:8)` of various
    !>                          functionalities of this module to be passed directly to the procedures under this generic interface.<br>
    !>                          (**optional**. It **can** be present <b>if and only if</b> all other input arguments are missing.)
    !>  \param[in]  year    :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the year of the Gregorian calendar.<br>
    !>                          (**optional**. It **can** be present <b>if and only if</b> the input argument `values` is missing.)
    !>  \param[in]  month   :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the month of the Gregorian calendar.<br>
    !>                          (**optional**. It **must** be present <b>if and only if</b> the input argument `year` is present.)
    !>  \param[in]  day     :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the day of the Gregorian calendar.<br>
    !>                          (**optional**. It **must** be present <b>if and only if</b> the input argument `month` is present.)
    !>
    !>  \return
    !>  `WeekDate`          :   The output array of shape `(1:3)` of type `integer` of default kind \IK containing the **ISO 8601 Week Date** triple `[week year, week number, week day]`.
    !>
    !>  \interface{getWeekDate}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getWeekDate
    !>      integer(IK) :: WeekDate(3)
    !>
    !>      WeekDate(1:3) = getWeekDate() ! always impure.
    !>      WeekDate(1:3) = getWeekDate(values(1:3))
    !>      WeekDate(1:3) = getWeekDate(year, month, day)
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The size of the input argument `values(:)` must be at least `3` and at most `8`.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>  The procedures under this generic interface are always `impure` when all input arguments are missing.
    !>
    !>  \see
    !>  [getWeekYear](@ref pm_dateTime::getWeekYear)<br>
    !>  [getDateTime](@ref pm_dateTime::getDateTime)<br>
    !>  [getDateAfter](@ref pm_dateTime::getDateAfter)<br>
    !>  [getDateBefore](@ref pm_dateTime::getDateBefore)<br>
    !>
    !>  \example{getWeekDate}
    !>  \include{lineno} example/pm_dateTime/getWeekDate/main.F90
    !>  \compilef{getWeekDate}
    !>  \output{getWeekDate}
    !>  \include{lineno} example/pm_dateTime/getWeekDate/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getWeekDate}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface getWeekDate

    impure module function getWeekDateCurrent() result(WeekDate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getWeekDateCurrent
#endif
        integer(IK)                                 :: WeekDate(3)
    end function

    PURE module function getWeekDateValues(values) result(WeekDate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getWeekDateValues
#endif
        integer(IK)     , intent(in), contiguous    :: values(:)
        integer(IK)                                 :: WeekDate(3)
    end function

    PURE module function getWeekDateTriple(year, month, day) result(WeekDate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getWeekDateTriple
#endif
        integer(IK)     , intent(in)                :: year, month, day
        integer(IK)                                 :: WeekDate(3)
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **week year** corresponding to the **ISO 8601 Week Date** triple `[week year, week number, week day]`
    !>  equivalent of the input Gregorian date and time `values(1:3)` or `(year, month, day)` triple.
    !>
    !>  \details
    !>  See the documentation of [getWeekDate](@ref pm_dateTime::getWeekDate) for more information on the ISO week date system.<br>
    !>
    !>  \note
    !>  When all input arguments are missing, the procedures under this generic interface return the current ISO Week Year.
    !>
    !>  \param[in]  values  :   The input `contiguous` array of shape `(:)`, of size `3` or larger, of type `integer` of default kind \IK, containing the `[year, month, day]` triple of the Gregorian calendar.<br>
    !>                          For the current local date, this triple can be obtained from the Fortran intrinsic `date_and_time()` or [getDateTime()](@ref pm_dateTime::getDateTime).<br>
    !>                          Only the first three elements (`values(1:3)`) are used to compute the output.<br>
    !>                          The ability to pass longer vectors as input is to allow the output `values(1:8)` of various
    !>                          functionalities of this module to be passed directly to the procedures under this generic interface.<br>
    !>                          (**optional**. It **can** be present <b>if and only if</b> all other input arguments are missing.)
    !>  \param[in]  year    :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the year of the Gregorian calendar.<br>
    !>                          (**optional**. It **can** be present <b>if and only if</b> the input argument `values` is missing.)
    !>  \param[in]  month   :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the month of the Gregorian calendar.<br>
    !>                          (**optional**. It **must** be present <b>if and only if</b> the input argument `year` is present.)
    !>  \param[in]  day     :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the day of the Gregorian calendar.<br>
    !>                          (**optional**. It **must** be present <b>if and only if</b> the input argument `month` is present.)
    !>
    !>  \return
    !>  `weekYear`          :   The output scalar (or array of the same rank as `year`, `month`, or day`) of shape `(1:3)` of type `integer` of default kind \IK containing the <b>ISO 8601 Week Year</b>.
    !>
    !>  \interface{getWeekYear}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getWeekYear
    !>      integer(IK) :: WeekDate(3)
    !>
    !>      WeekDate(1:3) = getWeekYear() ! always impure.
    !>      WeekDate(1:3) = getWeekYear(values(1:3))
    !>      WeekDate(1:3) = getWeekYear(year, month, day)
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The size of the input argument `values(:)` must be at least `3` and at most `8`.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>  The procedures under this generic interface are always `impure` when all input arguments are missing.<br>
    !>
    !>  \remark
    !>  The procedures under this generic interface are `elemental` when the input arguments `year, month, day` are present.<br>
    !>
    !>  \see
    !>  [getWeekDate](@ref pm_dateTime::getWeekDate)<br>
    !>  [getDateTime](@ref pm_dateTime::getDateTime)<br>
    !>  [getDateAfter](@ref pm_dateTime::getDateAfter)<br>
    !>  [getDateBefore](@ref pm_dateTime::getDateBefore)<br>
    !>
    !>  \example{getWeekYear}
    !>  \include{lineno} example/pm_dateTime/getWeekYear/main.F90
    !>  \compilef{getWeekYear}
    !>  \output{getWeekYear}
    !>  \include{lineno} example/pm_dateTime/getWeekYear/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getWeekYear}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface getWeekYear

    impure module function getWeekYearCurrent() result(weekYear)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getWeekYearCurrent
#endif
        integer(IK)                                 :: weekYear
    end function

    PURE module function getWeekYearValues(values) result(weekYear)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getWeekYearValues
#endif
        integer(IK)     , intent(in), contiguous    :: values(:)
        integer(IK)                                 :: weekYear
    end function

    PURE elemental module function getWeekYearTriple(year, month, day) result(weekYear)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getWeekYearTriple
#endif
        integer(IK)     , intent(in)                :: year, month, day
        integer(IK)                                 :: weekYear
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input date and time `values(:)` or `(year, month, day, zone, hour, second, millisecond)` octuple
    !>  corresponds to a valid Gregorian Calendar date and time.
    !>
    !>  \details
    !>  Returning the correct result requires taking into account the possibility of leap years and the varying day counts of months.<br>
    !>  A valid date requires,
    !>  -#  a `month` between `1` and `12`,
    !>  -#  a `day` between `1` and `31` and consistent with the specified month and year,
    !>  -#  a `zone` must be between `-12:00` and `+14:00` UTC, corresponding to `-720` and `+840` minutes as used here.
    !>  -#  an `hour` of day between `0` and `23`,
    !>  -#  a `minute` of hour between `0` and `59`,
    !>  -#  a `second` of minute between `0` and `59`,
    !>  -#  a `millisecond` of second between `0` and `999`.
    !>
    !>  \param[in]  values      :   The input `contiguous` array of shape `(:)` of type `integer` of default kind \IK, containing the
    !>                              octuple `[year, month, day, zone, hour, minute, second, millisecond]` or a subset of it that always starts with `year`
    !>                              representing a Gregorian calendar date and/or time.<br>
    !>                              The order of the elements of `values` follows that of the `values` argument of the Fortran intrinsic `date_and_time()`.<br>
    !>                              (**optional**. It must be present <b>if and only if</b> all other input arguments are missing.)
    !>  \param[in]  year        :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing a year of the Gregorian calendar.<br>
    !>                              (**optional**, default = any value that yields a valid date given the other input argument. It can be present only if the input argument `values` is missing.)
    !>  \param[in]  month       :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing a month of a year of the Gregorian calendar.<br>
    !>                              (**optional**, default = any value that yields a valid date given the other input argument. It can be present only if the input argument `year` is present.)
    !>  \param[in]  day         :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing a day of a month of a year of the Gregorian calendar.<br>
    !>                              (**optional**, default = any value that yields a valid date given the other input argument. It can be present only if the input argument `month` is present.)
    !>  \param[in]  zone        :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing a time zone of the Gregorian calendar.<br>
    !>                              (**optional**, default = any value that yields a valid date given the other input argument. It can be present only if the input argument `day` is present.)
    !>  \param[in]  hour        :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the hour of a day of the Gregorian calendar.<br>
    !>                              (**optional**, default = any value that yields a valid date given the other input argument. It can be present only if the input argument `zone` is present.)
    !>  \param[in]  minute      :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the minute of the hour of a day of the Gregorian calendar.<br>
    !>                              (**optional**, default = any value that yields a valid date given the other input argument. It can be present only if the input argument `hour` is present.)
    !>  \param[in]  second      :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the second of the minute of the hour of a day of the Gregorian calendar.<br>
    !>                              (**optional**, default = any value that yields a valid date given the other input argument. It can be present only if the input argument `minute` is present.)
    !>  \param[in]  millisecond :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the millisecond of the second of the minute of the hour of a day of the Gregorian calendar.<br>
    !>                              (**optional**, default = any value that yields a valid date given the other input argument. It can be present only if the input argument `second` is present.)
    !>
    !>  \return
    !>  `isValid`               :   The output scalar or array of the same shape as the input array-like arguments (except `Vector`), of type `logical` of default kind \LK.
    !>                              It is `.true.` <b>if and only if</b> the input argument(s) correspond to a valid date and/or month of the Gregorian Calendar.<br>
    !>                              Otherwise, it is `.false.` is the input date and time is invalid or if the condition `0 < size(values) < 9` does not hold.
    !>
    !>  \interface{isValidDateTime}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: isValidDateTime
    !>      logical(LK) :: isValid
    !>
    !>      isValid = isValidDateTime(year) ! elemental
    !>      isValid = isValidDateTime(year, month) ! elemental
    !>      isValid = isValidDateTime(year, month, day) ! elemental
    !>      isValid = isValidDateTime(year, month, day, zone) ! elemental
    !>      isValid = isValidDateTime(year, month, day, zone, hour) ! elemental
    !>      isValid = isValidDateTime(year, month, day, zone, hour, minute) ! elemental
    !>      isValid = isValidDateTime(year, month, day, zone, hour, minute, second) ! elemental
    !>      isValid = isValidDateTime(year, month, day, zone, hour, minute, second, millisecond) ! elemental
    !>      isValid = isValidDateTime(values(1:1)) ! values = [year]
    !>      isValid = isValidDateTime(values(1:2)) ! values = [year, month]
    !>      isValid = isValidDateTime(values(1:3)) ! values = [year, month, day]
    !>      isValid = isValidDateTime(values(1:4)) ! values = [year, month, day, zone]
    !>      isValid = isValidDateTime(values(1:5)) ! values = [year, month, day, zone, hour]
    !>      isValid = isValidDateTime(values(1:6)) ! values = [year, month, day, zone, hour, minute]
    !>      isValid = isValidDateTime(values(1:7)) ! values = [year, month, day, zone, hour, minute, second]
    !>      isValid = isValidDateTime(values(1:8)) ! values = [year, month, day, zone, hour, minute, second, millisecond]
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  An input argument `year = 0` corresponds to the historic 1 BC notation of the Gregorian calendar.<br>
    !>  This is in accordance with the convention in astronomical year numbering and the international standard date system, **ISO 8601**.<br>
    !>  In these systems, the year `0` is a leap year.
    !>
    !>  \warning
    !>  The size of the input argument `values(:)` must be non-zero and less than `8`, otherwise, the returned value is `.false.`.<br>
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>  The procedures under this generic interface are non-elemental when the input argument `values(:)` is present.
    !>
    !>  \see
    !>  [getDateAfter](@ref pm_dateTime::getDateAfter)<br>
    !>  [getDateTime](@ref pm_dateTime::getDateTime)<br>
    !>
    !>  \example{isValidDateTime}
    !>  \include{lineno} example/pm_dateTime/isValidDateTime/main.F90
    !>  \compilef{isValidDateTime}
    !>  \output{isValidDateTime}
    !>  \include{lineno} example/pm_dateTime/isValidDateTime/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{isValidDateTime}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface isValidDateTime

    PURE module function isValidDateTimeV(values) result(isValid)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isValidDateTimeV
#endif
        integer(IK)     , intent(in), contiguous    :: values(:)
        logical(LK)                                 :: isValid
    end function

    PURE elemental module function isValidDateTimeY(year) result(isValid)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isValidDateTimeY
#endif
        integer(IK)     , intent(in)                :: year
        logical(LK)                                 :: isValid
    end function

    PURE elemental module function isValidDateTimeYM(year, month) result(isValid)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isValidDateTimeYM
#endif
        integer(IK)     , intent(in)                :: year, month
        logical(LK)                                 :: isValid
    end function

    PURE elemental module function isValidDateTimeYMD(year, month, day) result(isValid)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isValidDateTimeYMD
#endif
        integer(IK)     , intent(in)                :: year, month, day
        logical(LK)                                 :: isValid
    end function

    PURE elemental module function isValidDateTimeYMDZ(year, month, day, zone) result(isValid)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isValidDateTimeYMDZ
#endif
        integer(IK)     , intent(in)                :: year, month, day, zone
        logical(LK)                                 :: isValid
    end function

    PURE elemental module function isValidDateTimeYMDZH(year, month, day, zone, hour) result(isValid)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isValidDateTimeYMDZH
#endif
        integer(IK)     , intent(in)                :: year, month, day, zone, hour
        logical(LK)                                 :: isValid
    end function

    PURE elemental module function isValidDateTimeYMDZHM(year, month, day, zone, hour, minute) result(isValid)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isValidDateTimeYMDZHM
#endif
        integer(IK)     , intent(in)                :: year, month, day, zone, hour, minute
        logical(LK)                                 :: isValid
    end function

    PURE elemental module function isValidDateTimeYMDZHMS(year, month, day, zone, hour, minute, second) result(isValid)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isValidDateTimeYMDZHMS
#endif
        integer(IK)     , intent(in)                :: year, month, day, zone, hour, minute, second
        logical(LK)                                 :: isValid
    end function

    PURE elemental module function isValidDateTimeYMDZHMSM(year, month, day, zone, hour, minute, second, millisecond) result(isValid)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isValidDateTimeYMDZHMSM
#endif
        integer(IK)     , intent(in)                :: year, month, day, zone, hour, minute, second, millisecond
        logical(LK)                                 :: isValid
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input date (`year`, `month`, `day`) triple corresponds to the last day of a Gregorian Calendar month.
    !>
    !>  \details
    !>  Returning the correct result requires taking into account the possibility of leap years.
    !>
    !>  \param[in]  values  :   The input `contiguous` array of shape `(:)` of size `3` or larger, of type `integer` of default kind \IK, containing the `[year, month, day]` triple of the Gregorian calendar.<br>
    !>                          For the current local date, this triple can be obtained from the Fortran intrinsic `date_and_time()` or [getDateTime()](@ref pm_dateTime::getDateTime).<br>
    !>                          Only the first three elements (`values(1:3)`) are used to compute the output.<br>
    !>                          The ability to pass longer vectors as input is to allow the output `values(1:8)` of various
    !>                          functionalities of this module to be passed directly to the procedures under this generic interface.<br>
    !>                          (**optional**. It **can** be present <b>if and only if</b> all other input arguments are missing.)
    !>  \param[in]  year    :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the year of the Gregorian calendar.<br>
    !>                          (**optional**. It **can** be present <b>if and only if</b> the input argument `values` is missing.)
    !>  \param[in]  month   :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the month of the Gregorian calendar.<br>
    !>                          (**optional**. It **must** be present <b>if and only if</b> the input argument `year` is present.)
    !>  \param[in]  day     :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the day of the Gregorian calendar.<br>
    !>                          (**optional**. It **must** be present <b>if and only if</b> the input argument `month` is present.)
    !>
    !>  \return
    !>  `lastDayInMonth`    :   The output scalar or array of the same shape as the input array-like arguments of type `logical` of default kind \LK.
    !>                          It is `.true.` <b>if and only if</b> the input date corresponds to the last day of a month in the Gregorian Calendar.
    !>
    !>  \interface{isLastDayInMonth}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: isLastDayInMonth
    !>      logical(LK) :: lastDayInMonth
    !>
    !>      lastDayInMonth = isLastDayInMonth() ! current date is used.
    !>      lastDayInMonth = isLastDayInMonth(values(1:3)) ! values = [year, month, day]
    !>      lastDayInMonth = isLastDayInMonth(year, month, day)
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The size of the input argument `values(:)` must be at least `3` and at most `8`.<br>
    !>  The input values for the `year`, `month`, and `day` must be valid values.<br>
    !>  The input `month` must be a number between `1` and `12`.<br>
    !>  The input `day` must be a number between `1` and `31`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>  The procedures under this generic interface are non-elemental when the input argument `values(:)` is present.
    !>
    !>  \see
    !>  [getDateAfter](@ref pm_dateTime::getDateAfter)<br>
    !>
    !>  \example{isLastDayInMonth}
    !>  \include{lineno} example/pm_dateTime/isLastDayInMonth/main.F90
    !>  \compilef{isLastDayInMonth}
    !>  \output{isLastDayInMonth}
    !>  \include{lineno} example/pm_dateTime/isLastDayInMonth/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{isLastDayInMonth}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface isLastDayInMonth

    module function isLastDayInMonthC() result(lastDayInMonth)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isLastDayInMonthC
#endif
        logical(LK)                         :: lastDayInMonth
    end function

    PURE module function isLastDayInMonthValues(values) result(lastDayInMonth)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isLastDayInMonthValues
#endif
        integer(IK), intent(in), contiguous :: values(:)
        logical(LK)                         :: lastDayInMonth
    end function

    PURE elemental module function isLastDayInMonthTriple(year, month, day) result(lastDayInMonth)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isLastDayInMonthTriple
#endif
        integer(IK), intent(in) :: year, month, day
        logical(LK)             :: lastDayInMonth
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the date in the Gregorian Calendar that appears after the input date.
    !>
    !>  \details
    !>  Returning the correct result requires taking into account the possibility of leap years.
    !>
    !>  \param[in]  values  :   The input `contiguous` array of shape `(:)`, of size `3` or larger, of type `integer` of default kind \IK, containing the `[year, month, day]` triple of the Gregorian calendar.<br>
    !>                          For the current local date, this triple can be obtained from the Fortran intrinsic `date_and_time()` or [getDateTime()](@ref pm_dateTime::getDateTime).<br>
    !>                          Only the first three elements (`values(1:3)`) are used to compute the output.<br>
    !>                          The ability to pass longer vectors as input is to allow the output `values(1:8)` of various
    !>                          functionalities of this module to be passed directly to the procedures under this generic interface.<br>
    !>                          (**optional**. It **can** be present <b>if and only if</b> all other input arguments are missing.)
    !>  \param[in]  year    :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the year of the Gregorian calendar.<br>
    !>                          (**optional**. It **can** be present <b>if and only if</b> the input argument `values` is missing.)
    !>  \param[in]  month   :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the month of the Gregorian calendar.<br>
    !>                          (**optional**. It **must** be present <b>if and only if</b> the input argument `year` is present.)
    !>  \param[in]  day     :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the day of the Gregorian calendar.<br>
    !>                          (**optional**. It **must** be present <b>if and only if</b> the input argument `month` is present.)
    !>
    !>  \return
    !>  `dateAfter`         :   The output vector of size `3` of type `integer` of default kind \IK,
    !>                          containing the date in the Gregorian Calendar in the format `[year, month, day]` that appears after the input date.
    !>
    !>  \interface{getDateAfter}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getDateAfter
    !>      integer(IK) :: dateAfter(3)
    !>
    !>      dateAfter(1:3) = getDateAfter() ! current date is used.
    !>      dateAfter(1:3) = getDateAfter(values(1:3)) ! values = [year, month, day]
    !>      dateAfter(1:3) = getDateAfter(year, month, day)
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The size of the input argument `values(:)` must be at least `3` and at most `8`.<br>
    !>  The input values for the `year`, `month`, and `day` must be valid values.<br>
    !>  The input `month` must be a number between `1` and `12`.<br>
    !>  The input `day` must be a number between `1` and `31`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>  The procedures under this generic interface are non-elemental when the input argument `values(:)` is present.
    !>
    !>  \see
    !>  [isLastDayInMonth](@ref pm_dateTime::isLastDayInMonth)<br>
    !>
    !>  \example{getDateAfter}
    !>  \include{lineno} example/pm_dateTime/getDateAfter/main.F90
    !>  \compilef{getDateAfter}
    !>  \output{getDateAfter}
    !>  \include{lineno} example/pm_dateTime/getDateAfter/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getDateAfter}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface getDateAfter

    module function getDateAfterC() result(dateAfter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateAfterC
#endif
        integer(IK)                         :: dateAfter(3)
    end function

    PURE module function getDateAfterValues(values) result(dateAfter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateAfterValues
#endif
        integer(IK), intent(in), contiguous :: values(:)
        integer(IK)                         :: dateAfter(3)
    end function

    PURE module function getDateAfterTriple(year, month, day) result(dateAfter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateAfterTriple
#endif
        integer(IK), intent(in) :: year, month, day
        integer(IK)             :: dateAfter(3)
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the date in the Gregorian Calendar that appears before the input date.
    !>
    !>  \details
    !>  Returning the correct result requires taking into account the possibility of leap years.
    !>
    !>  \param[in]  values  :   The input `contiguous` array of shape `(:)`, of size `3` or larger, of type `integer` of default kind \IK, containing the `[year, month, day]` triple of the Gregorian calendar.<br>
    !>                          For the current local date, this triple can be obtained from the Fortran intrinsic `date_and_time()` or [getDateTime()](@ref pm_dateTime::getDateTime).<br>
    !>                          Only the first three elements (`values(1:3)`) are used to compute the output.<br>
    !>                          The ability to pass longer vectors as input is to allow the output `values(1:8)` of various
    !>                          functionalities of this module to be passed directly to the procedures under this generic interface.<br>
    !>                          (**optional**. It **can** be present <b>if and only if</b> all other input arguments are missing.)
    !>  \param[in]  year    :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the year of the Gregorian calendar.<br>
    !>                          (**optional**. It **can** be present <b>if and only if</b> the input argument `values` is missing.)
    !>  \param[in]  month   :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the month of the Gregorian calendar.<br>
    !>                          (**optional**. It **must** be present <b>if and only if</b> the input argument `year` is present.)
    !>  \param[in]  day     :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the day of the Gregorian calendar.<br>
    !>                          (**optional**. It **must** be present <b>if and only if</b> the input argument `month` is present.)
    !>
    !>  \return
    !>  `dateBefore`        :   The output vector of size `3` of type `integer` of default kind \IK,
    !>                          containing the date in the Gregorian Calendar in the format `[year, month, day]` that appears before the input date.
    !>
    !>  \interface{getDateBefore}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getDateBefore
    !>      integer(IK) :: dateBefore(3)
    !>
    !>      dateBefore(1:3) = getDateBefore() ! current date is used.
    !>      dateBefore(1:3) = getDateBefore(values(1:3)) ! values = [year, month, day]
    !>      dateBefore(1:3) = getDateBefore(year, month, day)
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The size of the input argument `values(:)` must be at least `3` and at most `8`.<br>
    !>  The input values for the `year`, `month`, and `day` must be valid values.<br>
    !>  The input `month` must be a number between `1` and `12`.<br>
    !>  The input `day` must be a number between `1` and `31`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>  The procedures under this generic interface are non-elemental when the input argument `values(:)` is present.
    !>
    !>  \see
    !>  [isLastDayInMonth](@ref pm_dateTime::isLastDayInMonth)<br>
    !>
    !>  \example{getDateBefore}
    !>  \include{lineno} example/pm_dateTime/getDateBefore/main.F90
    !>  \compilef{getDateBefore}
    !>  \output{getDateBefore}
    !>  \include{lineno} example/pm_dateTime/getDateBefore/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getDateBefore}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface getDateBefore

    module function getDateBeforeC() result(dateBefore)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateBeforeC
#endif
        integer(IK)                         :: dateBefore(3)
    end function

    PURE module function getDateBeforeValues(values) result(dateBefore)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateBeforeValues
#endif
        integer(IK), intent(in), contiguous :: values(:)
        integer(IK)                         :: dateBefore(3)
    end function

    PURE module function getDateBeforeTriple(year, month, day) result(dateBefore)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateBeforeTriple
#endif
        integer(IK), intent(in) :: year, month, day
        integer(IK)             :: dateBefore(3)
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **ordinal day**, also knowns as **Day of Year (DOY)**, i.e., the number of days
    !>  since the beginning of the input `year` until and including the input Gregorian Calendar date `[year, month, day]`.<br>
    !>
    !>  \details
    !>  The **day of year** is important for forming **ISO Ordinal Dates**.
    !>  An **Ordinal Date** is a simple form for occasions when the arbitrary nature of week and month definitions are more of an impediment than an aid,
    !>  for instance, when comparing dates from different calendars. The ordinal date is comprised of the year `[YYYY]` and the day of that year `[DDD]`,
    !>  from `001` through `365` (`366` in leap years). For example, `1981-04-05` corresponds to the ordinal date `1981-095`.
    !>
    !>  This format is used with simple hardware systems that have a need for a date system, but where including full calendar calculation software may be a significant nuisance.
    !>  This system is sometimes referred to as **Julian Date**, but this can cause confusion with the astronomical **Julian day**, a sequential count of the number of days since
    !>  day 0 beginning 1 January 4713 BC Greenwich noon, Julian proleptic calendar (or noon on ISO date `−4713-11-24` which uses the Gregorian proleptic calendar with a year `0000`).
    !>
    !>  Returning the correct result requires taking into account the possibility of leap years.
    !>
    !>  \param[in]  values  :   The input `contiguous` array of shape `(:)`, of size `3` or larger, of type `integer` of default kind \IK, containing the `[year, month, day]` triple of the Gregorian calendar.<br>
    !>                          For the current local date, this triple can be obtained from the Fortran intrinsic `date_and_time()` or [getDateTime()](@ref pm_dateTime::getDateTime).<br>
    !>                          Only the first three elements (`values(1:3)`) are used to compute the output.<br>
    !>                          The ability to pass longer vectors as input is to allow the output `values(1:8)` of various
    !>                          functionalities of this module to be passed directly to the procedures under this generic interface.<br>
    !>                          (**optional**. It **can** be present <b>if and only if</b> all other input arguments are missing.)
    !>  \param[in]  year    :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the year of the Gregorian calendar.<br>
    !>                          (**optional**. It **can** be present <b>if and only if</b> the input argument `values` is missing.)
    !>  \param[in]  month   :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the month of the Gregorian calendar.<br>
    !>                          (**optional**. It **must** be present <b>if and only if</b> the input argument `year` is present.)
    !>  \param[in]  day     :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the day of the Gregorian calendar.<br>
    !>                          (**optional**. It **must** be present <b>if and only if</b> the input argument `month` is present.)
    !>
    !>  \return
    !>  `ordinalDay`        :   The output scalar of type `integer` of default kind \IK, containing the **ordinal day** of the specified Gregorian Calendar date.<br>
    !>                          If all input arguments are missing, the ordinal day corresponding to the current Gregorian date is returned.
    !>
    !>  \interface{getOrdinalDay}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getOrdinalDay
    !>      integer(IK) :: ordinalDay
    !>
    !>      ordinalDay = getOrdinalDay() ! use the current date.
    !>      ordinalDay = getOrdinalDay(values(1:3)) ! values = [year, month, day]
    !>      ordinalDay = getOrdinalDay(year, month, day) ! elemental
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The size of the input argument `values(:)` must be at least `3` and at most `8`.<br>
    !>  The input values for the `year`, `month`, and `day` must be valid values.<br>
    !>  The input `month` must be a number between `1` and `12`.<br>
    !>  The input `day` must be a number between `1` and `31`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>  The procedures under this generic interface are always `impure` when all input arguments are missing.
    !>
    !>  \elemental
    !>  The procedures under this generic interface are non-elemental when all input arguments are missing or only the input argument `values(:)` is present.
    !>
    !>  \see
    !>  [isLastDayInMonth](@ref pm_dateTime::isLastDayInMonth)<br>
    !>
    !>  \example{getOrdinalDay}
    !>  \include{lineno} example/pm_dateTime/getOrdinalDay/main.F90
    !>  \compilef{getOrdinalDay}
    !>  \output{getOrdinalDay}
    !>  \include{lineno} example/pm_dateTime/getOrdinalDay/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getOrdinalDay}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface getOrdinalDay

    impure module function getOrdinalDayCurrent() result(ordinalDay)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOrdinalDayCurrent
#endif
        integer(IK)                         :: ordinalDay
    end function

    PURE module function getOrdinalDayValues(values) result(ordinalDay)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOrdinalDayValues
#endif
        integer(IK), intent(in), contiguous :: values(:)
        integer(IK)                         :: ordinalDay
    end function

    PURE elemental module function getOrdinalDayTriple(year, month, day) result(ordinalDay)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOrdinalDayTriple
#endif
        integer(IK), intent(in) :: year, month, day
        integer(IK)             :: ordinalDay
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **Week Number**, i.e., the number of weeks of the input `year` until
    !>  and **including** the week containing the specified input Gregorian date.
    !>
    !>  \details
    !>  The **Week Number** is important for constructing the **ISO week date system**, a *leap week calendar* system that is part of the
    !>  \f$\ms{ISO 8601}\f$ date and time standard issued by the International Organization for Standardization (**ISO**) since 1988
    !>  It is (mainly) used in government and business for fiscal years, as well as in timekeeping.<br>
    !>  The system specifies a week year atop the Gregorian calendar by defining a notation for **ordinal weeks** of the year.<br>
    !>
    !>  The \f$\ms{ISO 8601}\f$ definition for week number `01` is the week with the **first Thursday of January of the Gregorian year in it.<br>
    !>  The following definitions for the first week of the year are mutually equivalent, since **the ISO week starts with Monday**:
    !>  -#  It is the first week with a majority (`4` or more) of its days in `January`.
    !>  -#  Its first day is the Monday nearest to January 1st.
    !>  -#  It has January 4th in it. Hence,
    !>      -#  the **earliest possible first week extends from Monday December 29th of the previous Gregorian year to Sunday January 4th**,
    !>      -#  the latest possible first week extends from Monday 4 January to Sunday January 10th.
    !>  -#  It has the year's first working day in it, if Saturdays, Sundays and January 1st are not working days.
    !>  -#  **If January 1st is on a Monday, Tuesday, Wednesday or Thursday, it is in `W01`**.
    !>  -#  **If January 1st is on a Friday, it is part of `W53` of the previous year.**
    !>  -#  **If January 1st is on a Saturday, it is part of the last week of the previous year which is numbered `W52` in a common year and `W53` in a leap year.**
    !>  -#  **If it is on a Sunday, it is part of `W52` of the previous year.**
    !>
    !>  Returning the correct result requires taking into account the possibility of leap years.
    !>
    !>  \param[in]  values  :   The input `contiguous` array of shape `(:)`, of size `3` or larger, of type `integer` of default kind \IK, containing the `[year, month, day]` triple of the Gregorian calendar.<br>
    !>                          For the current local date, this triple can be obtained from the Fortran intrinsic `date_and_time()` or [getDateTime()](@ref pm_dateTime::getDateTime).<br>
    !>                          Only the first three elements (`values(1:3)`) are used to compute the output.<br>
    !>                          The ability to pass longer vectors as input is to allow the output `values(1:8)` of various
    !>                          functionalities of this module to be passed directly to the procedures under this generic interface.<br>
    !>                          (**optional**. It **can** be present <b>if and only if</b> all other input arguments are missing.)
    !>  \param[in]  year    :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the year of the Gregorian calendar.<br>
    !>                          (**optional**. It **can** be present <b>if and only if</b> the input argument `values` is missing.)
    !>  \param[in]  month   :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the month of the Gregorian calendar.<br>
    !>                          (**optional**. It **must** be present <b>if and only if</b> the input argument `year` is present.)
    !>  \param[in]  day     :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the day of the Gregorian calendar.<br>
    !>                          (**optional**. It **must** be present <b>if and only if</b> the input argument `month` is present.)
    !>
    !>  \return
    !>  `weekNumber`        :   The output scalar of type `integer` of default kind \IK, containing the **ordinal week** of the specified Gregorian Calendar date.<br>
    !>                          If all input arguments are missing, the ordinal week corresponding to the current Gregorian date is returned.
    !>
    !>  \interface{getWeekNumber}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getWeekNumber
    !>      integer(IK) :: weekNumber
    !>
    !>      weekNumber = getWeekNumber() ! use the current date.
    !>      weekNumber = getWeekNumber(values(1:3)) ! values = [year, month, day]
    !>      weekNumber = getWeekNumber(year, month, day) ! elemental
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The size of the input argument `values(:)` must be at least `3` and at most `8`.<br>
    !>  The input values for the `year`, `month`, and `day` must be valid values.<br>
    !>  The input `month` must be a number between `1` and `12`.<br>
    !>  The input `day` must be a number between `1` and `31`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>  The procedures under this generic interface are always `impure` when all input arguments are missing.
    !>
    !>  \elemental
    !>  The procedures under this generic interface are non-elemental when all input arguments are missing or only the input argument `values(:)` is present.
    !>
    !>  \see
    !>  [isLastDayInMonth](@ref pm_dateTime::isLastDayInMonth)<br>
    !>
    !>  \example{getWeekNumber}
    !>  \include{lineno} example/pm_dateTime/getWeekNumber/main.F90
    !>  \compilef{getWeekNumber}
    !>  \output{getWeekNumber}
    !>  \include{lineno} example/pm_dateTime/getWeekNumber/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getWeekNumber}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface getWeekNumber

    impure module function getWeekNumberCurrent() result(weekNumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getWeekNumberCurrent
#endif
        integer(IK)                         :: weekNumber
    end function

    PURE module function getWeekNumberValues(values) result(weekNumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getWeekNumberValues
#endif
        integer(IK), intent(in), contiguous :: values(:)
        integer(IK)                         :: weekNumber
    end function

    PURE elemental module function getWeekNumberTriple(year, month, day) result(weekNumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getWeekNumberTriple
#endif
        integer(IK), intent(in) :: year, month, day
        integer(IK)             :: weekNumber
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the day number of the week of a Gregorian Calendar date, **assuming Sunday is the zeroth day of the week**.
    !>  If there is no input argument, then the current day number of the week is returned.
    !>
    !>  \details
    !>  The internationally recognized method of conveying the day of the week is governed by ISO 8601,
    !>  which uses the [Zeller congruence](https://en.wikipedia.org/wiki/Zeller%27s_congruence) algorithm for
    !>  calculating the day of a week in a particular month and year, invented by Christian Zeller.<br>
    !>  <b>Monday is the official first day of the week according to ISO 8601.</b><br>
    !>
    !>  <b>To get an ISO-compliant day of week, see [getWeekDayISO](@ref pm_dateTime::getWeekDayISO).</b><br>
    !>
    !>  \param[in]  values  :   The input `contiguous` array of shape `(:)`, of size `3` or larger, of type `integer` of default kind \IK, containing the `[year, month, day]` triple of the Gregorian calendar.<br>
    !>                          For the current local date, this triple can be obtained from the Fortran intrinsic `date_and_time()` or [getDateTime()](@ref pm_dateTime::getDateTime).<br>
    !>                          Only the first three elements (`values(1:3)`) are used to compute the output.<br>
    !>                          The ability to pass longer vectors as input is to allow the output `values(1:8)` of various
    !>                          functionalities of this module to be passed directly to the procedures under this generic interface.<br>
    !>                          (**optional**. It **can** be present <b>if and only if</b> all other input arguments are missing.)
    !>  \param[in]  year    :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the year of the Gregorian calendar.<br>
    !>                          (**optional**. It **can** be present <b>if and only if</b> the input argument `values` is missing.)
    !>  \param[in]  month   :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the month of the Gregorian calendar.<br>
    !>                          (**optional**. It **must** be present <b>if and only if</b> the input argument `year` is present.)
    !>  \param[in]  day     :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the day of the Gregorian calendar.<br>
    !>                          (**optional**. It **must** be present <b>if and only if</b> the input argument `month` is present.)
    !>
    !>  \return
    !>  `weekday`           :   The output scalar or array of the same shape as the input array-like arguments (except `values`), of type `integer` of default kind \IK,
    !>                          containing the day of the week, starting with Sunday as day number `0`.
    !>
    !>  \interface{getWeekDay}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getWeekDay
    !>      integer(IK) :: weekday
    !>
    !>      weekday = getWeekDay() ! current date it used.
    !>      weekday = getWeekDay(values(1:3)) ! values = [year, month, day]
    !>      weekday = getWeekDay(year, month, day) ! elemental
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The size of the input argument `values(:)` must be at least `3` and at most `8`.<br>
    !>  The input values for the `year`, `month`, and `day` must be valid values.<br>
    !>  The input `month` must be a number between `1` and `12`.<br>
    !>  The input `day` must be a number between `1` and `31`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>  The procedures under this generic interface are non-elemental when the input argument `values(:)` is present.
    !>
    !>  \see
    !>  [getWeekDayISO](@ref pm_dateTime::getWeekDayISO)<br>
    !>  [WEEKDAY_NAME_ISO](@ref pm_dateTime::WEEKDAY_NAME_ISO)<br>
    !>  [WEEKDAY_NAME](@ref pm_dateTime::WEEKDAY_NAME)<br>
    !>  [strftime](https://www.cplusplus.com/reference/ctime/strftime/)<br>
    !>
    !>  \example{getWeekDay}
    !>  \include{lineno} example/pm_dateTime/getWeekDay/main.F90
    !>  \compilef{getWeekDay}
    !>  \output{getWeekDay}
    !>  \include{lineno} example/pm_dateTime/getWeekDay/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getWeekDay}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface getWeekDay

    impure module function getWeekDayCurrent() result(weekday)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getWeekDayCurrent
#endif
        integer(IK)                         :: weekday
    end function

    PURE module function getWeekDayValues(values) result(weekday)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getWeekDayValues
#endif
        integer(IK), intent(in), contiguous :: values(:)
        integer(IK)                         :: weekday
    end function

    PURE elemental module function getWeekDayTriple(year, month, day) result(weekday)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getWeekDayTriple
#endif
        integer(IK), intent(in) :: year, month, day
        integer(IK)             :: weekday
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the day number of the week of a Gregorian Calendar date, **assuming Monday is the first day of the week**.
    !>  If there is no input argument, then the current day number of the week is returned.
    !>
    !>  \details
    !>  The internationally recognized method of conveying the day of the week is governed by ISO 8601,
    !>  which uses the [Zeller congruence](https://en.wikipedia.org/wiki/Zeller%27s_congruence) algorithm for
    !>  calculating the day of a week in a particular month and year, invented by Christian Zeller.<br>
    !>  **Monday is the official first day of the week according to ISO 8601.**<br><br>
    !>  **Fun Facts:**<br>
    !>  -#  In the Gregorian Calendar, over a period of four hundred years, there are `97` leap years and `303` normal years.
    !>  -#  In each normal year, the weekday of January `1` advances by `1`. For each leap year it advances by `2`.
    !>  -#  Consequently, January `1` of year `N` occurs on the same day of the week as January `1` of year `N + 400`.
    !>  -#  Because the leap year pattern also recurs with a four hundred year cycle, a simple table of four hundred elements,
    !>      and single modulus, suffices to determine the day of the week (in the Gregorian Calendar).
    !>  -#  Because `7` does not divide `400`, January `1` occurs more frequently on some days than others.
    !>  -#  Similarly, in the Mathematical Gazette, vol. 53, pp.127-129, it is shown that the 13th of the month is more likely
    !>      to be a Friday than any other day.
    !>
    !>  \param[in]  values  :   The input `contiguous` array of shape `(:)`, of size `3` or larger, of type `integer` of default kind \IK, containing the `[year, month, day]` triple of the Gregorian calendar.<br>
    !>                          For the current local date, this triple can be obtained from the Fortran intrinsic `date_and_time()` or [getDateTime()](@ref pm_dateTime::getDateTime).<br>
    !>                          Only the first three elements (`values(1:3)`) are used to compute the output.<br>
    !>                          The ability to pass longer vectors as input is to allow the output `values(1:8)` of various
    !>                          functionalities of this module to be passed directly to the procedures under this generic interface.<br>
    !>                          (**optional**. It **can** be present <b>if and only if</b> all other input arguments are missing.)
    !>  \param[in]  year    :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the year of the Gregorian calendar.<br>
    !>                          (**optional**. It **can** be present <b>if and only if</b> the input argument `values` is missing.)
    !>  \param[in]  month   :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the month of the Gregorian calendar.<br>
    !>                          (**optional**. It **must** be present <b>if and only if</b> the input argument `year` is present.)
    !>  \param[in]  day     :   The input scalar, or array of the same shape as other array-like arguments, of type `integer` of default kind \IK, containing the day of the Gregorian calendar.<br>
    !>                          (**optional**. It **must** be present <b>if and only if</b> the input argument `month` is present.)
    !>
    !>  \return
    !>  `weekday`           :   The output scalar or array of the same shape as the input array-like arguments of type `integer` of default kind \IK,
    !>                          containing the day of the week, starting with Monday as day number `1`.
    !>
    !>  \interface{getWeekDayISO}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getWeekDayISO
    !>      integer(IK) :: weekday
    !>
    !>      weekday = getWeekDayISO() ! current date is used.
    !>      weekday = getWeekDayISO(values(1:3)) ! values = [year, month, day]
    !>      weekday = getWeekDayISO(year, month, day) ! elemental
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The size of the input argument `values(:)` must be at least `3` and at most `8`.<br>
    !>  The input values for the `year`, `month`, and `day` must be valid values.<br>
    !>  The input `month` must be a number between `1` and `12`.<br>
    !>  The input `day` must be a number between `1` and `31`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>  The procedures under this generic interface are non-elemental when the input argument `values(:)` is present.
    !>
    !>  \see
    !>  [getWeekDay](@ref pm_dateTime::getWeekDay)<br>
    !>  [WEEKDAY_NAME](@ref pm_dateTime::WEEKDAY_NAME)<br>
    !>  [WEEKDAY_NAME_ISO](@ref pm_dateTime::WEEKDAY_NAME_ISO)<br>
    !>
    !>  \example{getWeekDayISO}
    !>  \include{lineno} example/pm_dateTime/getWeekDayISO/main.F90
    !>  \compilef{getWeekDayISO}
    !>  \output{getWeekDayISO}
    !>  \include{lineno} example/pm_dateTime/getWeekDayISO/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getWeekDayISO}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface getWeekDayISO

    impure module function getWeekDayISOCurrent() result(weekday)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getWeekDayISOCurrent
#endif
        integer(IK)                         :: weekday
    end function

    PURE module function getWeekDayISOValues(values) result(weekday)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getWeekDayISOValues
#endif
        integer(IK), intent(in), contiguous :: values(:)
        integer(IK)                         :: weekday
    end function

    PURE elemental module function getWeekDayISOTriple(year, month, day) result(weekday)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getWeekDayISOTriple
#endif
        integer(IK), intent(in) :: year, month, day
        integer(IK)             :: weekday
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the total number of days in the specified year or the month of the year of the Gregorian calendar.
    !>
    !>  \details
    !>  The correct computation of the number of days in a year or a number of years requires taking into account the leap years.
    !>
    !>  \param[in]  year    :   The input scalar or array of the same shape as other array-like input arguments,
    !>                          of type `integer` of default kind \IK, representing the year of the Gregorian Calendar.
    !>  \param[in]  month   :   The input scalar or array of the same shape as other array-like input arguments,
    !>                          of type `integer` of default kind \IK, representing the month of the year of the Gregorian Calendar.<br>
    !>                          (**optional**, if present, the number of days within the specified `month` of the `year` will be returned.)
    !>
    !>  \return
    !>  `countDays`         :   The output scalar or array of the same shape as other array-like arguments of type `integer` of default kind \IK,
    !>                          containing the number of days within the specified `year` (or within the specified `month` of the `year`).
    !>
    !>  \interface{getCountDays}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getCountDays
    !>      use pm_kind, only: IK
    !>      integer(IK) :: countDaysInYear
    !>      integer(IK) :: countDaysInMonthOfYear
    !>
    !>      countDaysInYear = getCountDays(year)
    !>      countDaysInMonthOfYear = getCountDays(year, month)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  An input argument `year = 0` corresponds to the historic 1 BC notation of the Gregorian calendar.<br>
    !>  This is in accordance with the convention in astronomical year numbering and the international standard date system, **ISO 8601**.<br>
    !>  The input argument `month` must be a number between `1` and `12`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [dateTimeInt_typer](@ref pm_dateTime::dateTimeInt_typer) (class constructor)<br>
    !>  [dateTimeStr_type](@ref pm_dateTime::dateTimeStr_type)<br>
    !>
    !>  \example{getCountDays}
    !>  \include{lineno} example/pm_dateTime/getCountDays/main.F90
    !>  \compilef{getCountDays}
    !>  \output{getCountDays}
    !>  \include{lineno} example/pm_dateTime/getCountDays/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getCountDays}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface getCountDays

    PURE elemental module function getCountDaysInYear(year) result(countDays)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountDaysInYear
#endif
        integer(IK), intent(in) :: year
        integer(IK)             :: countDays
    end function

    PURE elemental module function getCountDaysInMonth(year, month) result(countDays)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountDaysInMonth
#endif
        integer(IK), intent(in) :: year, month
        integer(IK)             :: countDays
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the number of weeks in the specified year or the month of the year of the Gregorian calendar
    !>  **whose majority of days** fall within the specified year or month.
    !>
    !>  \details
    !>  The correct computation of the number of weeks in a year or a number of years requires taking into account the leap years.<br>
    !>  Note the number of weeks with majority of days in a given year can be only `53` (for **long years**) or `52` (for **short years**).
    !>
    !>  \param[in]  year    :   The input scalar or array of the same shape as other array-like input arguments,
    !>                          of type `integer` of default kind \IK, representing the year of the Gregorian Calendar.
    !>  \param[in]  month   :   The input scalar or array of the same shape as other array-like input arguments,
    !>                          of type `integer` of default kind \IK, representing the month of the year of the Gregorian Calendar.<br>
    !>                          (**optional**, if present, the number of major weeks within the specified month of the year will be returned.)
    !>
    !>  \return
    !>  `countWeeks`        :   The output scalar or array of the same shape as other array-like arguments of type `integer` of default kind \IK,
    !>                          containing the number of weeks with majority (4 or more) of their days within the specified `year`
    !>                          or within the specified `month` of the `year`.
    !>
    !>  \interface{getCountWeeks}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getCountWeeks
    !>      use pm_kind, only: IK
    !>      integer(IK) :: countWeeksInYear
    !>      integer(IK) :: countWeeksInMonthOfYear
    !>
    !>      countWeeksInYear = getCountWeeks(year)
    !>      countWeeksInMonthOfYear = getCountWeeks(year, month)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  An input argument `year = 0` corresponds to the historic 1 BC notation of the Gregorian calendar.<br>
    !>  This is in accordance with the convention in astronomical year numbering and the international standard date system, **ISO 8601**.<br>
    !>  The input argument `month` must be a number between `1` and `12`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [dateTimeInt_typer](@ref pm_dateTime::dateTimeInt_typer) (class constructor)<br>
    !>  [dateTimeStr_type](@ref pm_dateTime::dateTimeStr_type)<br>
    !>
    !>  \example{getCountWeeks}
    !>  \include{lineno} example/pm_dateTime/getCountWeeks/main.F90
    !>  \compilef{getCountWeeks}
    !>  \output{getCountWeeks}
    !>  \include{lineno} example/pm_dateTime/getCountWeeks/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getCountWeeks}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface getCountWeeks

    PURE elemental module function getCountWeeksInYear(year) result(countWeeks)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountWeeksInYear
#endif
        integer(IK), intent(in) :: year
        integer(IK)             :: countWeeks
    end function

    PURE elemental module function getCountWeeksInMonth(year, month) result(countWeeks)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountWeeksInMonth
#endif
        integer(IK), intent(in) :: year, month
        integer(IK)             :: countWeeks
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the number of leap years within the closed year interval `[since, until]`
    !>  where `since` represents the origin year, and `until` represents a later year after `since`.
    !>
    !>  \param[in]  until   :   The input scalar or array of arbitrary rank of type `integer` of default kind \IK,
    !>                          containing the Gregorian Calendar year up to which the number of leap years since `since` must be returned.<br>
    !>                          If `until < 0`, then the number of leap years from the negative year until **January, 1, 1** is returned.
    !>  \param[in]  since   :   The input scalar or array of the same shape as `until` of the same type and kind as `until`,
    !>                          containing the Gregorian Calendar year that marks the beginning of the period within which the number of leap years must be counted.<br>
    !>                          (**optional**, default = `1_IK`, i.e., the origin of the Gregorian Calendar.)
    !>
    !>  \interface{getCountLeapYears}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getCountLeapYears
    !>      use pm_kind, only: IK
    !>      integer(IK) :: countLeapYear
    !>
    !>      countLeapYear = getCountLeapYears(until)
    !>      countLeapYear = getCountLeapYears(until, since)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < since < until` must hold for the input argument values.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  The interface of this generic interface is intentionally defined to have `until` as the first input argument (as opposed to `since`).<br>
    !>  This is to ensure the purity of the procedures under this generic interface.<br>
    !>  Making `until` optional would require a call to the `impure` Fortran intrinsic `date_and_time()`.<br>
    !>
    !>  \see
    !>  [getCountDays](@ref pm_dateTime::getCountDays)<br>
    !>
    !>  \example{getCountLeapYears}
    !>  \include{lineno} example/pm_dateTime/getCountLeapYears/main.F90
    !>  \compilef{getCountLeapYears}
    !>  \output{getCountLeapYears}
    !>  \include{lineno} example/pm_dateTime/getCountLeapYears/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getCountLeapYears}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface getCountLeapYears

    PURE elemental module function getCountLeapYears(until) result(countLeapYear)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLeapYears
#endif
        integer(IK), intent(in) :: until
        integer(IK)             :: countLeapYear
    end function

    PURE elemental module function getCountLeapYearsSince(until, since) result(countLeapYear)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLeapYearsSince
#endif
        integer(IK), intent(in) :: until
        integer(IK), intent(in) :: since
        integer(IK)             :: countLeapYear
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the hour in the current or the input time zone `(-720 : +840)` is **Ante Meridiem** (before noon).
    !>
    !>  \param[in]  zone        :   The input scalar or array of arbitrary shape of type `integer`
    !>                              of default kind \IK containing the time zone of interest in units of minutes.<br>
    !>                              (**optional**, default = the current local time zone.)
    !>  \param[in]  julianDay   :   The input scalar or array of arbitrary shape of type `real` of default
    !>                              kind \RK containing the (possibly proleptic) Julian Day.<br>
    !>                              (**optional**, default = the current locale Julian Day.)
    !>
    !>  \return
    !>  `isMorning`             :   The output scalar or array of the same rank as input array-like arguments, of type `logical` of default kind \LK.
    !>                              It is `.true.` <b>if and only if</b> the current time or the time corresponding to the input `zone` or `julianDay` or both represents represents morning.
    !>
    !>  \interface{isMorning}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: isMorning
    !>      real(RK)    :: julianDay
    !>      logical(LK) :: morning
    !>      integer(IK) :: zone
    !>
    !>      morning = isMorning() ! impure : current local hour
    !>      morning = isMorning(zone) ! elemental : current (realtime) hour in the specified time zone
    !>      morning = isMorning(julianDay) ! elemental
    !>      morning = isMorning(julianDay, zone) ! elemental
    !>      !
    !>  \endcode
    !>
    !>  \warnpure
    !>  The procedures under this generic interface are always `impure` when there is no input argument.
    !>
    !>  \warning
    !>  The conditions `-12 * 60 < zone` and `zone < +14 * 60` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  The procedures under this generic interface are non-elemental when there is no input argument.
    !>
    !>  \see
    !>  [getMillisecond](@ref pm_dateTime::getMillisecond)<br>
    !>  [getSecond](@ref pm_dateTime::getSecond)<br>
    !>  [getMinute](@ref pm_dateTime::getMinute)<br>
    !>  [getHour12](@ref pm_dateTime::getHour12)<br>
    !>  [getHour](@ref pm_dateTime::getHour)<br>
    !>  [getZone](@ref pm_dateTime::getZone)<br>
    !>  [getZoneAbbr](@ref pm_dateTime::getZoneAbbr)<br>
    !>  [getDay](@ref pm_dateTime::getDay)<br>
    !>  [getMonth](@ref pm_dateTime::getMonth)<br>
    !>  [getYear](@ref pm_dateTime::getYear)<br>
    !>  [isLeapYear](@ref pm_dateTime::isLeapYear)<br>
    !>
    !>  \example{isMorning}
    !>  \include{lineno} example/pm_dateTime/isMorning/main.F90
    !>  \compilef{isMorning}
    !>  \output{isMorning}
    !>  \include{lineno} example/pm_dateTime/isMorning/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{isMorning}
    !>
    !>  \author
    !>  \AmirShahmoradi, January 30, 2021, 5:36 AM, Dallas, TX
    interface isMorning

    impure module function isMorningCurrent() result(morning)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMorningCurrent
#endif
        use pm_kind, only: IKG => IK, LKG => LK
        logical(LKG)                            :: morning
    end function

    impure elemental module function isMorningZ(zone) result(morning)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMorningZ
#endif
        use pm_kind, only: IKG => IK, LKG => LK
        integer(IKG)    , intent(in)            :: zone
        logical(LKG)                            :: morning
    end function

    pure elemental module function isMorningJD(julianDay) result(morning)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMorningJD
#endif
        use pm_kind, only: IKG => IK, LKG => LK, RKG => RK
        real(RKG)       , intent(in)            :: julianDay
        logical(LKG)                            :: morning
    end function

    PURE elemental module function isMorningJDZ(julianDay, zone) result(morning)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMorningJDZ
#endif
        use pm_kind, only: IKG => IK, LKG => LK, RKG => RK
        real(RKG)       , intent(in)            :: julianDay
        integer(IKG)    , intent(in)            :: zone
        logical(LKG)                            :: morning
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` is the input time zone `zone` in **units of minutes** is valid (i.e., `-12 * 60 <= zone <= +14 * 60`).<br>
    !>
    !>  \param[in]  zone    :   The output scalar, or array of arbitrary shape, rank, and size, of type `integer` of default kind \IK, containing
    !>                          the time zone, that is, the time difference **in minutes** with respect to the Coordinated Universal Time (UTC).
    !>
    !>  \return
    !>  `isValid`           :   The output scalar or array of the same shape, rank, and size as the input `zone`, of type `logical` of default kind \LK,
    !>                          each element of which is `.true.` only if the condition `-12 * 60 <= zone <= +14 * 60` holds for the corresponding input `zone`.<br>
    !>
    !>  \interface{isValidZone}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: isValidZone
    !>      use pm_kind, only: IK, LK
    !>      logical(LK) :: isValid
    !>      integer(IK) :: zone
    !>
    !>      isValid = isValidZone(zone)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  This generic interface is meant to be primarily used for internal testing purposes of the ParaMonte library.<br>
    !>  The time zone settings are highly fluid and flexible across the world.<br>
    !>  The current minimum and maximum time zones can change in future.<br>
    !>
    !>  \see
    !>  [getMillisecond](@ref pm_dateTime::getMillisecond)<br>
    !>  [getSecond](@ref pm_dateTime::getSecond)<br>
    !>  [getMinute](@ref pm_dateTime::getMinute)<br>
    !>  [getHour12](@ref pm_dateTime::getHour12)<br>
    !>  [getHour](@ref pm_dateTime::getHour)<br>
    !>  [getZone](@ref pm_dateTime::getZone)<br>
    !>  [getZoneAbbr](@ref pm_dateTime::getZoneAbbr)<br>
    !>  [getDay](@ref pm_dateTime::getDay)<br>
    !>  [getMonth](@ref pm_dateTime::getMonth)<br>
    !>  [getYear](@ref pm_dateTime::getYear)<br>
    !>  [isLeapYear](@ref pm_dateTime::isLeapYear)<br>
    !>
    !>  \example{isValidZone}
    !>  \include{lineno} example/pm_dateTime/isValidZone/main.F90
    !>  \compilef{isValidZone}
    !>  \output{isValidZone}
    !>  \include{lineno} example/pm_dateTime/isValidZone/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{isValidZone}
    !>
    !>  \author
    !>  \AmirShahmoradi, January 30, 2021, 5:36 AM, Dallas, TX
    interface isValidZone

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure elemental module function isValidZone_IK(zone) result(isValid)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isValidZone_IK
#endif
        use pm_kind, only: LKG => LK, IKG => IK
        integer(IKG), intent(in)    :: zone
        logical(LKG)                :: isValid
    end function

!#if IK5_ENABLED
!    pure elemental module function isValidZone_IK5(zone) result(isValid)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: isValidZone_IK5
!#endif
!        use pm_kind, only: LKG => LK, IKG => IK5
!        integer(IKG), intent(in)    :: zone
!        logical(LKG)                :: isValid
!    end function
!#endif
!
!#if IK4_ENABLED
!    pure elemental module function isValidZone_IK4(zone) result(isValid)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: isValidZone_IK4
!#endif
!        use pm_kind, only: LKG => LK, IKG => IK4
!        integer(IKG), intent(in)    :: zone
!        logical(LKG)                :: isValid
!    end function
!#endif
!
!#if IK3_ENABLED
!    pure elemental module function isValidZone_IK3(zone) result(isValid)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: isValidZone_IK3
!#endif
!        use pm_kind, only: LKG => LK, IKG => IK3
!        integer(IKG), intent(in)    :: zone
!        logical(LKG)                :: isValid
!    end function
!#endif
!
!#if IK2_ENABLED
!    pure elemental module function isValidZone_IK2(zone) result(isValid)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: isValidZone_IK2
!#endif
!        use pm_kind, only: LKG => LK, IKG => IK2
!        integer(IKG), intent(in)    :: zone
!        logical(LKG)                :: isValid
!    end function
!#endif
!
!#if IK1_ENABLED
!    pure elemental module function isValidZone_IK1(zone) result(isValid)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: isValidZone_IK1
!#endif
!        use pm_kind, only: LKG => LK, IKG => IK1
!        integer(IKG), intent(in)    :: zone
!        logical(LKG)                :: isValid
!    end function
!#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

#if     CHECK_ENABLED
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
block; \
use pm_val2str, only: getStr; \
use pm_err, only: setAsserted, getFine; \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG); \
end block;
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the constructor of the [dateTimeInt_type](@ref pm_dateTime::dateTimeInt_type) class.<br>
    !>
    !>  \details
    !>  Upon return, the constructor initializes all components of the object to the current date and time.<br>
    !>  See also the documentation details of [dateTimeInt_type](@ref pm_dateTime::dateTimeInt_type).
    !>
    !>  \return
    !>  `dateTimeInt`   :   The output scalar object of class [dateTimeInt_type](@ref pm_dateTime::dateTimeInt_type).
    !>
    !>  \interface{dateTimeInt_typer}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: dateTimeInt_type
    !>      type(dateTimeInt_type) :: dateTimeInt
    !>
    !>      dateTimeInt = dateTimeInt_type()
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  All object components are set to `-huge(0_K)` if the processor does not provide date and time.<br>
    !>  \vericon
    !>
    !>  \remark
    !>  See the documentation of [dateTimeInt_type](@ref pm_dateTime::dateTimeInt_type) for example usage.
    !>
    !>  \final{dateTimeInt_typer}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    function dateTimeInt_typer() result(dateTimeInt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: dateTimeInt_typer
#endif
        type(dateTimeInt_type) :: dateTimeInt
        integer(IK) :: values(8)
        call date_and_time(values = values)
        CHECK_ASSERTION(__LINE__, all(values /= -huge(values)), SK_"@dateTimeInt_typer(): The processor does not provide date and time.") ! fpp
        dateTimeInt%year        = values(1)
        dateTimeInt%month       = values(2)
        dateTimeInt%day         = values(3)
        dateTimeInt%zone        = values(4)
        dateTimeInt%hour        = values(5)
        dateTimeInt%minute      = values(6)
        dateTimeInt%second      = values(7)
        dateTimeInt%millisecond = values(8)
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return an `integer` vector of length `8` containing all component values
    !>  of the parent object of type [dateTimeInt_type](@ref pm_dateTime::dateTimeInt_type).
    !>
    !>  \details
    !>  This is a dynamic method of the [dateTimeInt_type](@ref pm_dateTime::dateTimeInt_type) class.<br>
    !>  See also the documentation details of [dateTimeInt_type](@ref pm_dateTime::dateTimeInt_type).
    !>
    !>  \return
    !>  `values`    :   The output vector of type `integer` of default kind \IK containing the values of all components
    !>                  of the parent object of class [dateTimeInt_type](@ref pm_dateTime::dateTimeInt_type) as
    !>                  `[year, month, day, zone, hour, minute, millisecond]`.<br>
    !>                  The format of the output vector conforms with the format of the output `values`
    !>                  argument of the Fortran intrinsic `date_and_time()`.
    !>
    !>  \interface{getDateTimeIntValues}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: dateTimeInt_type
    !>      type(dateTimeInt_type) :: dateTimeInt
    !>
    !>      dateTimeInt = dateTimeInt_type()
    !>
    !>  \endcode
    !>
    !>  \remark
    !>  See the documentation of [dateTimeInt_type](@ref pm_dateTime::dateTimeInt_type) for example usage.
    !>
    !>  \final{getDateTimeIntValues}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    function getDateTimeIntValues(self) result(values)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDateTimeIntValues
#endif
        class(dateTimeInt_type), intent(in) :: self
        integer(IK) :: values(8)
        values(1) = self%year
        values(2) = self%month
        values(3) = self%day
        values(4) = self%zone
        values(5) = self%hour
        values(6) = self%minute
        values(7) = self%second
        values(8) = self%millisecond
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the constructor of the [dateTimeStr_type](@ref pm_dateTime::dateTimeStr_type) class.<br>
    !>
    !>  \details
    !>  Upon return, the constructor initializes all components of the object to the current date and time.<br>
    !>  See also the documentation details of [dateTimeStr_type](@ref pm_dateTime::dateTimeStr_type).
    !>
    !>  \return
    !>  `dateTimeStr`   :   The output scalar object of class [dateTimeStr_type](@ref pm_dateTime::dateTimeStr_type).
    !>
    !>  \interface{dateTimeStr_typer}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: dateTimeStr_type
    !>      type(dateTimeStr_type) :: dateTimeStr
    !>
    !>      dateTimeStr = dateTimeStr_type()
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  All object components are set to blanks if the processor does not provide date and time.<br>
    !>  \vericon
    !>
    !>  \remark
    !>  See the documentation of [dateTimeStr_type](@ref pm_dateTime::dateTimeStr_type) for example usage.
    !>
    !>  \final{dateTimeStr_typer}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    function dateTimeStr_typer() result(dateTimeStr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: dateTimeStr_typer
#endif
        type(dateTimeStr_type)  :: dateTimeStr
        character(10)           :: time
        character(8)            :: date
        character(5)            :: zone
        call date_and_time(date = date, time = time, zone = zone)
        CHECK_ASSERTION(__LINE__, len_trim(date) /= 0 .or. len_trim(time) /= 0 .or. len_trim(zone) /= 0, SK_"@dateTimeStr_typer(): The processor does not provide date and time.") ! fpp
       !dateTimeStr%century            = date(1:2)
        dateTimeStr%year               = date(1:4)
        dateTimeStr%month              = date(5:6)
        dateTimeStr%day                = date(7:8)
        dateTimeStr%zone               = zone
        dateTimeStr%hour               = time(1:2)
        dateTimeStr%minute             = time(3:4)
        dateTimeStr%second             = time(5:6)
        dateTimeStr%millisecond        = time(8:10)
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input year is a **leap year**.
    !>
    !>  \details
    !>  To be a leap year, the year number must be divisible by four – except for end-of-century years, which must be divisible by 400.<br>
    !>  For BC, counting starts at `1`, so there is no year `0`. This means that the leap years are offset by `1` and
    !>  can be calculated by the same method as above, but with the year number increased by 1.<br>
    !>  For example, the year 2000 was a leap year, although 1900 was not.<br>
    !>  The years 2020, 2024 and 2028 are all leap years.
    !>
    !>  \param[in]  year    :   The input scalar or array of arbitrary shape of type `integer` of default kind \IK
    !>                          containing the (possibly proleptic) Gregorian calendar year(s).
    !>
    !>  \return
    !>  `leapYear`          :   The output object of the same rank and shape as the input `year`, of type `logical` of default kind \LK.
    !>                          It is `.true.` <b>if and only if</b> the corresponding input year is a leap year.
    !>
    !>  \interface{isLeapYear}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: isLeapYear
    !>      logical(LK) :: leapYear
    !>      integer(IK) :: year
    !>
    !>      leapYear = isLeapYear(year)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \warning
    !>  An input argument `year = 0` corresponds to the historic 1 BC notation of the Gregorian calendar.<br>
    !>  This is in accordance with the convention in astronomical year numbering and the international standard date system, **ISO 8601**.<br>
    !>  In these systems, the year `0` is a leap year.
    !>
    !>  \see
    !>  [getMillisecond](@ref pm_dateTime::getMillisecond)<br>
    !>  [getSecond](@ref pm_dateTime::getSecond)<br>
    !>  [getMinute](@ref pm_dateTime::getMinute)<br>
    !>  [getHour12](@ref pm_dateTime::getHour12)<br>
    !>  [getHour](@ref pm_dateTime::getHour)<br>
    !>  [getZone](@ref pm_dateTime::getZone)<br>
    !>  [getZoneAbbr](@ref pm_dateTime::getZoneAbbr)<br>
    !>  [getDay](@ref pm_dateTime::getDay)<br>
    !>  [getMonth](@ref pm_dateTime::getMonth)<br>
    !>  [getYear](@ref pm_dateTime::getYear)<br>
    !>  [isLeapYear](@ref pm_dateTime::isLeapYear)<br>
    !>
    !>  \example{isLeapYear}
    !>  \include{lineno} example/pm_dateTime/isLeapYear/main.F90
    !>  \compilef{isLeapYear}
    !>  \output{isLeapYear}
    !>  \include{lineno} example/pm_dateTime/isLeapYear/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{isLeapYear}
    !>
    !>  \author
    !>  \AmirShahmoradi, January 30, 2021, 5:36 AM, Dallas, TX
    pure elemental function isLeapYear(year) result(leapYear)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isLeapYear
#endif
        integer(IK), intent(in) :: year
        logical(LK)             :: leapYear
        !check_assertion(__LINE__, year /= 0_IK, SK_"@isLeapYear(): The input year must be non-zero.") ! fpp
        leapYear = mod(year, 400_IK) == 0_IK .or. (mod(year, 4_IK) == 0_IK .and. mod(year, 100_IK) /= 0_IK)
        !if (year > 0_IK) then
        !    leapYear = mod(year, 400_IK) == 0_IK .or. (mod(year, 4_IK) == 0_IK .and. mod(year, 100_IK) /= 0_IK)
        !else
        !    leapYear = mod(year + 1_IK, 400_IK) == 0_IK .or. (mod(year + 1_IK, 4_IK) == 0_IK .and. mod(year + 1_IK, 100_IK) /= 0_IK)
        !end if
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the current millisecond of the current second of the current minute of
    !>  the local hour of the current day of the Gregorian calendar.
    !>
    !>  \return
    !>  `second`    :   The output scalar of type `integer` of default kind \IK containing
    !>                  the current millisecond of the current second of the current minute
    !>                  of the local hour of the current day of the Gregorian calendar.
    !>
    !>  \interface{getMillisecond}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getMillisecond
    !>      integer(IK) :: millisecond
    !>
    !>      millisecond = getSecond()
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The output value is set to `-huge(0_K)` if the processor does not provide date and time.<br>
    !>  \vericon
    !>
    !>  \see
    !>  [getMillisecond](@ref pm_dateTime::getMillisecond)<br>
    !>  [getSecond](@ref pm_dateTime::getSecond)<br>
    !>  [getMinute](@ref pm_dateTime::getMinute)<br>
    !>  [getHour12](@ref pm_dateTime::getHour12)<br>
    !>  [getHour](@ref pm_dateTime::getHour)<br>
    !>  [getZone](@ref pm_dateTime::getZone)<br>
    !>  [getZoneAbbr](@ref pm_dateTime::getZoneAbbr)<br>
    !>  [getDay](@ref pm_dateTime::getDay)<br>
    !>  [getMonth](@ref pm_dateTime::getMonth)<br>
    !>  [getYear](@ref pm_dateTime::getYear)<br>
    !>  [isLeapYear](@ref pm_dateTime::isLeapYear)<br>
    !>
    !>  \example{getMillisecond}
    !>  \include{lineno} example/pm_dateTime/getMillisecond/main.F90
    !>  \compilef{getMillisecond}
    !>  \output{getMillisecond}
    !>  \include{lineno} example/pm_dateTime/getMillisecond/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getMillisecond}
    !>
    !>  \author
    !>  \AmirShahmoradi, January 30, 2021, 5:36 AM, Dallas, TX
    function getMillisecond() result(millisecond)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMillisecond
#endif
        integer(IK) :: millisecond
        integer(IK) :: values(8)
        call date_and_time(values = values)
        CHECK_ASSERTION(__LINE__, all(values /= -huge(values)), SK_"@getMillisecond(): The processor does not provide date and time.") ! fpp
        millisecond = values(8)
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the current second of the current minute of the local hour of the current day of the Gregorian calendar.
    !>
    !>  \return
    !>  `second`    :   The output scalar of type `integer` of default kind \IK containing
    !>                  the current second of the current minute of the local hour of the current day of the Gregorian calendar.
    !>
    !>  \interface{getSecond}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getSecond
    !>      integer(IK) :: second
    !>
    !>      second = getSecond()
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The output value is set to `-huge(0_K)` if the processor does not provide date and time.<br>
    !>  \vericon
    !>
    !>  \see
    !>  [getMillisecond](@ref pm_dateTime::getMillisecond)<br>
    !>  [getSecond](@ref pm_dateTime::getSecond)<br>
    !>  [getMinute](@ref pm_dateTime::getMinute)<br>
    !>  [getHour12](@ref pm_dateTime::getHour12)<br>
    !>  [getHour](@ref pm_dateTime::getHour)<br>
    !>  [getZone](@ref pm_dateTime::getZone)<br>
    !>  [getZoneAbbr](@ref pm_dateTime::getZoneAbbr)<br>
    !>  [getDay](@ref pm_dateTime::getDay)<br>
    !>  [getMonth](@ref pm_dateTime::getMonth)<br>
    !>  [getYear](@ref pm_dateTime::getYear)<br>
    !>  [isLeapYear](@ref pm_dateTime::isLeapYear)<br>
    !>
    !>  \example{getSecond}
    !>  \include{lineno} example/pm_dateTime/getSecond/main.F90
    !>  \compilef{getSecond}
    !>  \output{getSecond}
    !>  \include{lineno} example/pm_dateTime/getSecond/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getSecond}
    !>
    !>  \author
    !>  \AmirShahmoradi, January 30, 2021, 5:36 AM, Dallas, TX
    function getSecond() result(second)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSecond
#endif
        integer(IK) :: second
        integer(IK) :: values(8)
        call date_and_time(values = values)
        CHECK_ASSERTION(__LINE__, all(values /= -huge(values)), SK_"@getSecond(): The processor does not provide date and time.") ! fpp
        second = values(7)
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the current minute of the local hour of the current day of the Gregorian calendar.
    !>
    !>  \return
    !>  `minute`    :   The output scalar of type `integer` of default kind \IK containing
    !>                  the current minute of the local hour of the current day of the Gregorian calendar.
    !>
    !>  \interface{getMinute}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getMinute
    !>      integer(IK) :: minute
    !>
    !>      minute = getMinute()
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The output value is set to `-huge(0_K)` if the processor does not provide date and time.<br>
    !>  \vericon
    !>
    !>  \see
    !>  [getMillisecond](@ref pm_dateTime::getMillisecond)<br>
    !>  [getSecond](@ref pm_dateTime::getSecond)<br>
    !>  [getMinute](@ref pm_dateTime::getMinute)<br>
    !>  [getHour12](@ref pm_dateTime::getHour12)<br>
    !>  [getHour](@ref pm_dateTime::getHour)<br>
    !>  [getZone](@ref pm_dateTime::getZone)<br>
    !>  [getZoneAbbr](@ref pm_dateTime::getZoneAbbr)<br>
    !>  [getDay](@ref pm_dateTime::getDay)<br>
    !>  [getMonth](@ref pm_dateTime::getMonth)<br>
    !>  [getYear](@ref pm_dateTime::getYear)<br>
    !>  [isLeapYear](@ref pm_dateTime::isLeapYear)<br>
    !>
    !>  \example{getMinute}
    !>  \include{lineno} example/pm_dateTime/getMinute/main.F90
    !>  \compilef{getMinute}
    !>  \output{getMinute}
    !>  \include{lineno} example/pm_dateTime/getMinute/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getMinute}
    !>
    !>  \author
    !>  \AmirShahmoradi, January 30, 2021, 5:36 AM, Dallas, TX
    function getMinute() result(minute)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinute
#endif
        integer(IK) :: minute
        integer(IK) :: values(8)
        call date_and_time(values = values)
        CHECK_ASSERTION(__LINE__, all(values /= -huge(values)), SK_"@getMinute(): The processor does not provide date and time.") ! fpp
        minute = values(6)
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the current local hour of the current day of the Gregorian calendar, or the hour at the specified time zone `zone` in **minutes**.
    !>
    !>  \param[in]  zone    :   The output scalar of type `integer` of default kind \IK containing the time
    !>                          difference **in minutes** with respect to the Coordinated Universal Time (UTC).<br>
    !>
    !>  \return
    !>  `hour`              :   The output scalar of type `integer` of default kind \IK containing
    !>                          the current local hour of the current day of the Gregorian calendar.
    !>
    !>  \interface{getHour}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getHour
    !>      integer(IK) :: hour
    !>
    !>      hour = getHour() ! current hour at local time zone
    !>      hour = getHour(zone) ! current hour at the specified time zone in minutes
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The output value is set to `-huge(0_K)` if the processor does not provide date and time.<br>
    !>  \vericon
    !>
    !>  \see
    !>  [getMillisecond](@ref pm_dateTime::getMillisecond)<br>
    !>  [getSecond](@ref pm_dateTime::getSecond)<br>
    !>  [getMinute](@ref pm_dateTime::getMinute)<br>
    !>  [getHour12](@ref pm_dateTime::getHour12)<br>
    !>  [getHour](@ref pm_dateTime::getHour)<br>
    !>  [getZone](@ref pm_dateTime::getZone)<br>
    !>  [getZoneAbbr](@ref pm_dateTime::getZoneAbbr)<br>
    !>  [getDay](@ref pm_dateTime::getDay)<br>
    !>  [getMonth](@ref pm_dateTime::getMonth)<br>
    !>  [getYear](@ref pm_dateTime::getYear)<br>
    !>  [isLeapYear](@ref pm_dateTime::isLeapYear)<br>
    !>
    !>  \example{getHour}
    !>  \include{lineno} example/pm_dateTime/getHour/main.F90
    !>  \compilef{getHour}
    !>  \output{getHour}
    !>  \include{lineno} example/pm_dateTime/getHour/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getHour}
    !>
    !>  \author
    !>  \AmirShahmoradi, January 30, 2021, 5:36 AM, Dallas, TX
    impure elemental function getHour(zone) result(hour)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHour
#endif
        use pm_kind, only: IKG => IK
        integer(IKG), intent(in), optional  :: zone
        integer(IKG)                        :: hour
        integer(IKG)                        :: values(8)
        if (present(zone)) then
            CHECK_ASSERTION(__LINE__, isValidZone(zone), SK_"@getHour(): The condition `isValidZone(zone)` must hold. zone = "//getStr(zone)) ! fpp
            values = getDateTimeNewZone(zone)
        else
            values = getDateTime()
        end if
        !check_assertion(__LINE__, all(values /= -huge(values)), SK_"@getHour(): The processor does not provide date and time.") ! fpp
        hour = values(5)
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the local time difference **in minutes** with respect to the Coordinated Universal Time (UTC).
    !>
    !>  \return
    !>  `zone`  :   The output scalar of type `integer` of default kind \IK containing the local
    !>              time difference **in minutes** with respect to the Coordinated Universal Time (UTC).
    !>
    !>  \interface{getZone}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getZone
    !>      integer(IK) :: zone
    !>
    !>      zone = getZone()
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The output value is set to `-huge(0_K)` if the processor does not provide date and time.<br>
    !>  \vericon
    !>
    !>  \see
    !>  [getMillisecond](@ref pm_dateTime::getMillisecond)<br>
    !>  [getSecond](@ref pm_dateTime::getSecond)<br>
    !>  [getMinute](@ref pm_dateTime::getMinute)<br>
    !>  [getHour12](@ref pm_dateTime::getHour12)<br>
    !>  [getHour](@ref pm_dateTime::getHour)<br>
    !>  [getZone](@ref pm_dateTime::getZone)<br>
    !>  [getZoneAbbr](@ref pm_dateTime::getZoneAbbr)<br>
    !>  [getDay](@ref pm_dateTime::getDay)<br>
    !>  [getMonth](@ref pm_dateTime::getMonth)<br>
    !>  [getYear](@ref pm_dateTime::getYear)<br>
    !>  [isLeapYear](@ref pm_dateTime::isLeapYear)<br>
    !>
    !>  \example{getZone}
    !>  \include{lineno} example/pm_dateTime/getZone/main.F90
    !>  \compilef{getZone}
    !>  \output{getZone}
    !>  \include{lineno} example/pm_dateTime/getZone/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getZone}
    !>
    !>  \author
    !>  \AmirShahmoradi, January 30, 2021, 5:36 AM, Dallas, TX
    function getZone() result(zone)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getZone
#endif
        integer(IK) :: zone
        integer(IK) :: values(8)
        call date_and_time(values = values)
        CHECK_ASSERTION(__LINE__, all(values /= -huge(values)), SK_"@getZone(): The processor does not provide date and time.") ! fpp
        zone = values(4)
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the current day of the Gregorian calendar.
    !>
    !>  \return
    !>  `day`   :   The output scalar of type `integer` of default kind \IK containing
    !>              the current day of the Gregorian calendar.
    !>
    !>  \interface{getDay}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getDay
    !>      integer(IK) :: day
    !>
    !>      day = getDay()
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The output value is set to `-huge(0_K)` if the processor does not provide date and time.<br>
    !>  \vericon
    !>
    !>  \see
    !>  [getMillisecond](@ref pm_dateTime::getMillisecond)<br>
    !>  [getSecond](@ref pm_dateTime::getSecond)<br>
    !>  [getMinute](@ref pm_dateTime::getMinute)<br>
    !>  [getHour12](@ref pm_dateTime::getHour12)<br>
    !>  [getHour](@ref pm_dateTime::getHour)<br>
    !>  [getZone](@ref pm_dateTime::getZone)<br>
    !>  [getZoneAbbr](@ref pm_dateTime::getZoneAbbr)<br>
    !>  [getDay](@ref pm_dateTime::getDay)<br>
    !>  [getMonth](@ref pm_dateTime::getMonth)<br>
    !>  [getYear](@ref pm_dateTime::getYear)<br>
    !>  [isLeapYear](@ref pm_dateTime::isLeapYear)<br>
    !>
    !>  \example{getDay}
    !>  \include{lineno} example/pm_dateTime/getDay/main.F90
    !>  \compilef{getDay}
    !>  \output{getDay}
    !>  \include{lineno} example/pm_dateTime/getDay/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getDay}
    !>
    !>  \author
    !>  \AmirShahmoradi, January 30, 2021, 5:36 AM, Dallas, TX
    function getDay() result(day)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDay
#endif
        integer(IK) :: day
        integer(IK) :: values(8)
        call date_and_time(values = values)
        CHECK_ASSERTION(__LINE__, all(values /= -huge(values)), SK_"@getDay(): The processor does not provide date and time.") ! fpp
        day = values(3)
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the current month of the Gregorian calendar.
    !>
    !>  \return
    !>  `month` :   The output scalar of type `integer` of default kind \IK containing
    !>              the current month of the Gregorian calendar.
    !>
    !>  \interface{getMonth}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getMonth
    !>      integer(IK) :: month
    !>
    !>      month = getMonth()
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The output value is set to `-huge(0_K)` if the processor does not provide date and time.<br>
    !>  \vericon
    !>
    !>  \see
    !>  [getMillisecond](@ref pm_dateTime::getMillisecond)<br>
    !>  [getSecond](@ref pm_dateTime::getSecond)<br>
    !>  [getMinute](@ref pm_dateTime::getMinute)<br>
    !>  [getHour12](@ref pm_dateTime::getHour12)<br>
    !>  [getHour](@ref pm_dateTime::getHour)<br>
    !>  [getZone](@ref pm_dateTime::getZone)<br>
    !>  [getZoneAbbr](@ref pm_dateTime::getZoneAbbr)<br>
    !>  [getDay](@ref pm_dateTime::getDay)<br>
    !>  [getMonth](@ref pm_dateTime::getMonth)<br>
    !>  [getYear](@ref pm_dateTime::getYear)<br>
    !>  [isLeapYear](@ref pm_dateTime::isLeapYear)<br>
    !>
    !>  \example{getMonth}
    !>  \include{lineno} example/pm_dateTime/getMonth/main.F90
    !>  \compilef{getMonth}
    !>  \output{getMonth}
    !>  \include{lineno} example/pm_dateTime/getMonth/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getMonth}
    !>
    !>  \author
    !>  \AmirShahmoradi, January 30, 2021, 5:36 AM, Dallas, TX
    function getMonth() result(month)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMonth
#endif
        integer(IK) :: month
        integer(IK) :: values(8)
        call date_and_time(values = values)
        CHECK_ASSERTION(__LINE__, all(values /= -huge(values)), SK_"@getMonth(): The processor does not provide date and time.") ! fpp
        month = values(2)
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the current year of the Gregorian calendar.
    !>
    !>  \return
    !>  `year`  :   The output scalar of type `integer` of default kind \IK containing
    !>              the current year of the Gregorian calendar.
    !>
    !>  \interface{getYear}
    !>  \code{.F90}
    !>
    !>      use pm_dateTime, only: getYear
    !>      integer(IK) :: year
    !>
    !>      year = getYear()
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The output value is set to `-huge(0_K)` if the processor does not provide date and time.<br>
    !>  \vericon
    !>
    !>  \see
    !>  [getMillisecond](@ref pm_dateTime::getMillisecond)<br>
    !>  [getSecond](@ref pm_dateTime::getSecond)<br>
    !>  [getMinute](@ref pm_dateTime::getMinute)<br>
    !>  [getHour12](@ref pm_dateTime::getHour12)<br>
    !>  [getHour](@ref pm_dateTime::getHour)<br>
    !>  [getZone](@ref pm_dateTime::getZone)<br>
    !>  [getZoneAbbr](@ref pm_dateTime::getZoneAbbr)<br>
    !>  [getDay](@ref pm_dateTime::getDay)<br>
    !>  [getMonth](@ref pm_dateTime::getMonth)<br>
    !>  [getYear](@ref pm_dateTime::getYear)<br>
    !>  [isLeapYear](@ref pm_dateTime::isLeapYear)<br>
    !>
    !>  \example{getYear}
    !>  \include{lineno} example/pm_dateTime/getYear/main.F90
    !>  \compilef{getYear}
    !>  \output{getYear}
    !>  \include{lineno} example/pm_dateTime/getYear/main.out.F90
    !>
    !>  \test
    !>  [test_pm_dateTime](@ref test_pm_dateTime)
    !>
    !>  \final{getYear}
    !>
    !>  \author
    !>  \AmirShahmoradi, January 30, 2021, 5:36 AM, Dallas, TX
    function getYear() result(year)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getYear
#endif
        integer(IK) :: year
        integer(IK) :: values(8)
        call date_and_time(values = values)
        CHECK_ASSERTION(__LINE__, all(values /= -huge(values)), SK_"@getMonth(): The processor does not provide date and time.") ! fpp
        year = values(1)
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! LCOV_EXCL_START
    !>  \cond excluded
    !>  \brief
    !>  \legacy
    !>  Return date and time in a nice format.
    !>
    !>  \return
    !>  A character vector containing date and time in a nice format.
    !>
    !>  \remark
    !>  This is an impure function due to its dependence on `date_and_time()` Fortran intrinsic function.
    function getNiceDateTime() result(niceDateTime)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNiceDateTime
#endif
        implicit none
        character(len=21)                        :: niceDateTime
        character(10)                            :: time
        character(8)                             :: date
        call date_and_time(date,time)
        niceDateTime = date(1:4)//'/'//date(5:6)//'/'//date(7:8)//' - '//time(1:2)//':'//time(3:4)//':'//time(5:6)
    end function getNiceDateTime
    !>  \endcond excluded
    ! LCOV_EXCL_STOP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef  CHECK_ASSERTION

end module pm_dateTime ! LCOV_EXCL_LINE