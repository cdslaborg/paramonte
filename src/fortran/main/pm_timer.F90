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
!>  This module contains the timer procedures and derived types to facilitate timing applications at runtime.<br>
!>
!>  \details
!>  Currently, five types of timers are available in this module in addition to a generic timer class.<br>
!>  <ol>
!>      <li>    [timer_type](@ref pm_timer::timer_type)
!>              <ol>
!>                  <li>    This is a generic timer class which defaults to the most appropriate timer type based on the ParaMonte library build.<br>
!>                          <ol>
!>                              <li>    If the ParaMonte library is built for serial applications,
!>                                      then [timerSYS_type](@ref pm_timer::timerSYS_type) is used.<br>
!>                              <li>    If the ParaMonte library is built for MPI-parallel applications,
!>                                      then [timerMPI_type](@ref pm_timer::timerMPI_type) is used.<br>
!>                              <li>    If the ParaMonte library is built for OpenMP-parallel applications,
!>                                      then [timerOMP_type](@ref pm_timer::timerOMP_type) is used.<br>
!>                          </ol>
!>                  <li>    Use this timer class if you are not sure which timer is suitable for the task.<br>
!>              </ol>
!>      <li>    [timerCPU_type](@ref pm_timer::timerCPU_type)
!>              <ol>
!>                  <li>    This timer relies on the Fortran intrinsic `cpu_time()`.<br>
!>                  <li>    This timer is **not ideal** for timing parallel or multi-threaded applications
!>                          because `cpu_time()` returns the sum of the total time spent by all processes.<br>
!>                  <li>    The timing precision of the Fortran intrinsic `cpu_time()` is processor-dependent.<br>
!>                  <li>    There is currently no direct way of measuring the time resolution of `cpu_time()`.<br>
!>                  <li>    This module uses a simple cycling method to compute the resolution of the time returned by `cpu_time()`.<br>
!>                          It computes the resolution as the least noticeable non-zero time difference between any two time points returned by `cpu_time()`.<br>
!>                  <li>    On processors with no clock, the time returned by `cpu_time()` is negative to indicate the lack of a processor clock.<br>
!>                          The clock existence is verified only if the library is built with the preprocessor macro `CHECK_ENABLED=1`.<br>
!>              </ol>
!>      <li>    [timerDAT_type](@ref pm_timer::timerDAT_type)
!>              <ol>
!>                  <li>    This timer relies on the Fortran intrinsic `date_and_time()`.<br>
!>                  <li>    This timer is **not ideal** for high resolution time measurements
!>                          because the resolution is limited to milliseconds.<br>
!>              </ol>
!>      <li>    [timerMPI_type](@ref pm_timer::timerMPI_type)
!>              <ol>
!>                  <li>    This timer relies on the MPI library intrinsic `mpi_wtime()`.<br>
!>                  <li>    The clock time resolution is computed by calling `mpi_wtick()`.<br>
!>                  <li>    This timer is **ideal** for high resolution (e.g., nanseconds) time measurements in distributed parallel applications.<br>
!>                  <li>    This timer is available only when the ParaMonte library is built with MPI support enabled.<br>
!>                  <li>    This timer requires the MPI library to have been already initialized via `mpi_init()`
!>                          and not have been finalized via `mpi_finalize()`.<br>
!>                          \vericons
!>              </ol>
!>      <li>    [timerOMP_type](@ref pm_timer::timerOMP_type)
!>              <ol>
!>                  <li>    This timer relies on the OpenMP library intrinsic `omp_get_wtime()`.<br>
!>                  <li>    The clock time resolution is computed by calling `omp_get_wtick()`.<br>
!>                  <li>    This timer is **ideal** for high resolution time measurements in multithreaded parallel applications.<br>
!>                  <li>    This timer is available only when the ParaMonte library is built with OpenMP support enabled.<br>
!>                  <li>    This timer measures elapsed wall clock time in seconds.<br>
!>                  <li>    The time is measured per thread.<br>
!>                  <li>    No guarantee can be made that two distinct threads measure the same time.<br>
!>                  <li>    Time is measured from **some time in the past**, which is an arbitrary time
!>                          guaranteed not to change during the execution of the program.<br>
!>              </ol>
!>      <li>    [timerSYS_type](@ref pm_timer::timerSYS_type)
!>              <ol>
!>                  <li>    This timer relies on the Fortran intrinsic `system_clock()`.<br>
!>                  <li>    This timer can be used in most parallel or serial timing scenarios.<br>
!>                  <li>    The resolution of this timer is processor-dependent.<br>
!>                  <li>    Time is measured from **some time in the past**, which is an arbitrary
!>                          time guaranteed not to change during the execution of the program.<br>
!>                  <li>    The time is measured per processor.<br>
!>                  <li>    On processors with no clock, `sysctem_clock()` returns a negative clock `count` to indicate the lack of a processor clock.<br>
!>                          The clock existence is verified only if the library is built with the preprocessor macro `CHECK_ENABLED=1`.<br>
!>              </ol>
!>  </ol>
!>
!>  \test
!>  [test_pm_timer](@ref test_pm_timer)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, March 22, 2012, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_timer

#if     CHECK_ENABLED
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
block; \
use pm_err, only: getFine; \
use pm_val2str, only: getStr; \
use pm_err, only: setAsserted; \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG); \
end block;
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_kind, only: SK, LK, IKD, RKD

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_timer"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the `abstract interface` of the [setTime](@ref pm_timer::timer_type) static type-bound procedure
    !>  component of [timer_type](@ref pm_timer::timer_type) `abstract` type that performs the timing since a user-specified or a processor-dependent origin.
    !>
    !>  \final{getTime_proc}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 23, 2012, 2:21 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    abstract interface
    function getTime_proc(since) result(timeInSec)
        use pm_kind, only: RKD
        real(RKD), intent(in), optional :: since
        real(RKD) :: timeInSec
    end function
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the `abstract interface` of the [wait](@ref pm_timer::timer_type) static type-bound procedure component of
    !>  [timer_type](@ref pm_timer::timer_type) `abstract` type that keeps the system waiting for the specified amount of seconds.<br>
    !>
    !>  \details
    !>  This procedure is functionally equivalent to the popular non-standard `sleep()` procedure offered by some compiler vendors.<br>
    !>
    !>  \note
    !>  Note that even though the system sleeps during the idle time, the timer clock is still ticking.<br>
    !>
    !>  \final{setIdle_proc}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 23, 2012, 2:21 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    abstract interface
    subroutine setIdle_proc(seconds)
        use pm_kind, only: RKD
        real(RKD), intent(in) :: seconds
    end subroutine
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the `abstract` base derived type that serves as a simple
    !>  container template for other timer classes in the ParaMonte library.
    !>
    !>  \details
    !>  Timer classes derived from based on this `abstract` type minimally contain components to store,
    !>  <ol>
    !>      <li>    the **start** time of a timer in seconds (since the processor-dependent epoch).
    !>      <li>    the **clock** in seconds (optionally computed since the start of the timer).
    !>      <li>    the **delta** elapsed time (optionally computed since the last call to the timer in seconds).
    !>      <li>    the **resolution** of the timer in seconds.
    !>  </ol>
    !>  See the documentation of [timer_typer](@ref pm_timer::timer_typer)
    !>  for the non-default class constructor interface for this type.<br>
    !>  See also the documentation details of [pm_timer](@ref pm_timer).
    !>
    !>  \interface{timer_type}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: timer_type
    !>      class(timer_type), allocatable :: timer
    !>
    !>      timer = timer_type()
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [timer_type](@ref pm_timer::timer_type)<br>
    !>  [timerDAT_type](@ref pm_timer::timerDAT_type)<br>
    !>  [timerMPI_type](@ref pm_timer::timerMPI_type)<br>
    !>  [timerOMP_type](@ref pm_timer::timerOMP_type)<br>
    !>  [timerSYS_type](@ref pm_timer::timerSYS_type)<br>
    !>
    !>  \example{timer_type}
    !>  \include{lineno} example/pm_timer/timer_type/main.F90
    !>  \compilef{timer_type}
    !>  \output{timer_type}
    !>  \include{lineno} example/pm_timer/timer_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \todo
    !>  \pvhigh
    !>  The `start` component of this derived type must become `protected` as soon as Fortran 2023 is supported by the compilers.<br>
    !>
    !>  \final{timer_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>
    type, abstract  :: timer_type
        real(RKD)   :: start    !<  \public The `protected` scalar `real` of kind `double precision` \RKD provided as a convenience for the user to contain the start time of the timer since the processor epoch (a processor-dependent past time) in seconds.<br>
                                !!          \warning The object of class [timer_type](@ref pm_timer::timer_type) must be initialized with the constructor (e.g., `timer = timer_type()`) to appropriately set the object `start` component to the processor epoch.<br>
        real(RKD)   :: clock    !<  \public The scalar `real` of kind `double precision` \RKD provided as a convenience for the user to contain the total time in seconds elapsed since the start of the timer.<br>
                                !!          The clock time can be readily can be computed as `timer%%clock = timer%%time(since = timer%%start)`.<br>
        real(RKD)   :: delta    !<  \public The scalar `real` of kind `double precision` \RKD provided as a convenience for the user to contain the delta time in seconds since the last timing (last timer call).<br>
                                !!          The delta time can be readily can be computed as `timer%%delta = timer%%time(since = timer%%start) - timer%%clock` where `timer%%clock` is the last clock read as `timer%%clock = timer%%time(since = timer%%start)`.<br>
        real(RKD)   :: resol    !<  \public The `protected` scalar `real` of kind `double precision` \RKD provided as a convenience for the user to contain the time in seconds between the timer clock tics.<br>
                                !!          \warning The object of class [timer_type](@ref pm_timer::timer_type) must be initialized with the constructor (e.g., `timer = timer_type()`) to appropriately set the object `resol` component to the clock time resolution.<br>
    contains
        procedure(getTime_proc), nopass, deferred   :: time !<  \public See [getTime_proc](@ref pm_timer::getTime_proc).
        procedure(setIdle_proc), nopass, deferred   :: wait !<  \public See [setIdle_proc](@ref pm_timer::setIdle_proc).
    end type

    interface timer_type
        module procedure :: timer_typer
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the `timerCPU_type` class, containing attributes and static methods for
    !>  setting up a timer based on the CPU-clock, using the Fortran intrinsic `cpu_time()`.
    !>
    !>  \details
    !>  See the documentation of [timerCPU_typer](@ref pm_timer::timerCPU_typer)
    !>  for the non-default constructor interface of this type.<br>
    !>  See also the documentation details of [pm_timer](@ref pm_timer).
    !>
    !>  \interface{timerCPU_type}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: timerCPU_type
    !>      type(timerCPU_type) :: timer
    !>
    !>      timer = timerCPU_type()
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [timer_type](@ref pm_timer::timer_type)<br>
    !>  [timerDAT_type](@ref pm_timer::timerDAT_type)<br>
    !>  [timerMPI_type](@ref pm_timer::timerMPI_type)<br>
    !>  [timerOMP_type](@ref pm_timer::timerOMP_type)<br>
    !>  [timerSYS_type](@ref pm_timer::timerSYS_type)<br>
    !>
    !>  \example{timerCPU_type}
    !>  \include{lineno} example/pm_timer/timerCPU_type/main.F90
    !>  \compilef{timerCPU_type}
    !>  \output{timerCPU_type}
    !>  \include{lineno} example/pm_timer/timerCPU_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \final{timerCPU_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    type, extends(timer_type)       :: timerCPU_type
    contains
        procedure, nopass           :: time => getTimeCPU   !<  \public See [getTime_proc](@ref pm_timer::getTime_proc).
        procedure, nopass           :: wait => setIdleCPU   !<  \public See [setIdle_proc](@ref pm_timer::setIdle_proc).
    end type

    interface timerCPU_type
        module procedure            :: timerCPU_typer
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the `timerDAT_type` class, containing attributes and static methods for
    !>  setting up a timer based on the system-clock, using the Fortran intrinsic `date_and_time()`.
    !>
    !>  \details
    !>  See the documentation of [timerDAT_typer](@ref pm_timer::timerDAT_typer)
    !>  for the non-default constructor interface of this type.<br>
    !>  See also the documentation details of [pm_timer](@ref pm_timer).
    !>
    !>  \interface{timerDAT_type}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: timerDAT_type
    !>      type(timerDAT_type) :: timer
    !>
    !>      timer = timerDAT_type()
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [timer_type](@ref pm_timer::timer_type)<br>
    !>  [timerDAT_type](@ref pm_timer::timerDAT_type)<br>
    !>  [timerMPI_type](@ref pm_timer::timerMPI_type)<br>
    !>  [timerOMP_type](@ref pm_timer::timerOMP_type)<br>
    !>  [timerSYS_type](@ref pm_timer::timerSYS_type)<br>
    !>
    !>  \example{timerDAT_type}
    !>  \include{lineno} example/pm_timer/timerDAT_type/main.F90
    !>  \compilef{timerDAT_type}
    !>  \output{timerDAT_type}
    !>  \include{lineno} example/pm_timer/timerDAT_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \final{timerDAT_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    type, extends(timer_type)       :: timerDAT_type
    contains
        procedure, nopass           :: time => getTimeDAT   !<  \public See [getTime_proc](@ref pm_timer::getTime_proc).
        procedure, nopass           :: wait => setIdleDAT   !<  \public See [setIdle_proc](@ref pm_timer::setIdle_proc).
    end type

    interface timerDAT_type
        module procedure            :: timerDAT_typer
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the `timerMPI_type` class, containing attributes and static
    !>  methods for setting up a timer based on the MPI intrinsic `MPI_Wtime()`.
    !>
    !>  \details
    !>  See the documentation of [timerMPI_typer](@ref pm_timer::timerMPI_typer)
    !>  for the non-default constructor interface of this type.<br>
    !>  See also the documentation details of [pm_timer](@ref pm_timer).
    !>
    !>  \interface{timerMPI_type}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: timerMPI_type
    !>      type(timerMPI_type) :: timer
    !>
    !>      timer = timerMPI_type()
    !>
    !>  \endcode
    !>
    !>  \note
    !>  This type is available only if the library is built with the preprocessor macro `MPI_ENABLED=1`.
    !>
    !>  \see
    !>  [timer_type](@ref pm_timer::timer_type)<br>
    !>  [timerDAT_type](@ref pm_timer::timerDAT_type)<br>
    !>  [timerMPI_type](@ref pm_timer::timerMPI_type)<br>
    !>  [timerOMP_type](@ref pm_timer::timerOMP_type)<br>
    !>  [timerSYS_type](@ref pm_timer::timerSYS_type)<br>
    !>
    !>  \example{timerMPI_type}
    !>  \include{lineno} example/pm_timer/timerMPI_type/main.F90
    !>  \compilef{timerMPI_type}
    !>  \output{timerMPI_type}
    !>  \include{lineno} example/pm_timer/timerMPI_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \final{timerMPI_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
#if MPI_ENABLED
    type, extends(timer_type)       :: timerMPI_type
    contains
        procedure, nopass           :: time => getTimeMPI   !<  \public See [getTime_proc](@ref pm_timer::getTime_proc).
        procedure, nopass           :: wait => setIdleMPI   !<  \public See [setIdle_proc](@ref pm_timer::setIdle_proc).
    end type

    interface timerMPI_type
        module procedure            :: timerMPI_typer
    end interface
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the `timerMPI_type` class, containing attributes and static
    !>  methods for setting up a timer based on the OpenMP intrinsic `omp_get_wtime()`.
    !>
    !>  \details
    !>  See the documentation of [timerOMP_typer](@ref pm_timer::timerOMP_typer)
    !>  for the non-default constructor interface of this type.<br>
    !>  See also the documentation details of [pm_timer](@ref pm_timer).
    !>
    !>  \interface{timerOMP_type}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: timerOMP_type
    !>      type(timerOMP_type) :: timer
    !>
    !>      timer = timerOMP_type()
    !>
    !>  \endcode
    !>
    !>  \note
    !>  This type is available only if the library is built with the preprocessor macro `OMP_ENABLED=1`.
    !>
    !>  \see
    !>  [timer_type](@ref pm_timer::timer_type)<br>
    !>  [timerDAT_type](@ref pm_timer::timerDAT_type)<br>
    !>  [timerMPI_type](@ref pm_timer::timerMPI_type)<br>
    !>  [timerOMP_type](@ref pm_timer::timerOMP_type)<br>
    !>  [timerSYS_type](@ref pm_timer::timerSYS_type)<br>
    !>
    !>  \example{timerOMP_type}
    !>  \include{lineno} example/pm_timer/timerOMP_type/main.F90
    !>  \compilef{timerOMP_type}
    !>  \output{timerOMP_type}
    !>  \include{lineno} example/pm_timer/timerOMP_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \final{timerOMP_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
#if OMP_ENABLED
    type, extends(timer_type)       :: timerOMP_type
    contains
        procedure, nopass           :: time => getTimeOMP   !<  \public See [getTime_proc](@ref pm_timer::getTime_proc).
        procedure, nopass           :: wait => setIdleOMP   !<  \public See [setIdle_proc](@ref pm_timer::setIdle_proc).
    end type

    interface timerOMP_type
        module procedure            :: timerOMP_typer
    end interface
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !!>  \brief
    !!>  This s the derived type for generating objects containing information about the system clock count.
    !!>
    !!>  \final{Count_type}
    !!>
    !!>  \author
    !!>  Amir Shahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    !type                            :: SysClockCount_type
    !    integer(IKD)                :: start = 0_IKD    !< The scalar `integer` of kind double precision \IKD containing the processor-dependent starting clock counts of the system.
    !    integer(IKD)                :: clock = 0_IKD    !< The scalar `integer` of kind double precision \IKD containing the total processor clock counts since `start`.
    !    integer(IKD)                :: delta = 0_IKD    !< The scalar `integer` of kind double precision \IKD containing total processor clock count since the last clock count measurement.
    !    integer(IKD)                :: max   = 0_IKD    !< The scalar `integer` of kind double precision \IKD containing maximum value that the processor count may take, or `0` if there is no clock.
    !    integer(IKD)                :: rate  = 0_IKD    !< The scalar `integer` of kind double precision \IKD containing number of clock counts per second, or `0` if there is no clock.
    !end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the `timerSYS_type` class, containing attributes and static methods for
    !>  setting up a timer based on the system-clock, using the Fortran intrinsic `system_clock()`.
    !>
    !>  \details
    !>  See the documentation of [timerSYS_typer](@ref pm_timer::timerSYS_typer)
    !>  for the non-default constructor interface of this type.<br>
    !>  See also the documentation details of [pm_timer](@ref pm_timer).
    !>
    !>  \interface{timerSYS_type}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: timerSYS_type
    !>      type(timerSYS_type) :: timer
    !>
    !>      timer = timerSYS_type()
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [timer_type](@ref pm_timer::timer_type)<br>
    !>  [timerDAT_type](@ref pm_timer::timerDAT_type)<br>
    !>  [timerMPI_type](@ref pm_timer::timerMPI_type)<br>
    !>  [timerOMP_type](@ref pm_timer::timerOMP_type)<br>
    !>  [timerSYS_type](@ref pm_timer::timerSYS_type)<br>
    !>
    !>  \example{timerSYS_type}
    !>  \include{lineno} example/pm_timer/timerSYS_type/main.F90
    !>  \compilef{timerSYS_type}
    !>  \output{timerSYS_type}
    !>  \include{lineno} example/pm_timer/timerSYS_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \final{timerSYS_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    type, extends(timer_type)       :: timerSYS_type
       !type(SysClockCount_type)    :: Count    !< An object of type [SysClockCount_type](@ref pm_timer::SysClockCount_type) containing information about the system clock.
    contains
        procedure, nopass           :: time => getTimeSYS   !<  \public See [getTime_proc](@ref pm_timer::getTime_proc).
        procedure, nopass           :: wait => setIdleSYS   !<  \public See [setIdle_proc](@ref pm_timer::setIdle_proc).
    end type

    interface timerSYS_type
        module procedure            :: timerSYS_typer
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the constructor of the [timer_type](@ref pm_timer::timer_type) class.<br>
    !>
    !>  \details
    !>  Upon return, the constructor initializes all components of the timer object.<br>
    !>  See also the documentation details of [pm_timer](@ref pm_timer).<br>
    !>
    !>  \return
    !>  `timer` :   The output scalar object of class [timer_type](@ref pm_timer::timer_type).<br>
    !>              Specifically, the output of type,<br>
    !>              <ol>
    !>                  <li>    [timerMPI_type](@ref pm_timer::timerMPI_type) if the ParaMonte library is built with preprocessor macro `MPI_ENABLED` defined.
    !>                  <li>    [timerOMP_type](@ref pm_timer::timerOMP_type) if the ParaMonte library is built with preprocessor macro `OMP_ENABLED` defined.
    !>                  <li>    [timerSYS_type](@ref pm_timer::timerSYS_type) if the ParaMonte library is built is serial mode.
    !>              </ol>
    !>
    !>  \interface{timer_typer}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: timer_type
    !>      class(timer_type), allocatable :: timer
    !>
    !>      timer = timer_type()
    !>
    !>  \endcode
    !>
    !>  \remark
    !>  See the documentation of [timer_type](@ref pm_timer::timer_type) for example usage.
    !>
    !>  \final{timer_typer}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    function timer_typer() result(timer)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: timer_typer
#endif
#if     MPI_ENABLED
        type(timerMPI_type) :: timer
        timer = timerMPI_type()
#elif   OMP_ENABLED
        type(timerOMP_type) :: timer
        timer = timerOMP_type()
#else
        type(timerSYS_type) :: timer
        timer = timerSYS_type()
#endif
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the constructor of the [timerCPU_type](@ref pm_timer::timerCPU_type) class.<br>
    !>
    !>  \details
    !>  Upon return, the constructor initializes all components of the timer object.<br>
    !>  See also the documentation details of [pm_timer](@ref pm_timer).
    !>
    !>  \return
    !>  `timer` :   The output scalar object of class [timerCPU_type](@ref pm_timer::timerCPU_type).
    !>
    !>  \interface{timerCPU_typer}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: timerCPU_type
    !>      type(timerCPU_type) :: timer
    !>
    !>      timer = timerCPU_type()
    !>
    !>  \endcode
    !>
    !>  \remark
    !>  See the documentation of [timerCPU_type](@ref pm_timer::timerCPU_type) for example usage.
    !>
    !>  \final{timerCPU_typer}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    function timerCPU_typer() result(timer)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: timerCPU_typer
#endif
        type(timerCPU_type) :: timer
        timer%resol = getResTimerCPU()
        timer%start = timer%time()
        timer%delta = 0._RKD
        timer%clock = 0._RKD
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the constructor of the [timerDAT_type](@ref pm_timer::timerDAT_type) class.<br>
    !>
    !>  \details
    !>  Upon return, the constructor initializes all components of the timer object.<br>
    !>  See also the documentation details of [pm_timer](@ref pm_timer).
    !>
    !>  \interface{timerDAT_typer}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: timerDAT_type
    !>      type(timerDAT_type) :: timer
    !>
    !>      timer = timerDAT_type()
    !>
    !>  \endcode
    !>
    !>  \return
    !>  `timer` :   The output scalar object of class [timerDAT_type](@ref pm_timer::timerDAT_type).
    !>
    !>  \remark
    !>  See the documentation of [timerDAT_type](@ref pm_timer::timerDAT_type) for example usage.
    !>
    !>  \final{timerDAT_typer}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    function timerDAT_typer() result(timer)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: timerDAT_typer
#endif
        type(timerDAT_type) :: timer
        timer%resol = 0.001_RKD
        timer%start = timer%time()
        timer%delta = 0._RKD
        timer%clock = 0._RKD
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the constructor of the [timerMPI_type](@ref pm_timer::timerMPI_type) class.<br>
    !>
    !>  \details
    !>  Upon return, the constructor initializes all components of the timer object.<br>
    !>  See also the documentation details of [pm_timer](@ref pm_timer).
    !>
    !>  \interface{timerMPI_typer}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: timerMPI_type
    !>      type(timerMPI_type) :: timer
    !>
    !>      timer = timerMPI_type()
    !>
    !>  \endcode
    !>
    !>  \return
    !>  `timer` :   The output scalar object of class [timerMPI_type](@ref pm_timer::timerMPI_type).
    !>
    !>  \warning
    !>  This constructor requires the MPI library to have been initialized and not finalized.<br>
    !>  \vericon
    !>
    !>  \warning
    !>  This constructor is available only if the library is built with the preprocessor macro `MPI_ENABLED=1`.<br>
    !>
    !>  \remark
    !>  See the documentation of [timerMPI_type](@ref pm_timer::timerMPI_type) for example usage.
    !>
    !>  \final{timerMPI_typer}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 02:51 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
#if MPI_ENABLED
    function timerMPI_typer() result(timer)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: timerMPI_typer
#endif
        use mpi !mpi_f08, only: mpi_initialized, mpi_init, mpi_wtick
        type(timerMPI_type) :: timer
        logical :: initialized
        integer :: ierrMPI
!#if     CHECK_ENABLED
!        block
!            use pm_err, only: setAsserted
!            use pm_val2str, only: getStr
!            logical :: initialized, finalized
!            integer :: ierrMPI
!            call mpi_initialized(initialized, ierrMPI)
!            call setAsserted(ierrMPI /= 0 .and. initialized, MODULE_NAME//SK_"@timerMPI_typer(): The MPI library must be initialized before attempting to call the MPI timer.")
!            call mpi_finalized(finalized, ierrMPI)
!            call setAsserted(ierrMPI /= 0 .and. finalized, MODULE_NAME//SK_"@timerMPI_typer(): The MPI library must not be finalized prior to calling the MPI timer.")
!        end block
!#endif
        call mpi_initialized(initialized, ierrMPI)
        if (.not. initialized .and. ierrMPI == 0) call mpi_init(ierrMPI)
        if (ierrMPI /= 0) error stop MODULE_NAME//SK_"@timerMPI_typer(): Failed to initialize the MPI library."
        timer%resol = real(mpi_wtick(), RKD)
        timer%start = timer%time()
        timer%delta = 0._RKD
        timer%clock = 0._RKD
    end function
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the constructor of the [timerOMP_type](@ref pm_timer::timerOMP_type) class.<br>
    !>
    !>  \details
    !>  Upon return, the constructor initializes all components of the timer object.<br>
    !>  See also the documentation details of [pm_timer](@ref pm_timer).
    !>
    !>  \interface{timerOMP_typer}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: timerOMP_type
    !>      type(timerOMP_type) :: timer
    !>
    !>      timer = timerOMP_type()
    !>
    !>  \endcode
    !>
    !>  \return
    !>  `timer` :   The output scalar object of class [timerOMP_type](@ref pm_timer::timerOMP_type).
    !>
    !>  \warning
    !>  This constructor is available only if the library is built with the preprocessor macro `OMP_ENABLED=1`.<br>
    !>
    !>  \remark
    !>  See the documentation of [timerOMP_type](@ref pm_timer::timerOMP_type) for example usage.
    !>
    !>  \final{timerOMP_typer}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
#if OMP_ENABLED
    function timerOMP_typer() result(timer)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: timerOMP_typer
#endif
        use omp_lib
        type(timerOMP_type) :: timer
        timer%resol = real(omp_get_wtick(), RKD)
        timer%start = timer%time()
        timer%delta = 0._RKD
        timer%clock = 0._RKD
    end function
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the constructor of the [timerSYS_type](@ref pm_timer::timerSYS_type) class.<br>
    !>
    !>  \details
    !>  Upon return, the constructor initializes all components of the timer object.<br>
    !>  See also the documentation details of [pm_timer](@ref pm_timer).
    !>
    !>  \interface{timerSYS_typer}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: timerSYS_type
    !>      type(timerSYS_type) :: timer
    !>
    !>      timer = timerSYS_type()
    !>
    !>  \endcode
    !>
    !>  \return
    !>  `timer` :   The output scalar object of class [timerSYS_type](@ref pm_timer::timerSYS_type).
    !>
    !>  \remark
    !>  See the documentation of [timerSYS_type](@ref pm_timer::timerSYS_type) for example usage.
    !>
    !>  \final{timerSYS_typer}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 03:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    function timerSYS_typer() result(timer)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: timerSYS_typer
#endif
        type(timerSYS_type) :: timer
        timer%resol = getResTimerSYS()
        timer%start = timer%time()
        timer%delta = 0._RKD
        timer%clock = 0._RKD
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the time in units of seconds since the
    !>  specified time `since` or since an arbitrary processor-dependent time.
    !>
    !>  \details
    !>  This function relies on either of the following lower-level module functions depending on the ParaMonte library build:<br>
    !>  <ol>
    !>      <li>    [getTimeMPI](@ref pm_timer::getTimeMPI) for MPI-enabled parallel builds of the ParaMonte library.<br>
    !>      <li>    [getTimeOMP](@ref pm_timer::getTimeOMP) for OpenMP-enabled parallel builds of the ParaMonte library.<br>
    !>      <li>    [getTimeSYS](@ref pm_timer::getTimeSYS) for all other parallel or serial builds of the ParaMonte library.<br>
    !>  </ol>
    !>
    !>  \param[in]  since   :   The input scalar of type `real` of kind double precision \RKD representing the time origin (epoch) in seconds.<br>
    !>                          (**optional**. If not present, its value is processor-dependent but constant at runtime.)
    !>
    !>  \return
    !>  `timeInSec`         :   The output scalar of type `real` of kind double precision \RKD, containing the time in
    !>                          seconds since a specified time `since` or since an arbitrary processor-dependent time.
    !>
    !>  \interface{getTime}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: getTime
    !>      real(RKD) :: timeInSec
    !>      logical(LK) :: failed
    !>
    !>      timeInSec = getTime(since = since)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [timer_type](@ref pm_timer::timer_type)<br>
    !>
    !>  \example{getTime}
    !>  \include{lineno} example/pm_timer/getTime/main.F90
    !>  \compilef{getTime}
    !>  \output{getTime}
    !>  \include{lineno} example/pm_timer/getTime/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \final{getTime}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    function getTime(since) result(timeInSec)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getTime
#endif
        real(RKD), intent(in), optional :: since
        real(RKD) :: timeInSec
#if     MPI_ENABLED
        timeInSec = getTimeMPI(since)
#elif   OMP_ENABLED
        timeInSec = getTimeOMP(since)
#else
        timeInSec = getTimeSYS(since)
#endif
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the CPU time in units of seconds since the
    !>  specified time `since` or since an arbitrary processor-dependent time.
    !>
    !>  \details
    !>  See also the documentation details of [pm_timer](@ref pm_timer).
    !>
    !>  \param[in]  since   :   The input scalar of type `real` of kind double precision \RKD representing the time origin (epoch) in seconds.<br>
    !>                          (**optional**. If not present, its value is processor-dependent but constant at runtime.)
    !>
    !>  \return
    !>  `timeInSec`         :   The output scalar of type `real` of kind double precision \RKD, containing the time in
    !>                          seconds since a specified time `since` or since an arbitrary processor-dependent time.
    !>
    !>  \interface{getTimeCPU}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: getTimeCPU
    !>      real(RKD) :: timeInSec
    !>      logical(LK) :: failed
    !>
    !>      timeInSec = getTimeCPU(since = since)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [timerCPU_type](@ref pm_timer::timerCPU_type)<br>
    !>
    !>  \example{getTimeCPU}
    !>  \include{lineno} example/pm_timer/getTimeCPU/main.F90
    !>  \compilef{getTimeCPU}
    !>  \output{getTimeCPU}
    !>  \include{lineno} example/pm_timer/getTimeCPU/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \final{getTimeCPU}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    function getTimeCPU(since) result(timeInSec)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getTimeCPU
#endif
        real(RKD), intent(in), optional :: since
        real(RKD) :: timeInSec
        call cpu_time(timeInSec)
        CHECK_ASSERTION(__LINE__, timeInSec >= 0._RKD, SK_"@getTimeCPU(): The CPU does not have a clock.") ! fpp
        if (present(since)) timeInSec = timeInSec - since
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the calendrical clock time in units of seconds since the
    !>  specified time `since` or since the Gregorian calendrical time origin.
    !>
    !>  \details
    !>  See also the documentation details of [pm_timer](@ref pm_timer).
    !>
    !>  \param[in]  since   :   The input scalar of type `real` of kind double precision \RKD representing the time origin (epoch) in seconds.<br>
    !>                          (**optional**. If not present, its value is the Gregorian calendrical time origin.)
    !>
    !>  \return
    !>  `timeInSec`         :   The output scalar of type `real` of kind double precision \RKD,
    !>                          containing the time in seconds since the specified time origin.
    !>
    !>  \interface{getTimeDAT}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: getTimeDAT
    !>      real(RKD) :: timeInSec
    !>
    !>      timeInSec = getTimeDAT(since = since)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [timerDAT_type](@ref pm_timer::timerDAT_type)<br>
    !>
    !>  \example{getTimeDAT}
    !>  \include{lineno} example/pm_timer/getTimeDAT/main.F90
    !>  \compilef{getTimeDAT}
    !>  \output{getTimeDAT}
    !>  \include{lineno} example/pm_timer/getTimeDAT/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \final{getTimeDAT}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    function getTimeDAT(since) result(timeInSec)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getTimeDAT
#endif
        use pm_kind, only: IK, RKD
        use pm_dateTime, only: getJulianDay
        use pm_dateTime, only: SECONDS_PER_DAY
        real(RKD), intent(in), optional :: since
        real(RKD) :: timeInSec
        integer(IK) :: values(8)
        call date_and_time(values = values)
        timeInSec = getJulianDay(values) * SECONDS_PER_DAY
        if (present(since)) timeInSec = timeInSec - since
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the MPI clock time in units of seconds since the
    !>  specified time `since` or since an arbitrary time origin set by the MPI library.
    !>
    !>  \details
    !>  See also the documentation details of [pm_timer](@ref pm_timer).
    !>
    !>  \param[in]  since   :   The input scalar of type `real` of kind double precision \RKD representing the time origin (epoch) in seconds.<br>
    !>                          (**optional**. If not present, its value is an arbitrary time origin set by the MPI library.)
    !>
    !>  \return
    !>  `timeInSec`         :   The output scalar of type `real` of kind double precision \RKD,
    !>                          containing the time in seconds since the specified time origin.
    !>
    !>  \interface{getTimeMPI}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: getTimeMPI
    !>      real(RKD) :: timeInSec
    !>
    !>      timeInSec = getTimeMPI(since = since)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  This procedure exists only if the library is built with the preprocessor macro `MPI_ENABLED=1`.
    !>
    !>  \impure
    !>
    !>  \see
    !>  [timerMPI_type](@ref pm_timer::timerMPI_type)<br>
    !>
    !>  \example{getTimeMPI}
    !>  \include{lineno} example/pm_timer/getTimeMPI/main.F90
    !>  \compilef{getTimeMPI}
    !>  \output{getTimeMPI}
    !>  \include{lineno} example/pm_timer/getTimeMPI/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \final{getTimeMPI}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
#if MPI_ENABLED
    function getTimeMPI(since) result(timeInSec)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getTimeMPI
#endif
        use mpi !mpi_f08, only: mpi_wtime
        use pm_kind, only: IK, RKD
        real(RKD), intent(in), optional :: since
        real(RKD) :: timeInSec
        timeInSec = mpi_wtime()
        if (present(since)) timeInSec = timeInSec - since
    end function
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the OpenMP clock time in units of seconds since the
    !>  specified time `since` or since an arbitrary time origin set by the OpenMP library.
    !>
    !>  \details
    !>  See also the documentation details of [pm_timer](@ref pm_timer).
    !>
    !>  \param[in]  since   :   The input scalar of type `real` of kind double precision \RKD representing the time origin (epoch) in seconds.<br>
    !>                          (**optional**. If not present, its value is an arbitrary time origin set by the OpenMP library.)
    !>
    !>  \return
    !>  `timeInSec`         :   The output scalar of type `real` of kind double precision \RKD,
    !>                          containing the time in seconds since the specified time origin.
    !>
    !>  \interface{getTimeOMP}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: getTimeOMP
    !>      real(RKD) :: timeInSec
    !>
    !>      timeInSec = getTimeOMP(since = since)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  This procedure exists only if the library is built with the preprocessor macro `OMP_ENABLED=1`.
    !>
    !>  \impure
    !>
    !>  \see
    !>  [timerOMP_type](@ref pm_timer::timerOMP_type)<br>
    !>
    !>  \example{getTimeOMP}
    !>  \include{lineno} example/pm_timer/getTimeOMP/main.F90
    !>  \compilef{getTimeOMP}
    !>  \output{getTimeOMP}
    !>  \include{lineno} example/pm_timer/getTimeOMP/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \final{getTimeOMP}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
#if OMP_ENABLED
    function getTimeOMP(since) result(timeInSec)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getTimeOMP
#endif
        use pm_kind, only: IK, RKD
        use omp_lib, only: omp_get_wtime
        real(RKD), intent(in), optional :: since
        real(RKD) :: timeInSec
        timeInSec = omp_get_wtime()
        if (present(since)) timeInSec = timeInSec - since
    end function
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the system clock time in units of seconds since the
    !>  specified time `since` or since an arbitrary processor-dependent time origin.
    !>
    !>  \details
    !>  See also the documentation details of [pm_timer](@ref pm_timer).
    !>
    !>  \param[in]  since   :   The input scalar of type `real` of kind double precision \RKD representing the time origin (epoch) in seconds.<br>
    !>                          (**optional**. If not present, its value is processor-dependent but constant at runtime.)
    !>
    !>  \return
    !>  `timeInSec`         :   The output scalar of type `real` of kind double precision \RKD,
    !>                          containing the time in seconds since a processor-defined origin of time.
    !>
    !>  \interface{getTimeSYS}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: getTimeSYS
    !>      real(RKD) :: timeInSec
    !>
    !>      timeInSec = getTimeSYS(since = since)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [timerSYS_type](@ref pm_timer::timerSYS_type)<br>
    !>
    !>  \example{getTimeSYS}
    !>  \include{lineno} example/pm_timer/getTimeSYS/main.F90
    !>  \compilef{getTimeSYS}
    !>  \output{getTimeSYS}
    !>  \include{lineno} example/pm_timer/getTimeSYS/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \todo
    !>  \pmed
    !>  The performance of this procedure could be improved by avoiding the unnecessary retrieval of
    !>  the count rate and computing its inverse to get the period (which is fixed in the entire runtime).<br>
    !>  This requires redefining this procedure as a type-bound procedure of a parent timer type.<br>
    !>  However, preliminary benchmarks indicate that the division (instead of multiplication) slows down
    !>  the computation by at most 25 percent. On modern architecture, the extra cost is likely on the
    !>  of 25 CPU cycles, approximately 25 nanoseconds. This is likely negligible in most practical
    !>  scenarios. Here is a benchmark script for the cost of division vs. multiplication,
    !>  \code{.F90}
    !>
    !>      implicit none
    !>      double precision :: t0, t1, x(2), summ
    !>      integer :: v0(8), v1(8)
    !>      integer :: i
    !>
    !>      call cpu_time(t0)
    !>      !call date_and_time(values = v0)
    !>      do i = 1, 10**7
    !>          call random_number(x)
    !>          summ = summ + x(2) * x(1)
    !>      end do
    !>      call cpu_time(t1)
    !>      print *, summ, "prod", t1 - t0
    !>      !call date_and_time(values = v1)
    !>      !print *, summ, "prod", v1(7) - v0(7) + (v1(8) - v0(8)) / 1.d3
    !>
    !>      summ = 0.d0
    !>      call cpu_time(t0)
    !>      !call date_and_time(values = v0)
    !>      do i = 1, 10**7
    !>          call random_number(x)
    !>          summ = summ + x(2) / x(1)
    !>      end do
    !>      !call date_and_time(values = v1)
    !>      call cpu_time(t1)
    !>      print *, summ, t1 - t0
    !>      !print *, summ, "divi", v1(7) - v0(7) + (v1(8) - v0(8)) / 1.d3
    !>
    !>      end
    !>
    !>  \endcode
    !>  The above code yields highly similar time costs for multiplication vs. division operation in most benchmarks.
    !>
    !>  \final{getTimeSYS}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    function getTimeSYS(since) result(timeInSec)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getTimeSYS
#endif
        real(RKD), intent(in), optional :: since
        real(RKD) :: timeInSec
        real(RKD) :: count_rate
        integer(IKD) :: count
       !integer(IKD) :: cmax
        call system_clock(count, count_rate)!, count_max = cmax)
        CHECK_ASSERTION(__LINE__, count > 0_IKD .and. count_rate > 0._RKD, SK_"@getTimeSYS(): The system does not have a clock.") ! .and. cmax > 0_IKD
        timeInSec = count / count_rate
        if (present(since)) timeInSec = timeInSec - since
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the time resolution by the build-time default timer of the ParaMonte library
    !>  as used in [timer_type](@ref pm_timer::timer_type), in units of seconds.
    !>
    !>  \details
    !>  This function relies on either of the following lower-level module functions depending on the ParaMonte library build:<br>
    !>  <ol>
    !>      <li>    [getResTimerMPI](@ref pm_timer::getResTimerMPI) for MPI-enabled parallel builds of the ParaMonte library.<br>
    !>      <li>    [getResTimerOMP](@ref pm_timer::getResTimerOMP) for OpenMP-enabled parallel builds of the ParaMonte library.<br>
    !>      <li>    [getResTimerSYS](@ref pm_timer::getResTimerSYS) for all other parallel or serial builds of the ParaMonte library.<br>
    !>  </ol>
    !>
    !>  \return
    !>  `resolution`    :   The output scalar of type `real` of kind double precision \RKD,
    !>                      containing the time resolution of the clock in units of seconds.
    !>
    !>  \interface{getResTimer}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: getResTimer
    !>      real(RKD) :: resolution
    !>
    !>      resolution = getResTimer()
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [timer_type](@ref pm_timer::timer_type)<br>
    !>
    !>  \example{getResTimer}
    !>  \include{lineno} example/pm_timer/getResTimer/main.F90
    !>  \compilef{getResTimer}
    !>  \output{getResTimer}
    !>  \include{lineno} example/pm_timer/getResTimer/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \final{getResTimer}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    function getResTimer() result(resolution)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getResTimer
#endif
        real(RKD) :: resolution
#if     MPI_ENABLED
        resolution = getResTimerMPI()
#elif   OMP_ENABLED
        resolution = getResTimerOMP()
#else
        resolution = getResTimerSYS()
#endif
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the time resolution of the Fortran intrinsic `cpu_time()`,
    !>  used in [timerCPU_type](@ref pm_timer::timerCPU_type), in units of seconds.
    !>
    !>  \details
    !>  See also the documentation details of [pm_timer](@ref pm_timer).
    !>
    !>  \return
    !>  `resolution`    :   The output scalar of type `real` of kind double precision \RKD,
    !>                      containing the time resolution of the CPU clock in units of seconds.
    !>
    !>  \interface{getResTimerCPU}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: getResTimerCPU
    !>      real(RKD) :: resolution
    !>
    !>      resolution = getResTimerCPU()
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getResTimerSYS](@ref pm_timer::getResTimerSYS)<br>
    !>
    !>  \example{getResTimerCPU}
    !>  \include{lineno} example/pm_timer/getResTimerCPU/main.F90
    !>  \compilef{getResTimerCPU}
    !>  \output{getResTimerCPU}
    !>  \include{lineno} example/pm_timer/getResTimerCPU/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \final{getResTimerCPU}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    function getResTimerCPU() result(resolution)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getResTimerCPU
#endif
        real(RKD) :: resolution
        real(RKD) :: resolution_
        call cpu_time(time = resolution_)
        CHECK_ASSERTION(__LINE__, resolution_ >= 0._RKD, SK_"@getResTimerCPU(): The CPU does not have a clock.") ! fpp
        do
            call cpu_time(time = resolution)
            if (resolution == resolution_) cycle
            resolution = resolution - resolution_
            exit
        end do
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the time resolution of the Fortran intrinsic `date_and_time()`,
    !>  used in [timerDAT_type](@ref pm_timer::timerDAT_type), in units of seconds.
    !>
    !>  \details
    !>  See also the documentation details of [pm_timer](@ref pm_timer).
    !>
    !>  \return
    !>  `resolution`    :   The output scalar of type `real` of kind double precision \RKD,
    !>                      containing the time resolution of `date_and_time()` in units of seconds.
    !>
    !>  \interface{getResTimerDAT}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: getResTimerDAT
    !>      real(RKD) :: resolution
    !>
    !>      resolution = getResTimerDAT()
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \see
    !>  [getResTimerSYS](@ref pm_timer::getResTimerSYS)<br>
    !>
    !>  \example{getResTimerDAT}
    !>  \include{lineno} example/pm_timer/getResTimerDAT/main.F90
    !>  \compilef{getResTimerDAT}
    !>  \output{getResTimerDAT}
    !>  \include{lineno} example/pm_timer/getResTimerDAT/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \final{getResTimerDAT}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    pure function getResTimerDAT() result(resolution)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getResTimerDAT
#endif
        real(RKD) :: resolution
        resolution = 0.001_RKD
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the time resolution of the MPI intrinsic timer,
    !>  used in [timerMPI_type](@ref pm_timer::timerMPI_type), in units of seconds.
    !>
    !>  \details
    !>  See also the documentation details of [pm_timer](@ref pm_timer).
    !>
    !>  \return
    !>  `resolution`    :   The output scalar of type `real` of kind double precision \RKD,
    !>                      containing the time resolution of `mpi_wtime()` in units of seconds.
    !>
    !>  \interface{getResTimerMPI}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: getResTimerMPI
    !>      real(RKD) :: resolution
    !>
    !>      resolution = getResTimerMPI()
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  This procedure exists only if the library is built with the preprocessor macro `MPI_ENABLED=1`.
    !>
    !>  \see
    !>  [getResTimerSYS](@ref pm_timer::getResTimerSYS)<br>
    !>
    !>  \example{getResTimerMPI}
    !>  \include{lineno} example/pm_timer/getResTimerMPI/main.F90
    !>  \compilef{getResTimerMPI}
    !>  \output{getResTimerMPI}
    !>  \include{lineno} example/pm_timer/getResTimerMPI/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \final{getResTimerMPI}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
#if MPI_ENABLED
    function getResTimerMPI() result(resolution)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getResTimerMPI
#endif
        use mpi !mpi_f08, only: mpi_wtick
        real(RKD) :: resolution
        resolution = real(mpi_wtick(), RKD)
    end function
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the time resolution of the OpenMP intrinsic timer,
    !>  used in [timerOMP_type](@ref pm_timer::timerOMP_type), in units of seconds.
    !>
    !>  \details
    !>  See also the documentation details of [pm_timer](@ref pm_timer).
    !>
    !>  \return
    !>  `resolution`    :   The output scalar of type `real` of kind double precision \RKD,
    !>                      containing the time resolution of `mpi_wtime()` in units of seconds.
    !>
    !>  \interface{getResTimerOMP}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: getResTimerOMP
    !>      real(RKD) :: resolution
    !>
    !>      resolution = getResTimerOMP()
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  This procedure exists only if the library is built with the preprocessor macro `OMP_ENABLED=1`.
    !>
    !>  \see
    !>  [getResTimerSYS](@ref pm_timer::getResTimerSYS)<br>
    !>
    !>  \example{getResTimerOMP}
    !>  \include{lineno} example/pm_timer/getResTimerOMP/main.F90
    !>  \compilef{getResTimerOMP}
    !>  \output{getResTimerOMP}
    !>  \include{lineno} example/pm_timer/getResTimerOMP/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \final{getResTimerOMP}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
#if OMP_ENABLED
    function getResTimerOMP() result(resolution)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getResTimerOMP
#endif
        use omp_lib, only: omp_get_wtick
        real(RKD) :: resolution
        resolution = real(omp_get_wtick(), RKD)
    end function
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the time resolution of the Fortran intrinsic `system_clock()`,
    !>  used in [timerSYS_type](@ref pm_timer::timerSYS_type), in units of seconds.
    !>
    !>  \details
    !>  See also the documentation details of [pm_timer](@ref pm_timer).
    !>
    !>  \return
    !>  `resolution`        :   The output scalar of type `real` of kind double precision \RKD,
    !>                          containing the time resolution of the system clock in units of seconds.
    !>
    !>  \interface{getResTimerSYS}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: getResTimerSYS
    !>      real(RKD) :: resolution
    !>
    !>      resolution = getResTimerSYS()
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [timerSYS_type](@ref pm_timer::timerSYS_type)<br>
    !>
    !>  \example{getResTimerSYS}
    !>  \include{lineno} example/pm_timer/getResTimerSYS/main.F90
    !>  \compilef{getResTimerSYS}
    !>  \output{getResTimerSYS}
    !>  \include{lineno} example/pm_timer/getResTimerSYS/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \final{getResTimerSYS}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    function getResTimerSYS() result(resolution)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getResTimerSYS
#endif
        real(RKD) :: resolution
        call system_clock(count_rate = resolution)
        CHECK_ASSERTION(__LINE__, resolution > 0._RKD, SK_"@getResTimerSYS(): The system does not have a clock.") ! fpp
        resolution = 1._RKD / resolution
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Set the processor in sleep mode for the input requested number of seconds via the default timer of the ParaMonte library.
    !>
    !>  \details
    !>  This function relies on either of the following lower-level module functions depending on the ParaMonte library build:<br>
    !>  <ol>
    !>      <li>    [setIdleMPI](@ref pm_timer::setIdleMPI) for MPI-enabled parallel builds of the ParaMonte library.<br>
    !>      <li>    [setIdleOMP](@ref pm_timer::setIdleOMP) for OpenMP-enabled parallel builds of the ParaMonte library.<br>
    !>      <li>    [setIdleSYS](@ref pm_timer::setIdleSYS) for all other parallel or serial builds of the ParaMonte library.<br>
    !>  </ol>
    !>
    !>  \param[in]  seconds :   The input scalar of type `real` of kind double precision \RKD
    !>                          representing the number of seconds to put the system into sleep.
    !>
    !>  \interface{setIdle}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: setIdle
    !>      real(RKD) :: seconds
    !>
    !>      call setIdle(seconds)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [setIdleCPU](@ref pm_timer::setIdleCPU)<br>
    !>  [setIdleDAT](@ref pm_timer::setIdleDAT)<br>
    !>  [setIdleMPI](@ref pm_timer::setIdleMPI)<br>
    !>  [setIdleOMP](@ref pm_timer::setIdleOMP)<br>
    !>  [setIdleSYS](@ref pm_timer::setIdleSYS)<br>
    !>  [getTimeCPU](@ref pm_timer::getTimeCPU)<br>
    !>  [getTimeDAT](@ref pm_timer::getTimeDAT)<br>
    !>  [getTimeMPI](@ref pm_timer::getTimeMPI)<br>
    !>  [getTimeOMP](@ref pm_timer::getTimeOMP)<br>
    !>  [getTimeSYS](@ref pm_timer::getTimeSYS)<br>
    !>
    !>  \example{setIdle}
    !>  \include{lineno} example/pm_timer/setIdle/main.F90
    !>  \compilef{setIdle}
    !>  \output{setIdle}
    !>  \include{lineno} example/pm_timer/setIdle/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \final{setIdle}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    subroutine setIdle(seconds)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setIdle
#endif
        use pm_kind, only: RKD
        real(RKD)   , intent(in)    :: seconds
#if     MPI_ENABLED
        call setIdleMPI(seconds)
#elif   OMP_ENABLED
        call setIdleOMP(seconds)
#else
        call setIdleSYS(seconds)
#endif
    end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Set the processor in sleep mode for the input requested number of seconds via the Fortran intrinsic `cpu_time()`.
    !>
    !>  \param[in]  seconds :   The input scalar of type `real` of kind double precision \RKD
    !>                          representing the number of seconds to put the system into sleep.
    !>
    !>  \interface{setIdleCPU}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: setIdleCPU
    !>      real(RKD) :: seconds
    !>
    !>      call setIdleCPU(seconds)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  If the system does not possess a clock, the results are unpredictable and the procedure may enter an indefinite loop.<br>
    !>  However, this condition rarely occurs and can be readily check only once on each new platform before a production run.<br>
    !>  The lack of a system clock will lead to runtime error only if the library is built with the preprocessor macro `CHECK_ENABLED=1`.
    !>
    !>  \impure
    !>
    !>  \see
    !>  [setIdleCPU](@ref pm_timer::setIdleCPU)<br>
    !>  [setIdleDAT](@ref pm_timer::setIdleDAT)<br>
    !>  [setIdleMPI](@ref pm_timer::setIdleMPI)<br>
    !>  [setIdleOMP](@ref pm_timer::setIdleOMP)<br>
    !>  [setIdleSYS](@ref pm_timer::setIdleSYS)<br>
    !>  [getTimeCPU](@ref pm_timer::getTimeCPU)<br>
    !>  [getTimeDAT](@ref pm_timer::getTimeDAT)<br>
    !>  [getTimeMPI](@ref pm_timer::getTimeMPI)<br>
    !>  [getTimeOMP](@ref pm_timer::getTimeOMP)<br>
    !>  [getTimeSYS](@ref pm_timer::getTimeSYS)<br>
    !>
    !>  \example{setIdleCPU}
    !>  \include{lineno} example/pm_timer/setIdleCPU/main.F90
    !>  \compilef{setIdleCPU}
    !>  \output{setIdleCPU}
    !>  \include{lineno} example/pm_timer/setIdleCPU/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \final{setIdleCPU}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    subroutine setIdleCPU(seconds)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setIdleCPU
#endif
        use pm_kind, only: RKD
        real(RKD)   , intent(in)    :: seconds
        real(RKD)                   :: since
        since = getTimeCPU()
        do
            if (getTimeCPU() - since > seconds) exit
        end do
    end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Set the processor in sleep mode for the input requested number of seconds via the Fortran intrinsic `date_and_time()`.
    !>
    !>  \param[in]  seconds :   The input scalar of type `real` of kind double precision \RKD
    !>                          representing the number of seconds to put the system into sleep.
    !>
    !>  \interface{setIdleDAT}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: setIdleDAT
    !>      real(RKD) :: seconds
    !>
    !>      call setIdleDAT(seconds)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  If the system does not possess a clock, the results are unpredictable and the procedure may enter an indefinite loop.<br>
    !>  However, this condition rarely occurs and can be readily check only once on each new platform before a production run.<br>
    !>  The lack of a system clock will lead to runtime error only if the library is built with the preprocessor macro `CHECK_ENABLED=1`.
    !>
    !>  \impure
    !>
    !>  \see
    !>  [setIdleCPU](@ref pm_timer::setIdleCPU)<br>
    !>  [setIdleDAT](@ref pm_timer::setIdleDAT)<br>
    !>  [setIdleMPI](@ref pm_timer::setIdleMPI)<br>
    !>  [setIdleOMP](@ref pm_timer::setIdleOMP)<br>
    !>  [setIdleSYS](@ref pm_timer::setIdleSYS)<br>
    !>  [getTimeCPU](@ref pm_timer::getTimeCPU)<br>
    !>  [getTimeDAT](@ref pm_timer::getTimeDAT)<br>
    !>  [getTimeMPI](@ref pm_timer::getTimeMPI)<br>
    !>  [getTimeOMP](@ref pm_timer::getTimeOMP)<br>
    !>  [getTimeSYS](@ref pm_timer::getTimeSYS)<br>
    !>
    !>  \example{setIdleDAT}
    !>  \include{lineno} example/pm_timer/setIdleDAT/main.F90
    !>  \compilef{setIdleDAT}
    !>  \output{setIdleDAT}
    !>  \include{lineno} example/pm_timer/setIdleDAT/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \final{setIdleDAT}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    subroutine setIdleDAT(seconds)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setIdleDAT
#endif
        use pm_kind, only: RKD
        real(RKD)   , intent(in)    :: seconds
        real(RKD)                   :: since
        since = getTimeDAT()
        do
            if (getTimeDAT() - since > seconds) exit
        end do
    end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Set the processor in sleep mode for the input requested number of seconds via the MPI intrinsic timer `MPI_Wtime()`.
    !>
    !>  \param[in]  seconds :   The input scalar of type `real` of kind double precision \RKD
    !>                          representing the number of seconds to put the system into sleep.
    !>
    !>  \interface{setIdleMPI}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: setIdleMPI
    !>      real(RKD) :: seconds
    !>
    !>      call setIdleMPI(seconds)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  This procedure exists only if the library is built with the preprocessor macro `MPI_ENABLED=1`.<br>
    !>
    !>  \warning
    !>  Calling this procedure requires the MPI library to have been initialized and not finalized.<br>
    !>  \vericon
    !>
    !>  \impure
    !>
    !>  \see
    !>  [setIdleCPU](@ref pm_timer::setIdleCPU)<br>
    !>  [setIdleDAT](@ref pm_timer::setIdleDAT)<br>
    !>  [setIdleMPI](@ref pm_timer::setIdleMPI)<br>
    !>  [setIdleOMP](@ref pm_timer::setIdleOMP)<br>
    !>  [setIdleSYS](@ref pm_timer::setIdleSYS)<br>
    !>  [getTimeCPU](@ref pm_timer::getTimeCPU)<br>
    !>  [getTimeDAT](@ref pm_timer::getTimeDAT)<br>
    !>  [getTimeMPI](@ref pm_timer::getTimeMPI)<br>
    !>  [getTimeOMP](@ref pm_timer::getTimeOMP)<br>
    !>  [getTimeSYS](@ref pm_timer::getTimeSYS)<br>
    !>
    !>  \example{setIdleMPI}
    !>  \include{lineno} example/pm_timer/setIdleMPI/main.F90
    !>  \compilef{setIdleMPI}
    !>  \output{setIdleMPI}
    !>  \include{lineno} example/pm_timer/setIdleMPI/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \final{setIdleMPI}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
#if MPI_ENABLED
    subroutine setIdleMPI(seconds)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setIdleMPI
#endif
        use pm_kind, only: RKD
        real(RKD)   , intent(in)    :: seconds
        real(RKD)                   :: since
        since = getTimeMPI()
        do
            if (getTimeMPI() - since > seconds) exit
        end do
    end subroutine
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Set the processor in sleep mode for the input requested number of seconds via the OpenMP intrinsic timer `omp_get_wtime()`.
    !>
    !>  \param[in]  seconds :   The input scalar of type `real` of kind double precision \RKD
    !>                          representing the number of seconds to put the system into sleep.
    !>
    !>  \interface{setIdleOMP}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: setIdleOMP
    !>      real(RKD) :: seconds
    !>
    !>      call setIdleOMP(seconds)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  This procedure exists only if the library is built with the preprocessor macro `OMP_ENABLED=1`.<br>
    !>
    !>  \impure
    !>
    !>  \see
    !>  [setIdleCPU](@ref pm_timer::setIdleCPU)<br>
    !>  [setIdleDAT](@ref pm_timer::setIdleDAT)<br>
    !>  [setIdleMPI](@ref pm_timer::setIdleMPI)<br>
    !>  [setIdleOMP](@ref pm_timer::setIdleOMP)<br>
    !>  [setIdleSYS](@ref pm_timer::setIdleSYS)<br>
    !>  [getTimeCPU](@ref pm_timer::getTimeCPU)<br>
    !>  [getTimeDAT](@ref pm_timer::getTimeDAT)<br>
    !>  [getTimeMPI](@ref pm_timer::getTimeMPI)<br>
    !>  [getTimeOMP](@ref pm_timer::getTimeOMP)<br>
    !>  [getTimeSYS](@ref pm_timer::getTimeSYS)<br>
    !>
    !>  \example{setIdleOMP}
    !>  \include{lineno} example/pm_timer/setIdleOMP/main.F90
    !>  \compilef{setIdleOMP}
    !>  \output{setIdleOMP}
    !>  \include{lineno} example/pm_timer/setIdleOMP/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \final{setIdleOMP}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
#if OMP_ENABLED
    subroutine setIdleOMP(seconds)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setIdleOMP
#endif
        use pm_kind, only: RKD
        real(RKD)   , intent(in)    :: seconds
        real(RKD)                   :: since
        since = getTimeOMP()
        do
            if (getTimeOMP() - since > seconds) exit
        end do
    end subroutine
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Set the processor in sleep mode for the input requested number of seconds via the Fortran intrinsic `system_clock()`.
    !>
    !>  \param[in]  seconds :   The input scalar of type `real` of kind double precision \RKD
    !>                          representing the number of seconds to put the system into sleep.
    !>
    !>  \interface{setIdleSYS}
    !>  \code{.F90}
    !>
    !>      use pm_timer, only: setIdleSYS
    !>      real(RKD) :: seconds
    !>
    !>      call setIdleSYS(seconds)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  If the system does not possess a clock, the results are unpredictable and the procedure may enter an indefinite loop.<br>
    !>  However, this condition rarely occurs and can be readily check only once on each new platform before a production run.<br>
    !>  The lack of a system clock will lead to runtime error only if the library is built with the preprocessor macro `CHECK_ENABLED=1`.
    !>
    !>  \impure
    !>
    !>  \see
    !>  [setIdleCPU](@ref pm_timer::setIdleCPU)<br>
    !>  [setIdleDAT](@ref pm_timer::setIdleDAT)<br>
    !>  [setIdleMPI](@ref pm_timer::setIdleMPI)<br>
    !>  [setIdleOMP](@ref pm_timer::setIdleOMP)<br>
    !>  [setIdleSYS](@ref pm_timer::setIdleSYS)<br>
    !>  [getTimeCPU](@ref pm_timer::getTimeCPU)<br>
    !>  [getTimeDAT](@ref pm_timer::getTimeDAT)<br>
    !>  [getTimeMPI](@ref pm_timer::getTimeMPI)<br>
    !>  [getTimeOMP](@ref pm_timer::getTimeOMP)<br>
    !>  [getTimeSYS](@ref pm_timer::getTimeSYS)<br>
    !>
    !>  \example{setIdleSYS}
    !>  \include{lineno} example/pm_timer/setIdleSYS/main.F90
    !>  \compilef{setIdleSYS}
    !>  \output{setIdleSYS}
    !>  \include{lineno} example/pm_timer/setIdleSYS/main.out.F90
    !>
    !>  \test
    !>  [test_pm_timer](@ref test_pm_timer)
    !>
    !>  \final{setIdleSYS}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    subroutine setIdleSYS(seconds)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setIdleSYS
#endif
        use pm_kind, only: RKD
        real(RKD)   , intent(in)    :: seconds
        real(RKD)                   :: since
        since = getTimeSYS()
        do
            if (getTimeSYS() - since > seconds) exit
        end do
    end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_timer ! LCOV_EXCL_LINE
