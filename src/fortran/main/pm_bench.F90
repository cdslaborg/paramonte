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
!>  This module contains abstract interfaces and types that facilitate benchmarking of different procedures.
!>
!>  \details
!>  The primary useful benchmarking object from this module are<br>
!>  <ol>
!>      <li>    [bench_type](@ref pm_bench::bench_type) for benchmarking single or collections of procedures/tasks.<br>
!>      <li>    [benchMulti_type](@ref pm_bench::benchMulti_type) for facilitating benchmarks of collections of procedures/tasks.<br>
!>  </ol>
!>  See the documentations of the respective objects for example usage.<br>
!>
!>  \see
!>  [pm_timer](@ref pm_timer)<br>
!>
!>  \test
!>  [test_pm_bench](@ref test_pm_bench)
!>
!>  \final{pm_bench}
!>
!>  \author
!>  \AmirShahmoradi, Wednesday 4:13 AM, August 13, 2016, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_bench

    use pm_kind, only: SK, IK, LK, RKD
    use pm_timer, only: TimerCPU_type
    use pm_timer, only: TimerDAT_type
#if MPI_ENABLED
    use pm_timer, only: timerMPI_type
#endif
#if OMP_ENABLED
    use pm_timer, only: timerOMP_type
#endif
    use pm_timer, only: timerSYS_type
    use pm_timer, only: timer_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = SK_"@pm_bench"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the `abstract interface` of the [exec](@ref pm_bench::bench_type) static type-bound procedure pointer
    !>  component of [bench_type](@ref pm_bench::bench_type) derived type that points to the user-supplied procedure to be timed.
    !>
    !>  \final{exec_proc}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 23, 2012, 2:21 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    abstract interface
        subroutine exec_proc()
        end subroutine exec_proc
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !!>  \brief
    !!>  This is the base class for creating objects that hold the statistics of the timing information of a benchmark.
    !!>
    !!>  \remark
    !!>  This type should generally not be explicitly used outside the [pm_bench](@ref pm_bench).
    !!>  It merely serves to create the `stat` component of the [timing_type](@ref timing_type) class below.
    !!>
    !!>  \see
    !!>  [timing_type](@ref timing_type)<br>
    !!>
    !!>  \test
    !!>  [test_pm_bench](@ref test_pm_bench)
    !!>
    !!>  \final{Stat_type}
    !!>
    !!>  \author
    !!>  Amir Shahmoradi, Wednesday 4:13 AM, August 13, 2016, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    !type                :: Stat_type
    !    real(RKD)       :: min      = -huge(0._RKD) !<  @public The minimum of the timing vector of the benchmark (in seconds).
    !    real(RKD)       :: max      = -huge(0._RKD) !<  @public The maximum of the timing vector of the benchmark (in seconds).
    !    real(RKD)       :: std      = -huge(0._RKD) !<  @public The standard deviation of the timing vector of the benchmark (in seconds).
    !    real(RKD)       :: mean     = -huge(0._RKD) !<  @public The mean of the timing vector of the benchmark (in seconds).
    !   !real(RKD)       :: median   = -huge(0._RKD) !<  @public The median of the timing vector of the benchmark (in seconds).
    !   !real(RKD)       :: skewness = -huge(0._RKD) !<  @public The skewness of the timing vector of the benchmark.
    !   !real(RKD)       :: kurtosis = -huge(0._RKD) !<  @public The kurtosis of the timing vector of the benchmark.
    !end type Stat_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the base class for creating objects that hold the timing information of a benchmark and the timing statistics.
    !>
    !>  \remark
    !>  This type should generally not be explicitly used outside the [pm_bench](@ref pm_bench).<br>
    !>  It merely serves to create the `timing` component of the [benchBase_type](@ref benchBase_type) class below.<br>
    !>
    !>  \see
    !>  [bench_type](@ref pm_bench::bench_type)<br>
    !>  [timer_type](@ref pm_timer::timer_type)<br>
    !>  [timerDAT_type](@ref pm_timer::timerDAT_type)<br>
    !>  [timerMPI_type](@ref pm_timer::timerMPI_type)<br>
    !>  [timerOMP_type](@ref pm_timer::timerOMP_type)<br>
    !>  [timerSYS_type](@ref pm_timer::timerSYS_type)<br>
    !>
    !>  \test
    !>  [test_pm_bench](@ref test_pm_bench)
    !>
    !>  \final{timing_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Wednesday 4:13 AM, August 13, 2016, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type                                :: timing_type
       !type(Stat_type)                 :: stat                     !<  @public The object of type [Stat_type](@ref Stat_type) containing the statistics of the `values` component.
        real(RKD)                       :: overhead = 0._RKD        !<  @public The average timing overhead in units of seconds.
        real(RKD)                       :: min      = -huge(0._RKD) !<  @public The minimum of the timing vector of the benchmark (in seconds).
        real(RKD)                       :: max      = -huge(0._RKD) !<  @public The maximum of the timing vector of the benchmark (in seconds).
        real(RKD)                       :: std      = -huge(0._RKD) !<  @public The standard deviation of the timing vector of the benchmark (in seconds).
        real(RKD)                       :: mean     = -huge(0._RKD) !<  @public The mean of the timing vector of the benchmark (in seconds).
       !real(RKD)                       :: median   = -huge(0._RKD) !<  @public The median of the timing vector of the benchmark (in seconds).
       !real(RKD)                       :: skewness = -huge(0._RKD) !<  @public The skewness of the timing vector of the benchmark.
       !real(RKD)                       :: kurtosis = -huge(0._RKD) !<  @public The kurtosis of the timing vector of the benchmark.
        real(RKD)       , allocatable   :: values(:)                !<  @public The vector of timing results in units of seconds, corrected by the average overhead time.
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the base class for creating low-level benchmark objects.
    !>
    !>  \details
    !>  See also the [constructor](@ref pm_bench::bench_typer) of this type.<br>
    !>
    !>  \interface{benchBase_type}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: SK, IK, RKD
    !>      use pm_bench, only: benchBase_type
    !>      use pm_timer, only: timerCPU_type
    !>      use pm_timer, only: timerDAT_type
    !>      use pm_timer, only: timerMPI_type
    !>      use pm_timer, only: timerOMP_type
    !>      use pm_timer, only: timerSYS_type
    !>      type(benchBase_type) :: benchBase
    !>      character(:, SK) :: name
    !>      integer(IK) :: miniter
    !>      real(RKD) :: minsec
    !>
    !>      benchBase = benchBase_type(name, minsec = minsec, miniter = miniter, timer = timerCPU_type())
    !>      benchBase = benchBase_type(name, minsec = minsec, miniter = miniter, timer = timerDAT_type())
    !>      benchBase = benchBase_type(name, minsec = minsec, miniter = miniter, timer = timerMPI_type())
    !>      benchBase = benchBase_type(name, minsec = minsec, miniter = miniter, timer = timerOMP_type())
    !>      benchBase = benchBase_type(name, minsec = minsec, miniter = miniter, timer = timerSYS_type())
    !>
    !>  \endcode
    !>
    !>  \remark
    !>  See [benchBase_typer](@ref pm_bench::benchBase_typer) for the non-default constructor of this type.<br>
    !>
    !>  \note
    !>  Although it is possible, **this type is not meant to be directly used for benchmarking**.<br>
    !>  Instead, <b>use the child class [bench_type](@ref pm_bench::bench_type) which greatly simplifies benchmarking</b>.<br>
    !>
    !>  \see
    !>  [bench_type](@ref pm_bench::bench_type)<br>
    !>  [benchBase_typer](@ref pm_bench::benchBase_typer) (constructor of the type)<br>
    !>  [timer_type](@ref pm_timer::timer_type)<br>
    !>  [timerDAT_type](@ref pm_timer::timerDAT_type)<br>
    !>  [timerMPI_type](@ref pm_timer::timerMPI_type)<br>
    !>  [timerOMP_type](@ref pm_timer::timerOMP_type)<br>
    !>  [timerSYS_type](@ref pm_timer::timerSYS_type)<br>
    !>
    !>  \example{benchBase_type}
    !>  \include{lineno} example/pm_bench/benchBase_type/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_bench/benchBase_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_bench](@ref test_pm_bench)
    !>
    !>  \final{benchBase_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Wednesday 4:13 AM, August 13, 2016, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type                                        :: benchBase_type
        real(RKD)                               :: minsec = 0.05_RKD    !<  @public The minimum time in units of seconds that the benchmark should last. It could take longer, but not less than `minsec`.
        integer(IK)                             :: miniter = 1_IK       !<  @public The minimum number of timing of the user-specified wrapper procedure. It could run more, but not fewer than `miniter`.
        type(timing_type)                       :: timing               !<  @public The object of type [timing_type](@ref pm_bench::timing_type) containing the timing information and statistics of the benchmark.
        character(:, SK)    , allocatable       :: name                 !<  @public The name of the procedure to be timed.
        class(timer_type)   , allocatable       :: timer                !<  @public The `allocatable` component of `abstract` class [timer_type](@ref pm_timer::timer_type) used for internal timing.<br>
                                                                        !!          Public access to this component provided is provided solely for the convenience of the user when access to a timer is needed.<br>
                                                                        !!          Otherwise, it not meant to be directly accessed or manipulated.<br>
                                                                        !!          The concrete type of this class component is set by the user at runtime depending on their choice of timer.<br>
    end type

    !>  \cond excluded
    interface benchBase_type
        module procedure :: benchBase_typer
    end interface
    !>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Construct and return an object of type [benchBase_type](@ref pm_bench::benchBase_type).
    !>
    !>  \details
    !>  This is the constructor of the type [benchBase_type](@ref pm_bench::benchBase_type)
    !>  for creating objects that facilitate the benchmarking of arbitrary procedures.
    !>
    !>  \param[in]  name        :   The input scalar `character` of arbitrary length of default kind \SK,
    !>                              containing the benchmark name (typically the name of the procedure to be timed).
    !>  \param[in]  minsec      :   The input scalar of type `real` of kind double precision \RKD representing
    !>                              the minimum time in units of seconds that the overall benchmark should last.<br>
    !>                              (**optional**, default = `0.05` seconds)
    !>  \param[in]  miniter     :   The input scalar of type `integer` of default kind \IK representing
    !>                              the minimum number of timing the user-specified wrapper procedure repeatedly.<br>
    !>                              (**optional**, default = `1`)
    !>  \param[in]  timer       :   The input object of `abstract` class [timer_type](@ref pm_timer::timer_type) representing the benchmark timer.<br>
    !>                              This object can be one of the available timers in [pm_timer](@ref pm_timer):
    !>                              <ol>
    !>                                  <li>    [timerCPU_type](@ref pm_timer::timerCPU_type),
    !>                                  <li>    [timerDAT_type](@ref pm_timer::timerDAT_type),
    !>                                  <li>    [timerMPI_type](@ref pm_timer::timerMPI_type),
    !>                                  <li>    [timerOMP_type](@ref pm_timer::timerOMP_type),
    !>                                  <li>    [timerSYS_type](@ref pm_timer::timerSYS_type),
    !>                              </ol>
    !>                              or any other user-defined subclass of the `abstract` [timer_type](@ref pm_timer::timer_type) class.<br>
    !>                              (**optional**, default = [timer_type](@ref pm_timer::timer_type))
    !>
    !>  \return
    !>  `benchBase`             :   The output scalar object of type [benchBase_type](@ref pm_bench::benchBase_type).
    !>
    !>  \interface{benchBase_typer}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: SK, IK, RKD
    !>      use pm_bench, only: benchBase_type
    !>      use pm_timer, only: timer_type
    !>      use pm_timer, only: timerCPU_type
    !>      use pm_timer, only: timerDAT_type
    !>      use pm_timer, only: timerMPI_type
    !>      use pm_timer, only: timerOMP_type
    !>      use pm_timer, only: timerSYS_type
    !>      type(benchBase_type) :: benchBase
    !>      character(:, SK) :: name
    !>      integer(IK) :: miniter
    !>      real(RKD) :: minsec
    !>
    !>      benchBase = benchBase_type(name, minsec = minsec, miniter = miniter, timer = timer_type())
    !>      benchBase = benchBase_type(name, minsec = minsec, miniter = miniter, timer = timerCPU_type())
    !>      benchBase = benchBase_type(name, minsec = minsec, miniter = miniter, timer = timerDAT_type())
    !>      benchBase = benchBase_type(name, minsec = minsec, miniter = miniter, timer = timerMPI_type())
    !>      benchBase = benchBase_type(name, minsec = minsec, miniter = miniter, timer = timerOMP_type())
    !>      benchBase = benchBase_type(name, minsec = minsec, miniter = miniter, timer = timerSYS_type())
    !>
    !>  \endcode
    !>
    !>  \remark
    !>  See also [benchBase_type](@ref pm_bench::benchBase_type) for possible calling interfaces.<br>
    !>
    !>  \remark
    !>  Note that the timing of the user-specified procedure wrapper is performed **until both conditions related to `minsec` and `miniter` are satisfied**.<br>
    !>  In other words, the timing is performed for **at least** `minsec` seconds **and at least** `miniter` number of iterations.<br>
    !>  The default values for these two variables are are set so that the benchmark typically ends within `50 ms`,
    !>  unless the runtime of the user-specified procedure is longer than `50 ms`.<br>
    !>
    !>  \note
    !>  If the specified benchmark routine is to be called certain number of times,
    !>  set `minsec = 0.` and `miniter` to desired number of times the routine must be timed.<br>
    !>
    !>  \see
    !>  [benchBase_type](@ref pm_bench::benchBase_type)<br>
    !>  [timerDAT_type](@ref pm_timer::timerDAT_type)<br>
    !>  [timerMPI_type](@ref pm_timer::timerMPI_type)<br>
    !>  [timerOMP_type](@ref pm_timer::timerOMP_type)<br>
    !>  [timerSYS_type](@ref pm_timer::timerSYS_type)<br>
    !>  [timer_type](@ref pm_timer::timer_type)<br>
    !>
    !>  \example{benchBase_typer}
    !>  \include{lineno} example/pm_bench/benchBase_type/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_bench/benchBase_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_bench](@ref test_pm_bench)
    !>
    !>  \final{benchBase_typer}
    !>
    !>  \author
    !>  \AmirShahmoradi, Wednesday 4:13 AM, August 13, 2016, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface benchBase_typer
    module function benchBase_typer(name, minsec, miniter, timer) result(benchBase)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: benchBase_typer
#endif
        use pm_kind, only: SK, IK, RKD
        use pm_timer, only: timer_type
        character(*, SK)    , intent(in)                :: name
        real(RKD)           , intent(in)    , optional  :: minsec
        integer(IK)         , intent(in)    , optional  :: miniter
        class(timer_type)   , intent(in)    , optional  :: timer
        type(benchBase_type)                            :: benchBase
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the class for creating benchmark and performance-profiling objects.
    !>
    !>  \details
    !>  This type has two `public` methods both of which accomplish the same task of timing the user-specified procedure wrapper.<br>
    !>  However, one ([getTiming](@ref pm_bench::getTiming)) has a `function` interface where the output of the method can be clearly specified,
    !>  while the other ([setTiming](@ref pm_bench::setTiming)) has a `subroutine` interface where the output of the timing is implicitly assigned
    !>  to the `timing` component of the parent object of type `bench_type`.<br>
    !>  The former has an explicit clear calling syntax.<br>
    !>  The latter has a more concise syntax with potentially faster runtime performance
    !>  (which is likely irrelevant in almost all practical scenarios).<br>
    !>
    !>  See also the [constructor](@ref pm_bench::bench_typer) of this type.<br>
    !>
    !>  \interface{bench_type}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: SK, IK, RKD
    !>      use pm_bench, only: benchBase_type
    !>      use pm_timer, only: timerCPU_type
    !>      use pm_timer, only: timerDAT_type
    !>      use pm_timer, only: timerMPI_type
    !>      use pm_timer, only: timerOMP_type
    !>      use pm_timer, only: timerSYS_type
    !>      type(bench_type) :: bench
    !>      character(:, SK) :: name
    !>      integer(IK) :: miniter
    !>      real(RKD) :: minsec
    !>
    !>      benchBase = bench_type(name, exec, overhead = overhead, minsec = minsec, miniter = miniter, timer = timerCPU_type())
    !>      benchBase = bench_type(name, exec, overhead = overhead, minsec = minsec, miniter = miniter, timer = timerDAT_type())
    !>      benchBase = bench_type(name, exec, overhead = overhead, minsec = minsec, miniter = miniter, timer = timerMPI_type())
    !>      benchBase = bench_type(name, exec, overhead = overhead, minsec = minsec, miniter = miniter, timer = timerOMP_type())
    !>      benchBase = bench_type(name, exec, overhead = overhead, minsec = minsec, miniter = miniter, timer = timerSYS_type())
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [benchMulti_type](@ref pm_bench::benchMulti_type)<br>
    !>  [benchBase_type](@ref pm_bench::benchBase_type)<br>
    !>  [timing_type](@ref pm_bench::timing_type)<br>
    !>  [bench_type](@ref pm_bench::bench_type)<br>
    !>  [timer_type](@ref pm_timer::timer_type)<br>
    !>  [timerDAT_type](@ref pm_timer::timerDAT_type)<br>
    !>  [timerMPI_type](@ref pm_timer::timerMPI_type)<br>
    !>  [timerOMP_type](@ref pm_timer::timerOMP_type)<br>
    !>  [timerSYS_type](@ref pm_timer::timerSYS_type)<br>
    !>
    !>  \example{bench_type}
    !>  \include{lineno} example/pm_bench/bench_type/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_bench/bench_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_bench](@ref test_pm_bench)
    !>
    !>  \final{bench_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Wednesday 4:13 AM, August 13, 2016, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(benchBase_type)                       :: bench_type
       !logical(LK)                          , private  :: isdone = .false._LK          !<  \private    The scalar `logical` of default kind \LK that is set to `.true.` once the benchmark is over (used internally to finalize the pointer components of the type).
        procedure(exec_proc), pointer, nopass, private  :: exec => null()               !<  \private    Procedure Pointer to the user-provided wrapper with explicit interface [exec_proc](@ref pm_bench::exec_proc).
                                                                                        !!              The wrapper must wrap the user-specified procedure that is to be timed.
        procedure(exec_proc), pointer, nopass, private  :: overhead => null()           !<  \public     Procedure Pointer to the (user-provided) overhead wrapper with explicit interface [exec_proc](@ref pm_bench::exec_proc).
                                                                                        !!              The overhead wrapper is used to measure the overhead of the `exec` wrapper component besides
                                                                                        !!              the cost of calling the wrapped procedure with the `exec` wrapper component.
    contains
        final                                           :: finalizeBench                !<  \private    The finalization method of the class.
        procedure, pass                                 :: getTiming => getTimingMethod !<  \public     The generic method name pointing to `function`   [getTimingMethod](@ref pm_bench::getTiming) be called by the user to initiate the timing of the wrapper for procedure of interest.
        procedure, pass                                 :: setTiming => setTimingMethod !<  \public     The generic method name pointing to `subroutine` [setTimingMethod](@ref pm_bench::setTiming) be called by the user to initiate the timing of the wrapper for procedure of interest.
    end type

    !>  \cond excluded
    interface bench_type
        module procedure :: bench_typer  !<  This is the [constructor](@ref pm_bench::bench_typer) of objects of type [bench_type](@ref pm_bench::bench_type).
    end interface
    !>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Construct and return an object of type [bench_type](@ref pm_bench::bench_type).
    !>
    !>  \details
    !>  This is the constructor of the type [bench_type](@ref pm_bench::bench_type) for creating objects to
    !>  perform benchmarking of the user-specified procedures within the user-provided wrapper function.
    !>
    !>  \param[in]  name        :   The input scalar `character` of arbitrary length of default kind \SK,
    !>                              containing the benchmark name (typically the name of the procedure to be timed).
    !>  \param[in]  exec        :   The input **procedure pointer** with the `abstract interface` [exec_proc](@ref exec_proc)
    !>                              pointing to the user-defined wrapper procedure that calls the arbitrary procedure to be timed.
    !>  \param[in]  overhead    :   The input **procedure pointer** with the `abstract interface` [exec_proc](@ref exec_proc)
    !>                              pointing to a user-defined wrapper procedure that calls everything executed within
    !>                              wrapper procedure `exec()`, except the call to the procedure that is being timed.<br>
    !>                              This procedure pointer will be used to measure the overhead due to calling
    !>                              or executing any non-relevant statements within the wrapper procedure `exec()`.<br>
    !>                              (**optional**. If missing, then a default empty wrapper procedure will be used to compute the overhead.<br>
    !>                              Note the default overhead is likely optimized away in production builds of the library,
    !>                              effectively yielding zero overhead.)
    !>  \param[in]  minsec      :   The input scalar of type `real` of kind double precision \RKD representing
    !>                              the minimum time in units of seconds that the overall benchmark should last.<br>
    !>                              (**optional**. The default is set by the constructor of the superclass [benchBase_type](@ref pm_bench::benchBase_type).)
    !>  \param[in]  miniter     :   The input scalar of type `integer` of default kind \IK representing
    !>                              the minimum number of timing the user-specified wrapper procedure repeatedly.<br>
    !>                              (**optional**. The default is set by the constructor of the superclass [benchBase_type](@ref pm_bench::benchBase_type).)
    !>  \param[in]  timer       :   The input object of `abstract` class [timer_type](@ref pm_timer::timer_type) representing the benchmark timer.<br>
    !>                              This object can be one of the available timers in [pm_timer](@ref pm_timer):
    !>                              <ol>
    !>                                  <li>    [timerCPU_type](@ref pm_timer::timerCPU_type),
    !>                                  <li>    [timerDAT_type](@ref pm_timer::timerDAT_type),
    !>                                  <li>    [timerMPI_type](@ref pm_timer::timerMPI_type),
    !>                                  <li>    [timerOMP_type](@ref pm_timer::timerOMP_type),
    !>                                  <li>    [timerSYS_type](@ref pm_timer::timerSYS_type),
    !>                              </ol>
    !>                              or any other user-defined subclass of the `abstract` [timer_type](@ref pm_timer::timer_type) class.<br>
    !>                              (**optional**, default = [timer_type](@ref pm_timer::timer_type))
    !>
    !>  \return
    !>  `bench`                 :   The output scalar object of type [bench_type](@ref pm_bench::bench_type).
    !>
    !>  \interface{bench_typer}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: SK, IK, RKD
    !>      use pm_bench, only: benchBase_type
    !>      use pm_timer, only: timerCPU_type
    !>      use pm_timer, only: timerDAT_type
    !>      use pm_timer, only: timerMPI_type
    !>      use pm_timer, only: timerOMP_type
    !>      use pm_timer, only: timerSYS_type
    !>      use pm_timer, only: timer_type
    !>      type(bench_type) :: bench
    !>      character(:, SK) :: name
    !>      integer(IK) :: miniter
    !>      real(RKD) :: minsec
    !>
    !>      benchBase = bench_type(name, exec, overhead = overhead, minsec = minsec, miniter = miniter, timer = timer_type())
    !>      benchBase = bench_type(name, exec, overhead = overhead, minsec = minsec, miniter = miniter, timer = timerCPU_type())
    !>      benchBase = bench_type(name, exec, overhead = overhead, minsec = minsec, miniter = miniter, timer = timerDAT_type())
    !>      benchBase = bench_type(name, exec, overhead = overhead, minsec = minsec, miniter = miniter, timer = timerMPI_type())
    !>      benchBase = bench_type(name, exec, overhead = overhead, minsec = minsec, miniter = miniter, timer = timerOMP_type())
    !>      benchBase = bench_type(name, exec, overhead = overhead, minsec = minsec, miniter = miniter, timer = timerSYS_type())
    !>
    !>  \endcode
    !>
    !>  \note
    !>  If the specified benchmark routine is to be called certain number of times,
    !>  set `minsec = 0.` and `miniter` to desired number of times the routine must be timed.<br>
    !>
    !>  \see
    !>  [benchBase_type](@ref benchBase_type)<br>
    !>  [benchMulti_type](@ref pm_bench::benchMulti_type)<br>
    !>  [timerDAT_type](@ref pm_timer::timerDAT_type)<br>
    !>  [timerMPI_type](@ref pm_timer::timerMPI_type)<br>
    !>  [timerOMP_type](@ref pm_timer::timerOMP_type)<br>
    !>  [timerSYS_type](@ref pm_timer::timerSYS_type)<br>
    !>  [timer_type](@ref pm_timer::timer_type)<br>
    !>
    !>  \example{bench_typer}
    !>  \include{lineno} example/pm_bench/bench_type/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_bench/bench_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_bench](@ref test_pm_bench)
    !>
    !>  \final{bench_typer}
    !>
    !>  \author
    !>  \AmirShahmoradi, Wednesday 4:13 AM, August 13, 2016, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface bench_typer
    module function bench_typer(name, exec, overhead, minsec, miniter, timer) result(bench)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: bench_typer
#endif
        use pm_timer, only: timer_type
        character(*, SK)    , intent(in)            :: name
        procedure(exec_proc)                        :: exec
        procedure(exec_proc)            , optional  :: overhead
        real(RKD)           , intent(in), optional  :: minsec
        integer(IK)         , intent(in), optional  :: miniter
        class(timer_type)   , intent(in), optional  :: timer
        type(bench_type)                            :: bench
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return an object of type [timing_type](@ref timing_type) containing the benchmark timing information and statistics.
    !>
    !>  \details
    !>  This procedure is a method of the class [bench_type](@ref pm_bench::bench_type).
    !>
    !>  \param[inout]   self    :   The parent object of class [bench_type](@ref pm_bench::bench_type) (passed implicitly to the method).<br>
    !>  \param[in]      minsec  :   The input scalar `real` of kind double precision \RKD
    !>                              representing the minimum time in seconds to spend on timing repeatedly.<br>
    !>                              (**optional**. The default is set by the constructor of the superclass [benchBase_type](@ref pm_bench::benchBase_type).)
    !>  \param[in]      miniter :   The input scalar `integer` of default kind \IK representing
    !>                              the minimum number of iterations (repetitions) of the timing to perform.<br>
    !>                              (**optional**. The default is set by the constructor of the superclass [benchBase_type](@ref pm_bench::benchBase_type).)
    !>
    !>  \return
    !>  `timing`                :   The output object of class [timing_type](@ref timing_type) containing the resulting timing vector and statistics.<br>
    !>                              For convenience, you can assign this output to the `timing` component of the same [bench_type](@ref pm_bench::bench_type)
    !>                              object to which this [getTiming](@ref pm_bench::getTiming) method belongs.<br>
    !>                              See [bench_type](@ref pm_bench::bench_type) for example usage.<br>
    !>
    !>  \interface{getTiming}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: RKD
    !>      type(bench_type) :: bench
    !>      real(RKD) :: minsec
    !>
    !>      bench = bench_type(name, exec, overhead = overhead, minsec = minsec, timer = timer)
    !>      bench%timing = bench%getTiming(minsec = minsec, miniter = miniter)
    !>
    !>  \endcode
    !>
    !>  \remark
    !>  This `function` method of type [bench_type](@ref pm_bench::bench_type) has the same functionality as the `subroutine`
    !>  interface [setTiming](@ref pm_bench::setTiming) with the only difference that the output is explicit.
    !>
    !>  \note
    !>  If the specified benchmark routine is to be called certain number of times,
    !>  set `minsec = 0.` and `miniter` to desired number of times the routine must be timed.<br>
    !>
    !>  \see
    !>  [setTiming](@ref pm_bench::setTiming)<br>
    !>  [bench_type](@ref pm_bench::bench_type)<br>
    !>
    !>  \remark
    !>  See [bench_type](@ref pm_bench::bench_type) for example usage.
    !>
    !>  \todo
    !>  \phigh
    !>  The computation of the median, skewness, and kurtosis of the timing vector in the `stat` component must be implemented.
    !>
    !>  \test
    !>  [test_pm_bench](@ref test_pm_bench)
    !>
    !>  \final{getTiming}
    !>
    !>  \author
    !>  \AmirShahmoradi, Wednesday 4:13 AM, August 13, 2016, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getTiming
    module function getTimingMethod(self, minsec, miniter) result(timing)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getTimingMethod
#endif
        use pm_kind, only: IK, RKD
        class(bench_type)   , intent(inout)             :: self
        integer(IK)         , intent(in)    , optional  :: miniter
        real(RKD)           , intent(in)    , optional  :: minsec
        type(timing_type)                               :: timing
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Time the user-specified procedure wrapper in the parent object of type [bench_type](@ref pm_bench::bench_type)
    !>  and store the output benchmark timing information and statistics implicitly in the `timing` component of
    !>  the input/output [bench_type](@ref pm_bench::bench_type) object.<br>
    !>
    !>  \details
    !>  This procedure is a method of the class [bench_type](@ref pm_bench::bench_type).<br>
    !>
    !>  \details
    !>
    !>  \param[inout]   self    :   The parent object of class [bench_type](@ref pm_bench::bench_type) (passed implicitly to the method).<br>
    !>  \param[in]      miniter :   The input scalar `integer` of default kind \IK representing the minimum number of iterations (repetitions) of the timing to perform.<br>
    !>                              (**optional**. The default is set by the constructor of the superclass [benchBase_type](@ref pm_bench::benchBase_type).)
    !>  \param[in]      minsec  :   The input scalar `real` of kind double precision \RKD representing the minimum time in seconds to spend on timing repeatedly.<br>
    !>                              (**optional**. The default is set by the constructor of the superclass [benchBase_type](@ref pm_bench::benchBase_type).)
    !>
    !>  \interface{setTiming}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: RKD
    !>      type(bench_type) :: bench
    !>      real(RKD) :: minsec
    !>
    !>      bench = bench_type(name, exec, overhead = overhead, minsec = minsec, timer = timer)
    !>      call bench%setTiming(minsec = minsec, miniter = miniter)
    !>
    !>  \endcode
    !>
    !>  \remark
    !>  This `subroutine` method of type [bench_type](@ref pm_bench::bench_type) has the same functionality as the `function` interface
    !>  [getTiming](@ref pm_bench::getTiming) with the only difference that the result output is implicit and potentially slightly faster.
    !>
    !>  \note
    !>  If the specified benchmark routine is to be called certain number of times,
    !>  set `minsec = 0.` and `miniter` to desired number of times the routine must be timed.<br>
    !>
    !>  \see
    !>  [getTiming](@ref pm_bench::getTiming)<br>
    !>  [bench_type](@ref pm_bench::bench_type)<br>
    !>
    !>  \remark
    !>  See [bench_type](@ref pm_bench::bench_type) for example usage.
    !>
    !>  \todo
    !>  \phigh
    !>  The computation of the median, skewness, and kurtosis of the timing vector in the `stat` component must be implemented.<br>
    !>
    !>  \test
    !>  [test_pm_bench](@ref test_pm_bench)
    !>
    !>  \final{setTiming}
    !>
    !>  \author
    !>  \AmirShahmoradi, Wednesday 4:13 AM, August 13, 2016, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setTiming
    module subroutine setTimingMethod(self, minsec, miniter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setTimingMethod
#endif
        use pm_kind, only: IK, RKD
        class(bench_type)   , intent(inout)             :: self
        integer(IK)         , intent(in)    , optional  :: miniter
        real(RKD)           , intent(in)    , optional  :: minsec
        !RUN_TIMING(self%timing)
    end subroutine
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !!>  \brief
    !!>  This is the base class for creating objects that hold the timing information of multiple benchmark objects and the timing statistics.
    !!>
    !!>  \remark
    !!>  This type should generally not be explicitly used outside the [pm_bench](@ref pm_bench).<br>
    !!>  It merely serves to create the `timing` component of the [benchMulti_type](@ref benchMulti_type) class below.<br>
    !!>
    !!>  \see
    !!>  [benchMulti_type](@ref pm_bench::benchMulti_type)<br>
    !!>  [timerDAT_type](@ref pm_timer::timerDAT_type)<br>
    !!>  [timerMPI_type](@ref pm_timer::timerMPI_type)<br>
    !!>  [timerOMP_type](@ref pm_timer::timerOMP_type)<br>
    !!>  [timerSYS_type](@ref pm_timer::timerSYS_type)<br>
    !!>  [timer_type](@ref pm_timer::timer_type)<br>
    !!>
    !!>  \test
    !!>  [test_pm_bench](@ref test_pm_bench)
    !!>
    !!>  \final{timing_type}
    !!>
    !!>  \author
    !!>  Amir Shahmoradi, Wednesday 4:13 AM, August 13, 2016, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    !type                                :: timing_type
    !    real(RKD)                       :: overhead = 0._RKD        !<  \public The average timing overhead in units of seconds.
    !    real(RKD)                       :: min      = -huge(0._RKD) !<  \public The minimum of the timing vector of the benchmark (in seconds).
    !    real(RKD)                       :: max      = -huge(0._RKD) !<  \public The maximum of the timing vector of the benchmark (in seconds).
    !    real(RKD)                       :: std      = -huge(0._RKD) !<  \public The standard deviation of the timing vector of the benchmark (in seconds).
    !    real(RKD)                       :: mean     = -huge(0._RKD) !<  \public The mean of the timing vector of the benchmark (in seconds).
    !   !real(RKD)                       :: median   = -huge(0._RKD) !<  \public The median of the timing vector of the benchmark (in seconds).
    !   !real(RKD)                       :: skewness = -huge(0._RKD) !<  \public The skewness of the timing vector of the benchmark.
    !   !real(RKD)                       :: kurtosis = -huge(0._RKD) !<  \public The kurtosis of the timing vector of the benchmark.
    !    type(container) , allocatable   :: case(:)                  !<  \public The matrix of timing results in units of seconds, corrected by the average overhead time.
    !end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the class for creating object to perform multiple benchmarks and performance-profiling.
    !>
    !>  \details
    !>  This type facilitates comparison of performances of **multiple** user-specified procedure wrappers
    !>  by reducing the amount of code to be written and automating and randomizing the timing schemes.<br>
    !>  See also [benchMulti_typer](@ref pm_bench::benchMulti_typer), the type constructor.<br>
    !>
    !>  \interface{benchMulti_type}
    !>  \code{.F90}
    !>
    !>      type(benchMulti_type) :: self(:)
    !>      type(bench_type), allocatable :: case(:)
    !>
    !>      self = benchMulti_type(case(:))
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [benchMulti_typer](@ref pm_bench::benchMulti_typer)<br>
    !>  [showsum](@ref pm_bench::showsum)<br>
    !>  [bench_type](@ref pm_bench::bench_type)<br>
    !>  [timer_type](@ref pm_timer::timer_type)<br>
    !>  [timerDAT_type](@ref pm_timer::timerDAT_type)<br>
    !>  [timerMPI_type](@ref pm_timer::timerMPI_type)<br>
    !>  [timerOMP_type](@ref pm_timer::timerOMP_type)<br>
    !>  [timerSYS_type](@ref pm_timer::timerSYS_type)<br>
    !>
    !>  \example{benchMulti_type}
    !>  \include{lineno} example/pm_bench/benchMulti_type/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_bench/benchMulti_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_bench](@ref test_pm_bench)
    !>
    !>  \final{benchMulti_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Wednesday 4:13 AM, August 13, 2016, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type                                    :: benchMulti_type
        integer(IK)                         :: ncase                        !<  The scalar `integer` of default kind \IK
                                                                            !!  containing the number of user-specified benchmark cases.
        character(:, SK)    , allocatable   :: name                         !<  The `allocatable` scalar `character` of default kind \SK
                                                                            !!  containing the names of benchmark cases separated with `" vs."`.
        type(bench_type)    , allocatable   :: case(:)                      !<  The vector component object of type [bench_type](@ref pm_bench::bench_type)
                                                                            !!  of size `ncase * repeat` containing the Benchmark cases.<br>
       !type(timing_type)                   :: timing                       !<  The object of type [timing_type](@ref pm_bench::timing_type) containing
       !                                                                    !!  the timing information and statistics of the benchmark.
    contains
        procedure, pass                     :: showsum => showsum_          !<  The generic object method name that points to [showsum](@ref pm_bench::showsum).
    end type

    interface benchMulti_type
        module procedure :: benchMulti_typer
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Construct, perform multiple benchmarking, and return the multiple benchmarking
    !>  results as an object of type [benchMulti_type](@ref pm_bench::benchMulti_type).
    !>
    !>  \details
    !>  This is the constructor of the type [benchMulti_type](@ref benchMulti_type) for creating objects that automate,
    !>  benchmark, and compare the performances of **multiple** user-specified procedure wrappers.<br>
    !>
    !>  \param[inout]   case    :   The input `contiguous` vector of arbitrary size of type [bench_type](@ref pm_bench::bench_type)
    !>                              containing the multiple benchmark object instances corresponding to each user-specified procedure wrapper.<br>
    !>  \param[in]      repeat  :   The input scalar `integer` of default kind \IK, representing the number of times each benchmark must be done to obtain reliable results.<br>
    !>                              Setting `repeat` to a large number can significantly increase the overall runtime of the benchmark list, since each one will be repeated `repeat` times.<br>
    !>                              However, it can also lead to potentially less-biased and most accurate benchmark results and comparison.<br>
    !>                              (**optional**, default = `2_IK`, i.e., each benchmark will be repeated twice and the results will be averaged for each benchmark instance.)
    !>  \param[in]      sorted  :   The input scalar `logical` of default kind \LK.<br>
    !>                              If `.true.`, the input benchmark instances `case` will be called and timed in the same order as given.<br>
    !>                              Otherwise, the order of the calls to the benchmarks instances are randomized.<br>
    !>                              A `.false.` input value aids the removal of potential biases due to the orderly calling of benchmark instances.<br>
    !>                              (**optional**, default = `.false._LK`)
    !>
    !>  \return
    !>  `self`                  :   The output scalar object of type [benchMulti_type](@ref pm_bench::benchMulti_type).
    !>
    !>  \interface{benchMulti_typer}
    !>  \code{.F90}
    !>
    !>      type(benchMulti_type) :: self(:)
    !>      type(bench_type), allocatable :: case(:)
    !>
    !>      self = benchMulti_type(case(:), sorted = sorted, repeat = repeat)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [showsum](@ref pm_bench::showsum)<br>
    !>  [benchMulti_type](@ref pm_bench::benchMulti_type)<br>
    !>  [bench_type](@ref pm_bench::bench_type)<br>
    !>  [timer_type](@ref pm_timer::timer_type)<br>
    !>  [timerDAT_type](@ref pm_timer::timerDAT_type)<br>
    !>  [timerMPI_type](@ref pm_timer::timerMPI_type)<br>
    !>  [timerOMP_type](@ref pm_timer::timerOMP_type)<br>
    !>  [timerSYS_type](@ref pm_timer::timerSYS_type)<br>
    !>
    !>  \example{benchMulti_typer}
    !>  \include{lineno} example/pm_bench/benchMulti_type/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_bench/benchMulti_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_bench](@ref test_pm_bench)
    !>
    !>  \todo
    !>  \pvlow
    !>  The current construction of the `name` component of the output object relies on repeated allocation of `name`.<br>
    !>  This can be improved by removing the redundant allocation in future, although any performance benefits are questionable.<br>
    !>
    !>  \final{benchMulti_typer}
    !>
    !>  \author
    !>  \AmirShahmoradi, Wednesday 4:13 AM, August 13, 2016, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface benchMulti_typer
    module function benchMulti_typer(case, sorted, repeat) result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: benchMulti_typer
#endif
        use pm_kind, only: IK, LK
        type(benchMulti_type)                           :: self
        type(bench_type)    , intent(in), contiguous    :: case(:)
        logical(LK)         , intent(in), optional      :: sorted
        integer(IK)         , intent(in), optional      :: repeat
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Time the user-specified procedure wrappers in the `case` vector
    !>  component of the parent object of type [benchMulti_type](@ref pm_bench::benchMulti_type)
    !>  and store the output benchmark timing information and statistics implicitly in the `timing` component of the object.
    !>
    !>  \details
    !>  This procedure is a method of the class [benchMulti_type](@ref pm_bench::benchMulti_type).
    !>
    !>  \param[inout]   self        :   The input/output object of class [benchMulti_type](@ref pm_bench::benchMulti_type) (passed implicitly to the method).
    !>  \param[out]     unit        :   The input scalar of type `integer` of default kind \IK containing the output unit
    !>                                  (e.g., output_unit, or an external file unit) where the results should be displayed.<br>
    !>                                  (**optional**, default = `output_unit` taken from `iso_fortran_env` Fortran intrinsic module.)
    !>  \param[in]      tabular     :   The input scalar `logical` of default kind \LK. If `.true.`, then the results
    !>                                  will be output to the supplied file `unit` in simple ASCII tabular format.<br>
    !>                                  The non-tabular format is specially desirable for scenarios where the output
    !>                                  should be postprocessed by other software or in other programming languages.<br>
    !>                                  (**optional**, default = `.true.` if the input `unit` corresponds to `output_unit`
    !>                                  from `iso_fortran_env` intrinsic Fortran module and `.false.` otherwise.)
    !>
    !>  \interface{showsum}
    !>  \code{.F90}
    !>
    !>      type(benchMulti_type) :: self
    !>
    !>      call self%showsum(unit = unit, tabular = tabular)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [benchMulti_type](@ref pm_bench::benchMulti_type)<br>
    !>
    !>  \remark
    !>  See [benchMulti_type](@ref pm_bench::benchMulti_type) for example usage.
    !>
    !>  \test
    !>  [test_pm_bench](@ref test_pm_bench)
    !>
    !>  \final{showsum}
    !>
    !>  \author
    !>  \AmirShahmoradi, Wednesday 4:13 AM, August 13, 2016, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface showsum
    module subroutine showsum_(self, unit, tabular)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: showsum_
#endif
        class(benchMulti_type)  , intent(in)            :: self
        integer(IK)             , intent(in), optional  :: unit
        logical(LK)             , intent(in), optional  :: tabular
    end subroutine
    end interface

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine finalizeBench(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: finalizeBench
#endif
        type(bench_type), intent(inout) :: self
        !if (self%isdone) then
        nullify(self%overhead)
        nullify(self%exec)
        !end if
    end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Take nothing, [do nothing](https://www.youtube.com/watch?v=w0jcXD1m0W4), and return nothing.<br>
    !>
    !>  \details
    !>  This generic interface is simply an empty wrapper to a serve as proxy for subroutine call overhead in benchmarks.<br>
    !>
    !>  \impure
    !>
    !>  \see
    !>  [benchMulti_type](@ref pm_bench::benchMulti_type)<br>
    !>
    !>  \remark
    !>  See [benchMulti_type](@ref pm_bench::benchMulti_type) for example usage.
    !>
    !>  \test
    !>  [test_pm_bench](@ref test_pm_bench)
    !>
    !>  \final{doNothing}
    !>
    !>  \author
    !>  \AmirShahmoradi, Wednesday 4:13 AM, August 13, 2016, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    impure subroutine doNothing()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: doNothing
#endif
    end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_bench ! LCOV_EXCL_LINE