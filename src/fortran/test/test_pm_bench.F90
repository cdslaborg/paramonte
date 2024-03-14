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
!>  This module contains tests of the module [pm_bench](@ref pm_bench).
!>
!>  \todo
!>  \phigh 
!>  The following tests should be extended to include allocatable arrays of arbitrary bounds.
!>
!>  \fintest
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_bench

    use pm_bench
    use pm_test, only: test_type, LK
    use pm_timer, only: timerCPU_type
    use pm_timer, only: timerDAT_type
    use pm_timer, only: timerSYS_type
    use pm_kind, only: IK, RK, RKS, RKD, RKQ

    implicit none

    private
    public :: setTest
    type(test_type) :: test

    real(RKS)       :: unifrnd_RKS
    real(RKS)       :: unifsum_RKS = 0._RKS ! dummy calculation to prevent aggressive compiler optimizations.
    real(RKD)       :: unifrnd_RKD
    real(RKD)       :: unifsum_RKD = 0._RKD ! dummy calculation to prevent aggressive compiler optimizations.
    real(RKQ)       :: unifrnd_RKQ
    real(RKQ)       :: unifsum_RKQ = 0._RKQ ! dummy calculation to prevent aggressive compiler optimizations.

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)
        call test%run(test_constructBench, SK_"test_constructBench")
        call test%run(test_constructBenchMulti, SK_"test_constructBenchMulti")
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructBench() result(assertion)
        use pm_timer, only: timerCPU_type
        use pm_timer, only: timerDAT_type
        use pm_timer, only: timerSYS_type
        logical(LK) :: assertion
        type(timerCPU_type) :: TimerCPU
        type(timerDAT_type) :: TimerDAT
        type(timerSYS_type) :: TimerSYS
        type(bench_type) :: bench
        integer(IK) :: miniter_def
        integer(IK) :: miniter
        real(RKD) :: minsec_def
        real(RKD) :: minsec
        
        miniter = 100_IK
        minsec = 0.02_RKD
        miniter_def = 1
        minsec_def = 0.05_RKD
        assertion = .true._LK

        TimerCPU = timerCPU_type()
        TimerDAT = timerDAT_type()
        TimerSYS = timerSYS_type()

        bench = bench_type("testing", exec)
        assertion = assertion .and. logical(bench%name == SK_"testing", LK)
        call test%assert(assertion, line = int(__LINE__, IK))

        TimerSYS%start = TimerSYS%time()
        call bench%setTiming()
        TimerSYS%delta = TimerSYS%time(since = TimerSYS%start)
        !assertion = assertion .and. logical(TimerSYS%delta >= minsec_def, LK)
        call test%assert(assertion, line = int(__LINE__, IK))
        assertion = assertion .and. logical(size(bench%timing%values, kind = IK) >= miniter_def, LK)
        call test%assert(assertion, line = int(__LINE__, IK))

        bench = bench_type("testing", exec, overhead)
        assertion = assertion .and. logical(bench%name == SK_"testing", LK)
        call test%assert(assertion, line = int(__LINE__, IK))

        call bench%setTiming()
        assertion = assertion .and. logical(size(bench%timing%values, kind = IK) >= miniter_def, LK)
        call test%assert(assertion, line = int(__LINE__, IK))

        bench = bench_type("testing", exec, overhead, minsec = minsec)
        assertion = assertion .and. logical(bench%name == SK_"testing", LK)
        call test%assert(assertion, line = int(__LINE__, IK))
        assertion = assertion .and. logical(bench%minsec == minsec, LK)
        call test%assert(assertion, line = int(__LINE__, IK))

        call bench%setTiming()
        assertion = assertion .and. logical(size(bench%timing%values, kind = IK) >= miniter_def, LK)
        call test%assert(assertion, line = int(__LINE__, IK))

        bench = bench_type("testing", exec, overhead, minsec = minsec, miniter = miniter)
        assertion = assertion .and. logical(bench%name == SK_"testing", LK)
        call test%assert(assertion, line = int(__LINE__, IK))
        assertion = assertion .and. logical(bench%minsec == minsec, LK)
        call test%assert(assertion, line = int(__LINE__, IK))
        assertion = assertion .and. logical(bench%miniter == miniter, LK)
        call test%assert(assertion, line = int(__LINE__, IK))

        bench = bench_type("testing", exec, overhead, minsec = minsec, miniter = miniter, timer = timerCPU_type())
        assertion = assertion .and. logical(bench%name == SK_"testing", LK)
        call test%assert(assertion, line = int(__LINE__, IK))
        assertion = assertion .and. logical(bench%minsec == minsec, LK)
        call test%assert(assertion, line = int(__LINE__, IK))
        assertion = assertion .and. logical(bench%miniter == miniter, LK)
        call test%assert(assertion, line = int(__LINE__, IK))

        TimerCPU%start = TimerCPU%time()
        call bench%setTiming()
        TimerCPU%delta = TimerCPU%time(since = TimerCPU%start)
        !assertion = assertion .and. logical(TimerCPU%delta >= minsec, LK)
        call test%assert(assertion, line = int(__LINE__, IK))

        bench = bench_type("testing", exec, overhead, minsec = minsec, miniter = miniter, timer = timerDAT_type())
        assertion = assertion .and. logical(bench%name == SK_"testing", LK)
        call test%assert(assertion, line = int(__LINE__, IK))
        assertion = assertion .and. logical(bench%minsec == minsec, LK)
        call test%assert(assertion, line = int(__LINE__, IK))
        assertion = assertion .and. logical(bench%miniter == miniter, LK)
        call test%assert(assertion, line = int(__LINE__, IK))

        TimerDAT%start = TimerDAT%time()
        call bench%setTiming(miniter = miniter)
        TimerDAT%delta = TimerDAT%time(since = TimerDAT%start)
        !assertion = assertion .and. logical(TimerDAT%delta >= minsec, LK)
        !print *, TimerDAT%delta, minsec
        !print *, TimerDAT%resol, bench%timer%resol
        call test%assert(assertion, line = int(__LINE__, IK))
        assertion = assertion .and. logical(size(bench%timing%values, kind = IK) >= miniter, LK)
        call test%assert(assertion, line = int(__LINE__, IK))

        bench = bench_type("testing", exec, overhead, minsec = minsec, miniter = miniter, timer = timerSYS_type())
        assertion = assertion .and. logical(bench%name == SK_"testing", LK)
        call test%assert(assertion, line = int(__LINE__, IK))
        assertion = assertion .and. logical(bench%minsec == minsec, LK)
        call test%assert(assertion, line = int(__LINE__, IK))
        assertion = assertion .and. logical(bench%miniter == miniter, LK)
        call test%assert(assertion, line = int(__LINE__, IK))

        TimerSYS%start = TimerSYS%time()
        call bench%setTiming(miniter = miniter)
        TimerSYS%delta = TimerSYS%time(since = TimerSYS%start)
        !assertion = assertion .and. logical(TimerSYS%delta >= minsec_def, LK)
        call test%assert(assertion, line = int(__LINE__, IK))
        assertion = assertion .and. logical(size(bench%timing%values, kind = IK) >= miniter, LK)
        call test%assert(assertion, line = int(__LINE__, IK))

        bench%timing = bench%getTiming(miniter = miniter)
        assertion = assertion .and. logical(size(bench%timing%values, kind = IK) >= miniter, LK)
        call test%assert(assertion, line = int(__LINE__, IK))

        bench%timing = bench%getTiming(miniter = miniter)
        assertion = assertion .and. logical(size(bench%timing%values, kind = IK) >= miniter, LK)
        call test%assert(assertion, line = int(__LINE__, IK))

        bench%timing = bench%getTiming(minsec = minsec, miniter = miniter)
        assertion = assertion .and. logical(size(bench%timing%values, kind = IK) >= miniter, LK)
        call test%assert(assertion, line = int(__LINE__, IK))

        bench%timing = bench%getTiming(miniter = miniter)
        assertion = assertion .and. logical(size(bench%timing%values, kind = IK) >= miniter, LK)
        call test%assert(assertion, line = int(__LINE__, IK))

        call bench%setTiming(miniter = miniter)
        assertion = assertion .and. logical(size(bench%timing%values, kind = IK) >= miniter, LK)
        call test%assert(assertion, line = int(__LINE__, IK))

        call bench%setTiming(minsec = minsec)
        call test%assert(assertion, line = int(__LINE__, IK))

        call bench%setTiming(minsec = minsec, miniter = miniter)
        assertion = assertion .and. logical(size(bench%timing%values, kind = IK) >= miniter, LK)
        call test%assert(assertion, line = int(__LINE__, IK))

        miniter = 10001_IK
        call bench%setTiming(miniter = miniter)
        assertion = assertion .and. logical(size(bench%timing%values, kind = IK) >= miniter, LK)
        call test%assert(assertion, line = int(__LINE__, IK))

    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructBenchMulti() result(assertion)
        use pm_io, only: display_type
        logical(LK)             :: assertion
        type(display_type)      :: disp
        type(benchMulti_type)   :: benchMulti
     
        assertion = .true._LK

        disp = display_type(file = "test_constructBenchMulti.tmp")

        benchMulti= benchMulti_type([ bench_type(SK_'rng_single', procwrap_RKS, overhead_RKS) & ! LCOV_EXCL_LINE
                                    , bench_type(SK_'rng_double', procwrap_RKD, overhead_RKD) & ! LCOV_EXCL_LINE
                                    , bench_type(SK_'rng_quadro', procwrap_RKQ, overhead_RKQ) & ! LCOV_EXCL_LINE
                                    ])
        call test%assert(assertion, line = int(__LINE__, IK))
        call benchMulti%showsum(unit = disp%unit) ! Display the results in a nice tabular format.
        call test%assert(assertion, line = int(__LINE__, IK))
        call benchMulti%showsum(unit = disp%unit, tabular = .true._LK) ! Display the results in a nice tabular format.
        call test%assert(assertion, line = int(__LINE__, IK))
        call benchMulti%showsum(unit = disp%unit, tabular = .false._LK) ! Display the results in a nice tabular format.
        call test%assert(assertion, line = int(__LINE__, IK))
        benchMulti= benchMulti_type([ bench_type(SK_'rng_single', procwrap_RKS, overhead_RKS) & ! LCOV_EXCL_LINE
                                    , bench_type(SK_'rng_double', procwrap_RKD, overhead_RKD) & ! LCOV_EXCL_LINE
                                    , bench_type(SK_'rng_quadro', procwrap_RKQ, overhead_RKQ) & ! LCOV_EXCL_LINE
                                    ] & ! LCOV_EXCL_LINE
                                    , sorted = .false. & ! LCOV_EXCL_LINE
                                    )
        benchMulti= benchMulti_type([ bench_type(SK_'rng_single', procwrap_RKS, overhead_RKS) & ! LCOV_EXCL_LINE
                                    , bench_type(SK_'rng_double', procwrap_RKD, overhead_RKD) & ! LCOV_EXCL_LINE
                                    , bench_type(SK_'rng_quadro', procwrap_RKQ, overhead_RKQ) & ! LCOV_EXCL_LINE
                                    ] & ! LCOV_EXCL_LINE
                                    , sorted = .true. & ! LCOV_EXCL_LINE
                                    )
        benchMulti= benchMulti_type([ bench_type(SK_'rng_single', procwrap_RKS, overhead_RKS) & ! LCOV_EXCL_LINE
                                    , bench_type(SK_'rng_double', procwrap_RKD, overhead_RKD) & ! LCOV_EXCL_LINE
                                    , bench_type(SK_'rng_quadro', procwrap_RKQ, overhead_RKQ) & ! LCOV_EXCL_LINE
                                    ] & ! LCOV_EXCL_LINE
                                    , repeat = 1_IK & ! LCOV_EXCL_LINE
                                    )

    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine overhead_RKS(); call random_number(unifrnd_RKS); unifsum_RKS = unifsum_RKS + unifrnd_RKS; end
    subroutine overhead_RKD(); call random_number(unifrnd_RKD); unifsum_RKD = unifsum_RKD + unifrnd_RKD; end
    subroutine overhead_RKQ(); call random_number(unifrnd_RKQ); unifsum_RKQ = unifsum_RKQ + unifrnd_RKQ; end
    subroutine procwrap_RKS(); call random_number(unifrnd_RKS); unifsum_RKS = unifsum_RKS + 2._RKS**unifrnd_RKS; end
    subroutine procwrap_RKD(); call random_number(unifrnd_RKD); unifsum_RKD = unifsum_RKD + 2._RKD**unifrnd_RKD; end
    subroutine procwrap_RKQ(); call random_number(unifrnd_RKQ); unifsum_RKQ = unifsum_RKQ + 2._RKQ**unifrnd_RKQ; end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine exec()
        real :: summ
        summ = wasteTime()
    end subroutine

    subroutine overhead()
    end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function wasteTime() result(summ)
        real :: summ, Matrix(100,100)
        call random_number(Matrix)
        summ = sum(Matrix)
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_bench