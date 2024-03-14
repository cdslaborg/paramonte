program example

    use pm_kind, only: IK, LK, SK, RKD
    use pm_io, only: display_type
    use pm_timer, only: timerCPU_type
    use pm_timer, only: timerDAT_type
    use pm_timer, only: timerSYS_type
    use pm_bench, only: bench_type

    implicit none

    type(bench_type)    :: bench
    real(RKD)           :: unifrnd
    real(RKD)           :: unifsum = 0._RKD ! dummy calculation to prevent aggressive compiler optimizations.

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Benchmark the uniform random number generation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("subroutine wrapper(); call random_number(unifrnd); unifsum = unifsum + unifrnd; end")
    call disp%show("bench = bench_type(name = SK_'random_number', exec = wrapper)")
                    bench = bench_type(name = SK_'random_number', exec = wrapper)
    call disp%show("bench%timing = bench%getTiming() ! same as below")
                    bench%timing = bench%getTiming()
    call disp%show("bench%timing%mean")
    call disp%show( bench%timing%mean )
    call disp%show("bench%timing%std")
    call disp%show( bench%timing%std )
    call disp%show("call bench%setTiming(minsec = 0.1_RKD) ! same as above but with a non-default minimum overall repetitive timing for 0.1 seconds.")
                    call bench%setTiming(minsec = 0.1_RKD)
    call disp%show("bench%timing%mean")
    call disp%show( bench%timing%mean )
    call disp%show("bench%timing%std")
    call disp%show( bench%timing%std )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Benchmark the uniform random number generation while excluding the overhead of redundant operations.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("subroutine wrapper(); call random_number(unifrnd); unifsum = unifsum + unifrnd; end")
    call disp%show("subroutine overhead(); unifsum = unifsum + unifrnd; end")
    call disp%show("bench = bench_type(name = SK_'random_number', exec = wrapper, overhead = overhead)")
                    bench = bench_type(name = SK_'random_number', exec = wrapper, overhead = overhead)
    call disp%show("bench%timing = bench%getTiming() ! same as below")
                    bench%timing = bench%getTiming()
    call disp%show("bench%timing%mean")
    call disp%show( bench%timing%mean )
    call disp%show("bench%timing%std")
    call disp%show( bench%timing%std )
    call disp%show("bench%timer%resol")
    call disp%show( bench%timer%resol )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Benchmark the uniform random number generation while excluding the overhead of redundant operations using a non-default CPU timer.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("subroutine wrapper(); call random_number(unifrnd); unifsum = unifsum + unifrnd; end")
    call disp%show("subroutine overhead(); unifsum = unifsum + unifrnd; end")
    call disp%show("bench = bench_type(name = SK_'random_number', exec = wrapper, overhead = overhead, timer = timerCPU_type())")
                    bench = bench_type(name = SK_'random_number', exec = wrapper, overhead = overhead, timer = timerCPU_type())
    call disp%show("bench%timing = bench%getTiming() ! same as below")
                    bench%timing = bench%getTiming()
    call disp%show("bench%timing%mean")
    call disp%show( bench%timing%mean )
    call disp%show("bench%timing%std")
    call disp%show( bench%timing%std )
    call disp%show("bench%timer%resol")
    call disp%show( bench%timer%resol )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Benchmark the uniform random number generation while excluding the overhead of redundant operations using a non-default date_and_time() timer.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("subroutine wrapper(); call random_number(unifrnd); unifsum = unifsum + unifrnd; end")
    call disp%show("subroutine overhead(); unifsum = unifsum + unifrnd; end")
    call disp%show("bench = bench_type(name = SK_'random_number', exec = wrapper, overhead = overhead, timer = timerDAT_type())")
                    bench = bench_type(name = SK_'random_number', exec = wrapper, overhead = overhead, timer = timerDAT_type())
    call disp%show("bench%timing = bench%getTiming() ! same as below")
                    bench%timing = bench%getTiming()
    call disp%show("bench%timing%mean")
    call disp%show( bench%timing%mean )
    call disp%show("bench%timing%std")
    call disp%show( bench%timing%std )
    call disp%show("bench%timer%resol")
    call disp%show( bench%timer%resol )
    call disp%skip()

contains

    impure subroutine showTimerComponents()
        call disp%show("bench%name")
        call disp%show( bench%name , deliml = SK_"""" )
        call disp%show("bench%minsec")
        call disp%show( bench%minsec )
        call disp%show("bench%miniter")
        call disp%show( bench%miniter )
    end

    subroutine wrapper()
        call random_number(unifrnd)
        unifsum = unifsum + unifrnd
    end

    subroutine overhead()
        unifsum = unifsum + unifrnd
    end

end program example