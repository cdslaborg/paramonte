program example

    use iso_fortran_env, only: error_unit
    use pm_kind, only: IK, LK, SK, RKS, RKD, RKQ
    use pm_bench, only: benchMulti_type
    use pm_io, only: display_type
    use pm_bench, only: bench_type

    implicit none

    type(benchMulti_type)   :: benchMulti
    real(RKS)               :: unifrnd_RKS
    real(RKS)               :: unifsum_RKS = 0._RKS ! dummy calculation to prevent aggressive compiler optimizations.
    real(RKD)               :: unifrnd_RKD
    real(RKD)               :: unifsum_RKD = 0._RKD ! dummy calculation to prevent aggressive compiler optimizations.
    real(RKQ)               :: unifrnd_RKQ
    real(RKQ)               :: unifsum_RKQ = 0._RKQ ! dummy calculation to prevent aggressive compiler optimizations.

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Benchmark the runtime cost of single, double, and quad precision exponentiaion.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("benchMulti = benchMulti_type([bench_type(SK_'rng_single', procwrap_RKS, overhead_RKS), bench_type(SK_'rng_double', procwrap_RKD, overhead_RKD), bench_type(SK_'rng_quadro', procwrap_RKQ, overhead_RKQ)])")
                    benchMulti = benchMulti_type([bench_type(SK_'rng_single', procwrap_RKS, overhead_RKS), bench_type(SK_'rng_double', procwrap_RKD, overhead_RKD), bench_type(SK_'rng_quadro', procwrap_RKQ, overhead_RKQ)])
                    !benchMulti = benchMulti_type([bench_type(SK_'rng_single', procwrap_RKS), bench_type(SK_'rng_double', procwrap_RKD), bench_type(SK_'rng_quadro', procwrap_RKQ)])
    call disp%show("call benchMulti%showsum(unit = disp%unit, tabular = .true._LK) ! Display the results in a nice tabular format.")
                    call benchMulti%showsum(unit = disp%unit, tabular = .true._LK) ! Display the results in a nice tabular format.
    call disp%show("call benchMulti%showsum(unit = disp%unit, tabular = .false._LK) ! Display the results in the default CSV format.")
                    call benchMulti%showsum(unit = disp%unit, tabular = .false._LK) ! Display the results in the default CSV format.
    call disp%skip()

    ! Print the dummy results to present potential aggressive compiler optimizations.

    write(error_unit, *) unifsum_RKS + unifsum_RKD + unifsum_RKQ

contains

    subroutine overhead_RKS(); call random_number(unifrnd_RKS); unifsum_RKS = unifsum_RKS + unifrnd_RKS; end
    subroutine overhead_RKD(); call random_number(unifrnd_RKD); unifsum_RKD = unifsum_RKD + unifrnd_RKD; end
    subroutine overhead_RKQ(); call random_number(unifrnd_RKQ); unifsum_RKQ = unifsum_RKQ + unifrnd_RKQ; end
    subroutine procwrap_RKS(); call random_number(unifrnd_RKS); unifsum_RKS = unifsum_RKS + 2._RKS**unifrnd_RKS; end
    subroutine procwrap_RKD(); call random_number(unifrnd_RKD); unifsum_RKD = unifsum_RKD + 2._RKD**unifrnd_RKD; end
    subroutine procwrap_RKQ(); call random_number(unifrnd_RKQ); unifsum_RKQ = unifsum_RKQ + 2._RKQ**unifrnd_RKQ; end

end program example