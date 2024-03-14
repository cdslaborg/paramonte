program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysPath, only: getPathMatch

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

#if __GFORTRAN__
    ! \bug
    call disp%warn%show("GNU gfortran compiler as of version 13.1, has a bug with the use of `spread()` for containers below, leading to content corruption.", width = 72_IK)
#endif

    call disp%skip()
    call disp%show("spread(getPathMatch(), 2, 1)")
    call disp%show( spread(getPathMatch(), 2, 1) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("spread(getPathMatch(key = SK_''), 2, 1)")
    call disp%show( spread(getPathMatch(key = SK_''), 2, 1) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("spread(getPathMatch(key = SK_'intel:oneapi'), 2, 1)")
    call disp%show( spread(getPathMatch(key = SK_'intel:oneapi'), 2, 1) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("spread(getPathMatch(key = SK_'intel:oneapi', sep = SK_':'), 2, 1)")
    call disp%show( spread(getPathMatch(key = SK_'intel:oneapi', sep = SK_':'), 2, 1) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("spread(getPathMatch(key = SK_'intel:oneapi', inc = SK_'mpiexec', sep = SK_':'), 2, 1)")
    call disp%show( spread(getPathMatch(key = SK_'intel:oneapi', inc = SK_'mpiexec', sep = SK_':'), 2, 1) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("spread(getPathMatch(key = SK_'intel:oneapi', inc = SK_'vars.sh', sep = SK_':'), 2, 1)")
    call disp%show( spread(getPathMatch(key = SK_'intel:oneapi', inc = SK_'vars.sh', sep = SK_':'), 2, 1) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("spread(getPathMatch(key = SK_'open:mpi', sep = SK_':'), 2, 1)")
    call disp%show( spread(getPathMatch(key = SK_'open:mpi', sep = SK_':'), 2, 1) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("spread(getPathMatch(key = SK_'mpich', sep = SK_':'), 2, 1)")
    call disp%show( spread(getPathMatch(key = SK_'mpich', sep = SK_':'), 2, 1) , deliml = SK_"""" )
    call disp%skip()

end program example