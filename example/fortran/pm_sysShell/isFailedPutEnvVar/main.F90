program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysShell, only: isFailedPutEnvVar
    use pm_sysShell, only: isFailedGetEnvVar

    implicit none

    logical(LK) :: failed
    character(:, SK), allocatable :: name, value

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("name = SK_''; value = SK_''")
                    name = SK_''; value = SK_''
    call disp%show("isFailedPutEnvVar(name, value)")
    call disp%show( isFailedPutEnvVar(name, value) )
    call disp%skip()

    call disp%skip()
    call disp%show("name = SK_'PM'; value = SK_''")
                    name = SK_'PM'; value = SK_''
    call disp%show("isFailedGetEnvVar(name, value)")
    call disp%show( isFailedGetEnvVar(name, value) )
    call disp%show("value ! `value` should be empty.")
    call disp%show( value , deliml = SK_"""" )
    call disp%show("isFailedPutEnvVar(name, SK_'ParaMonte')")
    call disp%show( isFailedPutEnvVar(name, SK_'ParaMonte') )
    call disp%show("isFailedGetEnvVar(name, value)")
    call disp%show( isFailedGetEnvVar(name, value) )
    call disp%show("value ! `value` should be now `ParaMonte`.")
    call disp%show( value , deliml = SK_"""" )
    call disp%skip()

end program example