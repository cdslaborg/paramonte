program example

    use pm_kind, only: SK ! all processor types and kinds are supported.
    use pm_io, only: display_type
    use pm_str, only: getCharVec

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getCharVec(SK_'ParaMonte')")
    call disp%show( getCharVec(SK_'ParaMonte') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getCharVec(SK_'ParaMonte Library')")
    call disp%show( getCharVec(SK_'ParaMonte Library') , deliml = SK_"""" )
    call disp%skip()

end program example