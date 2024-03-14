program example

    use pm_kind, only: SK ! all processor types and kinds are supported.
    use pm_io, only: display_type
    use pm_str, only: getCharSeq

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getCharSeq([character(1,SK) :: 'P','a','r','a','M','o','n','t','e'])")
    call disp%show( getCharSeq([character(1,SK) :: 'P','a','r','a','M','o','n','t','e']) , deliml = SK_"""" )
    call disp%skip()

end program example