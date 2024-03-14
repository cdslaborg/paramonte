program example

    use pm_kind, only: LK
    use pm_kind, only: IK
    use pm_kind, only: SK ! all other string  kinds are also supported.
    use pm_io, only: display_type
    use pm_val2logical, only: getLogical

    implicit none

    integer(IK) :: iostat(3)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Convert string values to logical.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getLogical('F')")
    call disp%show( getLogical('F') )
    call disp%skip()

    call disp%skip()
    call disp%show("getLogical('f')")
    call disp%show( getLogical('f') )
    call disp%skip()

    call disp%skip()
    call disp%show("getLogical('T')")
    call disp%show( getLogical('T') )
    call disp%skip()

    call disp%skip()
    call disp%show("getLogical('t')")
    call disp%show( getLogical('t') )
    call disp%skip()

    call disp%skip()
    call disp%show("getLogical('.F.')")
    call disp%show( getLogical('.F.') )
    call disp%skip()

    call disp%skip()
    call disp%show("getLogical('.T.')")
    call disp%show( getLogical('.T.') )
    call disp%skip()

    call disp%skip()
    call disp%show("getLogical('.false.')")
    call disp%show( getLogical('.false.') )
    call disp%skip()

    call disp%skip()
    call disp%show("getLogical('.true.')")
    call disp%show( getLogical('.true.') )
    call disp%skip()

    call disp%skip()
    call disp%show("getLogical('false')")
    call disp%show( getLogical('false') )
    call disp%skip()

    call disp%skip()
    call disp%show("getLogical('true')")
    call disp%show( getLogical('true') )
    call disp%skip()

    call disp%skip()
    call disp%show("getLogical('true, false')")
    call disp%show( getLogical('true, false') )
    call disp%skip()

    call disp%skip()
    call disp%show("getLogical([character(5,SK) :: 'true', 'false'])")
    call disp%show( getLogical([character(5,SK) :: 'true', 'false']) )
    call disp%skip()

    call disp%skip()
    call disp%show("getLogical([character(5,SK) :: 't', 'F', 'TRUE'])")
    call disp%show( getLogical([character(5,SK) :: 't', 'F', 'TRUE']) )
    call disp%skip()

    call disp%skip()
    call disp%show("getLogical('true, paramonte')")
    call disp%show( getLogical('true, paramonte') )
    call disp%skip()

end program example