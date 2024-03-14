program example

    use pm_kind, only: IK
    use pm_kind, only: LK, SK ! all processor types and kinds are supported.
    use pm_io, only: display_type
    use pm_val2real, only: getReal

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Convert logical values to real.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getReal(.true._LK)")
    call disp%show( getReal(.true._LK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getReal(.false._LK)")
    call disp%show( getReal(.false._LK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getReal([.true._LK, .false._LK, .true._LK])")
    call disp%show( getReal([.true._LK, .false._LK, .true._LK]) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Convert string values to real.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getReal(SK_'1.e0')")
    call disp%show( getReal(SK_'1.e0') )
    call disp%skip()

    call disp%skip()
    call disp%show("getReal(SK_'123456.d0')")
    call disp%show( getReal(SK_'123456.d0') )
    call disp%skip()

    call disp%skip()
    call disp%show("getReal(SK_'1.d0 -1.d0')")
    call disp%show( getReal(SK_'1.d0 -1.d0') )
    call disp%skip()

    call disp%skip()
    call disp%show("getReal(SK_'1.d0, -1.d0')")
    call disp%show( getReal(SK_'1.d0, -1.d0') )
    call disp%skip()

    call disp%skip()
    call disp%show("getReal(SK_'3., paramonte, (1.d0, -1.e0)')")
    call disp%show( getReal(SK_'3., paramonte, (1.d0, -1.e0)') )
    call disp%skip()

    call disp%skip()
    call disp%show("getReal(SK_'-3., (1., -1.), paramonte')")
    call disp%show( getReal(SK_'-3., (1., -1.), paramonte') )
    call disp%skip()

    call disp%skip()
    call disp%show("getReal([character(10,SK) :: '1, -1', '1, paramonte', '100, -100'])")
    call disp%show( getReal([character(10,SK) :: '1, -1', '1, paramonte', '100, -100']) )
    call disp%skip()

    call disp%skip()
    call disp%show("getReal(SK_'1., -1., paramonte')")
    call disp%show( getReal(SK_'1., -1., paramonte') )
    call disp%skip()

    call disp%skip()
    call disp%show("getReal(SK_'1., paramonte')")
    call disp%show( getReal(SK_'1., paramonte') )
    call disp%skip()

    call disp%skip()
    call disp%show("getReal(SK_'1.d0')")
    call disp%show( getReal(SK_'1.d0') )
    call disp%skip()

    call disp%skip()
    call disp%show("getReal([character(10,SK) :: '1, -1', '10', '100, -100'])")
    call disp%show( getReal([character(10,SK) :: '1, -1', '10', '100, -100']) )
    call disp%skip()

end program example