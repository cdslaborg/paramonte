program example

    use pm_kind, only: IK
    use pm_kind, only: LK, SK ! all processor types and kinds are supported.
    use pm_io, only: display_type
    use pm_val2int, only: getInt

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Convert logical values to integer.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getInt(.true._LK)")
    call disp%show( getInt(.true._LK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getInt(.false._LK)")
    call disp%show( getInt(.false._LK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getInt([.true._LK, .false._LK, .true._LK])")
    call disp%show( getInt([.true._LK, .false._LK, .true._LK]) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Convert string values to integer.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getInt(SK_'1')")
    call disp%show( getInt(SK_'1') )
    call disp%skip()

    call disp%skip()
    call disp%show("getInt(SK_'123456')")
    call disp%show( getInt(SK_'123456') )
    call disp%skip()

    call disp%skip()
    call disp%show("getInt(SK_'1 -1')")
    call disp%show( getInt(SK_'1 -1') )
    call disp%skip()

    call disp%skip()
    call disp%show("getInt(SK_'1, -1')")
    call disp%show( getInt(SK_'1, -1') )
    call disp%skip()

    call disp%skip()
    call disp%show("getInt(SK_'3, paramonte, (1, -1)')")
    call disp%show( getInt(SK_'3, paramonte, (1, -1)') )
    call disp%skip()

    call disp%skip()
    call disp%show("getInt(SK_'-3, (1, -1), paramonte')")
    call disp%show( getInt(SK_'-3, (1, -1), paramonte') )
    call disp%skip()

    call disp%skip()
    call disp%show("getInt([character(10,SK) :: '1, -1', '1, paramonte', '100, -100'])")
    call disp%show( getInt([character(10,SK) :: '1, -1', '1, paramonte', '100, -100']) )
    call disp%skip()

    call disp%skip()
    call disp%show("getInt(SK_'1, -1, paramonte')")
    call disp%show( getInt(SK_'1, -1, paramonte') )
    call disp%skip()

    call disp%skip()
    call disp%show("getInt(SK_'1, paramonte')")
    call disp%show( getInt(SK_'1, paramonte') )
    call disp%skip()

    call disp%skip()
    call disp%show("getInt(SK_'1')")
    call disp%show( getInt(SK_'1') )
    call disp%skip()

    call disp%skip()
    call disp%show("getInt([character(10,SK) :: '1, -1', '10', '100, -100'])")
    call disp%show( getInt([character(10,SK) :: '1, -1', '10', '100, -100']) )
    call disp%skip()

end program example