program example

    use pm_kind, only: SK
    use pm_kind, only: IK ! All processor integer kinds are supported.
    use pm_io, only: display_type
    use pm_mathNumSys, only: getDecimal, getNumeral

    implicit none

    integer(IK) :: i

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate the numeral equivalent of the input decimal in the specified `base`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getDecimal(SK_'FFFFFF', 16_IK)")
    call disp%show( getDecimal(SK_'FFFFFF', 16_IK) )
    call disp%show("getNumeral(getDecimal(SK_'FFFFFF', 16_IK), 16_IK)")
    call disp%show( getNumeral(getDecimal(SK_'FFFFFF', 16_IK), 16_IK) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("huge(1_IK)")
    call disp%show( huge(1_IK) )
    call disp%show("getDecimal(SK_'1010011010', 2_IK)")
    call disp%show( getDecimal(SK_'1010011010', 2_IK) )
    call disp%show("getNumeral(getDecimal(SK_'1010011010', 2_IK), 2_IK)")
    call disp%show( getNumeral(getDecimal(SK_'1010011010', 2_IK), 2_IK) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getDecimal(getNumeral(huge(1_IK), base = 36_IK), 36_IK)")
    call disp%show( getDecimal(getNumeral(huge(1_IK), base = 36_IK), 36_IK) )
    call disp%show("getNumeral(huge(1_IK), base = 36_IK)")
    call disp%show( getNumeral(huge(1_IK), base = 36_IK) , deliml = SK_"""" )
    call disp%skip()

end program example