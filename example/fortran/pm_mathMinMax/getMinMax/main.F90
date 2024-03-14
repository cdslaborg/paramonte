program example

    use pm_kind, only: SK, IK, LK, CK, RK ! all intrinsic types and kinds are supported.
    use pm_io, only: display_type
    use pm_mathMinMax, only: getMinMax

    implicit none

    type(display_type) :: disp

    character(:, SK), allocatable   :: a_SK, b_SK, Pair_SK(:)
    integer(IK)                     :: a_IK, b_IK, Pair_IK(2)
    logical(LK)                     :: a_LK, b_LK, Pair_LK(2)
    complex(CK)                     :: a_CK, b_CK, Pair_CK(2)
    real(RK)                        :: a_RK, b_RK, Pair_RK(2)

    a_SK = "Hell"           ; b_SK = "Heaven"
    a_IK = 10_IK            ; b_IK = 5_IK
    a_LK = .true._LK        ; b_LK = .false._LK
    a_CK = (+10._CK, -5._CK); b_CK = (+5._CK, +15._CK)
    a_RK = 10._RK           ; b_RK = 5._RK

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Get the min/max of string values")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Pair_SK = [character(max(len(a_SK),len(b_SK)), SK) :: a_SK, b_SK]")
                    Pair_SK = [character(max(len(a_SK),len(b_SK)), SK) :: a_SK, b_SK]
    call disp%show("Pair_SK")
    call disp%show( Pair_SK , deliml = SK_"""" )
    call disp%show("getMinMax(Pair_SK)")
    call disp%show( getMinMax(Pair_SK) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("Pair_SK = [character(max(len(a_SK),len(b_SK)), SK) :: a_SK, b_SK]")
                    Pair_SK = [character(max(len(a_SK),len(b_SK)), SK) :: a_SK, b_SK]
    call disp%show("Pair_SK")
    call disp%show( Pair_SK , deliml = SK_"""" )
    call disp%show("getMinMax(Pair_SK(1), Pair_SK(2))")
    call disp%show( getMinMax(Pair_SK(1), Pair_SK(2)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Get the min/max of integer values")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Pair_IK = [a_IK, b_IK]")
                    Pair_IK = [a_IK, b_IK]
    call disp%show("Pair_IK")
    call disp%show( Pair_IK )
    call disp%show("getMinMax(Pair_IK)")
    call disp%show( getMinMax(Pair_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("Pair_IK = [a_IK, b_IK]")
                    Pair_IK = [a_IK, b_IK]
    call disp%show("Pair_IK")
    call disp%show( Pair_IK )
    call disp%show("getMinMax(Pair_IK(1), Pair_IK(2))")
    call disp%show( getMinMax(Pair_IK(1), Pair_IK(2)) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Get the min/max of logical values")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Pair_LK = [a_LK, b_LK]")
                    Pair_LK = [a_LK, b_LK]
    call disp%show("Pair_LK")
    call disp%show( Pair_LK )
    call disp%show("getMinMax(Pair_LK)")
    call disp%show( getMinMax(Pair_LK) )
    call disp%skip()

    call disp%skip()
    call disp%show("Pair_LK = [a_LK, b_LK]")
                    Pair_LK = [a_LK, b_LK]
    call disp%show("Pair_LK")
    call disp%show( Pair_LK )
    call disp%show("getMinMax(Pair_LK(1), Pair_LK(2))")
    call disp%show( getMinMax(Pair_LK(1), Pair_LK(2)) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Get the min/max of complex values")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Pair_CK = [a_CK, b_CK]")
                    Pair_CK = [a_CK, b_CK]
    call disp%show("Pair_CK")
    call disp%show( Pair_CK )
    call disp%show("getMinMax(Pair_CK)")
    call disp%show( getMinMax(Pair_CK) )
    call disp%skip()

    call disp%skip()
    call disp%show("Pair_CK = [a_CK, b_CK]")
                    Pair_CK = [a_CK, b_CK]
    call disp%show("Pair_CK")
    call disp%show( Pair_CK )
    call disp%show("getMinMax(Pair_CK(1), Pair_CK(2))")
    call disp%show( getMinMax(Pair_CK(1), Pair_CK(2)) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Get the min/max of real values")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Pair_RK = [a_RK, b_RK]")
                    Pair_RK = [a_RK, b_RK]
    call disp%show("Pair_RK")
    call disp%show( Pair_RK )
    call disp%show("getMinMax(Pair_RK)")
    call disp%show( getMinMax(Pair_RK) )
    call disp%skip()

    call disp%skip()
    call disp%show("Pair_RK = [a_RK, b_RK]")
                    Pair_RK = [a_RK, b_RK]
    call disp%show("Pair_RK")
    call disp%show( Pair_RK )
    call disp%show("getMinMax(Pair_RK(1), Pair_RK(2))")
    call disp%show( getMinMax(Pair_RK(1), Pair_RK(2)) )
    call disp%skip()

end program example