program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: isStrReal

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isStrReal(SK_'1')")
    call disp%show( isStrReal(SK_'1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrReal(SK_'1 ')")
    call disp%show( isStrReal(SK_'1 ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrReal(SK_' 1')")
    call disp%show( isStrReal(SK_' 1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrReal(SK_'+1')")
    call disp%show( isStrReal(SK_'+1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrReal(SK_'-1')")
    call disp%show( isStrReal(SK_'-1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrReal(SK_'.2')")
    call disp%show( isStrReal(SK_'.2') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrReal(SK_'-.0')")
    call disp%show( isStrReal(SK_'-.0') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrReal([character(3,SK) :: '1.0', '+1.', '-1.'])")
    call disp%show( isStrReal([character(3,SK) :: '1.0', '+1.', '-1.']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrReal([character(7,SK) :: '1.0e-01', '+1.E+10', '-1.12D0', '-1.1d-1'])")
    call disp%show( isStrReal([character(7,SK) :: '1.0e-01', '+1.E+10', '-1.12D0', '-1.1d-1']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrReal([character(10,SK) :: ' 1.e-01', '+1.E+10  ', '-1.12D0 ', '-1.1d-1'])")
    call disp%show( isStrReal([character(10,SK) :: ' 1.e-01', '+1.E+10  ', '-1.12D0 ', '-1.1d-1']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrReal(SK_' a1. ')")
    call disp%show( isStrReal(SK_' a1. ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrReal(SK_' 1.a ')")
    call disp%show( isStrReal(SK_' 1.a ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrReal(SK_' a ')")
    call disp%show( isStrReal(SK_' a ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrReal(SK_' a')")
    call disp%show( isStrReal(SK_' a') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrReal(SK_' ')")
    call disp%show( isStrReal(SK_' ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrReal(SK_'a')")
    call disp%show( isStrReal(SK_'a') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrReal(SK_'')")
    call disp%show( isStrReal(SK_'') )
    call disp%skip()

end program example