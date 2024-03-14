program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: isStrNumber

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isStrNumber(SK_'(0,0)')")
    call disp%show( isStrNumber(SK_'(0,0)') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrNumber(SK_'(+1,-1)')")
    call disp%show( isStrNumber(SK_'(+1,-1)') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrNumber(SK_' +1.d123)')")
    call disp%show( isStrNumber(SK_' +1.d123)') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrNumber(SK_' -5.e0 )')")
    call disp%show( isStrNumber(SK_' -5.e0 )') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrNumber(SK_' +1 ')")
    call disp%show( isStrNumber(SK_' +1 ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrNumber(SK_'-1 ')")
    call disp%show( isStrNumber(SK_'-1 ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrNumber([character(15,SK) :: '-.1', '12345', '(-1.D0,0)'])")
    call disp%show( isStrNumber([character(15,SK) :: '-.1', '12345', '(-1.D0,0)']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrNumber([character(15,SK) :: '1 2 3', 'paramonte', '.1.1'])")
    call disp%show( isStrNumber([character(15,SK) :: '1 2 3', 'paramonte', '.1.1']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrNumber(SK_'-')")
    call disp%show( isStrNumber(SK_'-') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrNumber(SK_'+')")
    call disp%show( isStrNumber(SK_'+') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrNumber(SK_' ')")
    call disp%show( isStrNumber(SK_' ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrNumber(SK_'')")
    call disp%show( isStrNumber(SK_'') )
    call disp%skip()

end program example