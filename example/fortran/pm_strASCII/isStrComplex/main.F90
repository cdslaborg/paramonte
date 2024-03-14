program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: isStrComplex

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isStrComplex(SK_'(0,0)')")
    call disp%show( isStrComplex(SK_'(0,0)') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrComplex(SK_'(+1,-1)')")
    call disp%show( isStrComplex(SK_'(+1,-1)') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrComplex(SK_'( +1,-1)')")
    call disp%show( isStrComplex(SK_'( +1,-1)') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrComplex(SK_'( +1 ,-1)')")
    call disp%show( isStrComplex(SK_'( +1 ,-1)') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrComplex(SK_'( +1 , -1)')")
    call disp%show( isStrComplex(SK_'( +1 , -1)') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrComplex(SK_'( +1 , -1 )')")
    call disp%show( isStrComplex(SK_'( +1 , -1 )') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrComplex(SK_'  ( +1 , -1 )')")
    call disp%show( isStrComplex(SK_'  ( +1 , -1 )') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrComplex(SK_'  ( +1 , -1  )  ')")
    call disp%show( isStrComplex(SK_'  ( +1 , -1  )  ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrComplex(SK_'( +1.d-0 , -.1E+5 )')")
    call disp%show( isStrComplex(SK_'( +1.d-0 , -.1E+5 )') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrComplex([character(15,SK) :: '(1.0,-.1)', '(1e0,+1.)', '(-1.D0,0)'])")
    call disp%show( isStrComplex([character(15,SK) :: '(1.0,-.1)', '(1e0,+1.)', '(-1.D0,0)']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrComplex([character(10,SK) :: '(1.0,-0.1)', '(1e0,+1.0)', '(-1.D0,0.)'])")
    call disp%show( isStrComplex([character(10,SK) :: '(1.0,-0.1)', '(1e0,+1.0)', '(-1.D0,0.)']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrComplex([character(5,SK) :: '(1.0,-.1)', '(1e0,+1.)', '(-1.D0,0)'])")
    call disp%show( isStrComplex([character(5,SK) :: '(1.0,-.1)', '(1e0,+1.)', '(-1.D0,0)']) )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrComplex(SK_'  ( + 1 , - 1 )  ')")
    call disp%show( isStrComplex(SK_'  ( + 1 , - 1 )  ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrComplex(SK_'  ( + 1 , - 1 )')")
    call disp%show( isStrComplex(SK_'  ( + 1 , - 1 )') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrComplex(SK_'( + 1 , - 1 )')")
    call disp%show( isStrComplex(SK_'( + 1 , - 1 )') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrComplex(SK_'-.0')")
    call disp%show( isStrComplex(SK_'-.0') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrComplex(SK_' 1')")
    call disp%show( isStrComplex(SK_' 1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrComplex(SK_'.2')")
    call disp%show( isStrComplex(SK_'.2') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrComplex(SK_'-1')")
    call disp%show( isStrComplex(SK_'-1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrComplex(SK_'+1')")
    call disp%show( isStrComplex(SK_'+1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrComplex(SK_'1 ')")
    call disp%show( isStrComplex(SK_'1 ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrComplex(SK_'1')")
    call disp%show( isStrComplex(SK_'1') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrComplex(SK_' ')")
    call disp%show( isStrComplex(SK_' ') )
    call disp%skip()

    call disp%skip()
    call disp%show("isStrComplex(SK_'')")
    call disp%show( isStrComplex(SK_'') )
    call disp%skip()

end program example