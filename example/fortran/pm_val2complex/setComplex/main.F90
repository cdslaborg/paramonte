program example

    use pm_kind, only: IK
    use pm_kind, only: CK => CK32, LK, SK ! all processor types and kinds are supported.
    use pm_io, only: display_type
    use pm_val2complex, only: setComplex

    implicit none

    integer(IK) :: iostat(10)
    complex(CK) :: Conversion(10)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Convert logical values to complex.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setComplex(Conversion(1), .true._LK)")
                    call setComplex(Conversion(1), .true._LK)
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setComplex(Conversion(1), .false._LK)")
                    call setComplex(Conversion(1), .false._LK)
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setComplex(Conversion(1:3), [.true._LK, .false._LK, .true._LK])")
                    call setComplex(Conversion(1:3), [.true._LK, .false._LK, .true._LK])
    call disp%show("Conversion(1:3)")
    call disp%show( Conversion(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Convert string values to complex.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setComplex(Conversion(1), SK_'1.e0')")
                    call setComplex(Conversion(1), SK_'1.e0')
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setComplex(Conversion(1), SK_'1.d0')")
                    call setComplex(Conversion(1), SK_'1.d0')
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setComplex(Conversion(1:2), SK_'1.d0 -1.d0')")
                    call setComplex(Conversion(1:2), SK_'1.d0 -1.d0')
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setComplex(Conversion(1), SK_'1.d0, -1.d0')")
                    call setComplex(Conversion(1), SK_'1.d0, -1.d0')
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setComplex(Conversion(1), SK_'(1.d0, -1.d0)')")
                    call setComplex(Conversion(1), SK_'(1.d0, -1.d0)')
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setComplex(Conversion(1:3), [character(10,SK) :: '(1, -1)', '10', '100, -100'])")
                    call setComplex(Conversion(1:3), [character(10,SK) :: '(1, -1)', '10', '100, -100'])
    call disp%show("Conversion(1:3)")
    call disp%show( Conversion(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setComplex(Conversion(1), SK_'(1.d0, -1.e0)', iostat(1))")
                    call setComplex(Conversion(1), SK_'(1.d0, -1.e0)', iostat(1))
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%show("iostat(1)")
    call disp%show( iostat(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setComplex(Conversion(1), SK_'(1., -1.), paramonte', iostat(1))")
                    call setComplex(Conversion(1), SK_'(1., -1.), paramonte', iostat(1))
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%show("iostat(1)")
    call disp%show( iostat(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setComplex(Conversion(1), SK_'1., -1., paramonte', iostat(1))")
                    call setComplex(Conversion(1), SK_'1., -1., paramonte', iostat(1))
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%show("iostat(1)")
    call disp%show( iostat(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setComplex(Conversion(1), SK_'1., paramonte', iostat(1))")
                    call setComplex(Conversion(1), SK_'1., paramonte', iostat(1))
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%show("iostat(1)")
    call disp%show( iostat(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setComplex(Conversion(1), SK_'1.d0', iostat(1))")
                    call setComplex(Conversion(1), SK_'1.d0', iostat(1))
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%show("iostat(1)")
    call disp%show( iostat(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setComplex(Conversion(1:4), [character(10,SK) :: '1, -1', '1, paramonte', '100, -100', 'paramonte, 1'], iostat(1:4))")
                    call setComplex(Conversion(1:4), [character(10,SK) :: '1, -1', '1, paramonte', '100, -100', 'paramonte, 1'], iostat(1:4))
    call disp%show("Conversion(1:3)")
    call disp%show( Conversion(1:3) )
    call disp%show("iostat(1:4)")
    call disp%show( iostat(1:4) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setComplex(Conversion(1), SK_'paramonte', iostat(1))")
                    call setComplex(Conversion(1), SK_'paramonte', iostat(1))
    call disp%show("iostat(1)")
    call disp%show( iostat(1) )
    call disp%skip()

end program example