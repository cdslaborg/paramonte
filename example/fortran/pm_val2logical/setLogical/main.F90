program example

    use pm_kind, only: IK
    use pm_kind, only: SK ! all other string  kinds are also supported.
    use pm_kind, only: LK ! all other logical kinds are also supported.
    use pm_io, only: display_type
    use pm_val2logical, only: setLogical

    implicit none

    logical(LK) :: Conversion(10)
    integer(IK) :: iostat(10)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Convert string values to logical.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogical(Conversion(1), 'F')")
                    call setLogical(Conversion(1), 'F')
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogical(Conversion(1), 'f')")
                    call setLogical(Conversion(1), 'f')
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogical(Conversion(1), 'T')")
                    call setLogical(Conversion(1), 'T')
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogical(Conversion(1), 't')")
                    call setLogical(Conversion(1), 't')
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogical(Conversion(1), '.F.')")
                    call setLogical(Conversion(1), '.F.')
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogical(Conversion(1), '.T.')")
                    call setLogical(Conversion(1), '.T.')
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogical(Conversion(1), '.false.')")
                    call setLogical(Conversion(1), '.false.')
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogical(Conversion(1), '.true.')")
                    call setLogical(Conversion(1), '.true.')
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogical(Conversion(1), 'false')")
                    call setLogical(Conversion(1), 'false')
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogical(Conversion(1), 'true')")
                    call setLogical(Conversion(1), 'true')
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogical(Conversion(1), 'true, false')")
                    call setLogical(Conversion(1), 'true, false')
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogical(Conversion(1:2), [character(5,SK) :: 'true', 'false'])")
                    call setLogical(Conversion(1:2), [character(5,SK) :: 'true', 'false'])
    call disp%show("Conversion(1:2)")
    call disp%show( Conversion(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogical(Conversion(1), 'true, paramonte', iostat(1))")
                    call setLogical(Conversion(1), 'true, paramonte', iostat(1))
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%show("iostat(1)")
    call disp%show( iostat(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogical(Conversion(1), 'paramonte, true', iostat(1))")
                    call setLogical(Conversion(1), 'paramonte, true', iostat(1))
    call disp%show("Conversion(1)")
    call disp%show( Conversion(1) )
    call disp%show("iostat(1)")
    call disp%show( iostat(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogical(Conversion(1:3), [character(5,SK) :: 't', 'F', 'TRUE'])")
                    call setLogical(Conversion(1:3), [character(5,SK) :: 't', 'F', 'TRUE'])
    call disp%show("Conversion(1:3)")
    call disp%show( Conversion(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogical(Conversion(1:3), [character(10,SK) :: 'F', 'True', 'paramonte'], iostat(1:3))")
                    call setLogical(Conversion(1:3), [character(10,SK) :: 'F', 'True', 'paramonte'], iostat(1:3))
    call disp%show("Conversion(1:2) ! Do not print the last element since it was an illegal conversion.")
    call disp%show( Conversion(1:2) )
    call disp%show("iostat(1:3)")
    call disp%show( iostat(1:3) )
    call disp%skip()

end program example