program example

    use pm_kind, only: LK, SK
    use pm_io, only: display_type
    use pm_arrayRange, only: setRange

    implicit none


    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("!%%%%%%%%%%")
    call disp%show("!real range")
    call disp%show("!%%%%%%%%%%")
    call disp%skip

    block
        real :: range(10) ! all non-default kinds are also supported.

        call disp%skip
        call disp%show("size(range)")
        call disp%show( size(range) )
        call disp%show("call setRange(range, 0.)")
                        call setRange(range, 0.)
        call disp%show("range")
        call disp%show( range )
        call disp%skip

        call disp%skip
        call disp%show("size(range)")
        call disp%show( size(range) )
        call disp%show("call setRange(range, 1.)")
                        call setRange(range, 1.)
        call disp%show("range")
        call disp%show( range )
        call disp%skip

        call disp%skip
        call disp%show("size(range)")
        call disp%show( size(range) )
        call disp%show("call setRange(range, -huge(1.))")
                        call setRange(range, -huge(1.))
        call disp%show("range")
        call disp%show( range )
        call disp%skip

        call disp%skip
        call disp%show("size(range)")
        call disp%show( size(range) )
        call disp%show("call setRange(range, 0., -2.5)")
                        call setRange(range, 0., -2.5)
        call disp%show("range")
        call disp%show( range )
        call disp%skip

        call disp%skip
        call disp%show("size(range)")
        call disp%show( size(range) )
        call disp%show("call setRange(range, 0., +2.5)")
                        call setRange(range, 0., +2.5)
        call disp%show("range")
        call disp%show( range )
        call disp%skip

    end block

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%show("!character range")
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%skip

    block
        character(5) :: range ! all non-default kinds are also supported.
        call disp%skip
        call disp%show("call setRange(range, 'A')")
                        call setRange(range, 'A')
        call disp%show("range")
        call disp%show( range )
        call disp%skip
        
        call disp%skip
        call disp%show("call setRange(range, 'A', 2)")
                        call setRange(range, 'A', 2)
        call disp%show("range")
        call disp%show( range )
        call disp%skip

        call disp%skip
        call disp%show("call setRange(range, 'A', 3)")
                        call setRange(range, 'A', 3)
        call disp%show("range")
        call disp%show( range )
        call disp%skip

        call disp%skip
        call disp%show("call setRange(range, 'A', 3)")
                        call setRange(range, 'A', 3)
        call disp%show("range")
        call disp%show( range )
        call disp%skip

        call disp%skip
        call disp%show("call setRange(range, 'Z', -2)")
                        call setRange(range, 'Z', -2)
        call disp%show("range")
        call disp%show( range )
        call disp%skip
    end block

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%")
    call disp%show("!integer range")
    call disp%show("!%%%%%%%%%%%%%")
    call disp%skip

    block
        integer :: range(5) ! all non-default kinds are also supported.
        call disp%skip
        call disp%show("call setRange(range, 0)")
                        call setRange(range, 0)
        call disp%show("range")
        call disp%show( range )
        call disp%skip
        
        call disp%skip
        call disp%show("call setRange(range, 0, 2)")
                        call setRange(range, 0, 2)
        call disp%show("range")
        call disp%show( range )
        call disp%skip
    
        call disp%skip
        call disp%show("call setRange(range, 0, 3)")
                        call setRange(range, 0, 3)
        call disp%show("range")
        call disp%show( range )
        call disp%skip
    
        call disp%skip
        call disp%show("call setRange(range, 0, 3)")
                        call setRange(range, 0, 3)
        call disp%show("range")
        call disp%show( range )
        call disp%skip
    
        call disp%skip
        call disp%show("call setRange(range, 0, -2)")
                        call setRange(range, 0, -2)
        call disp%show("range")
        call disp%show( range )
        call disp%skip
    end block

end program example