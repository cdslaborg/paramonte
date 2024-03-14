program example

    use pm_kind, only: LK
    use pm_kind, only: SKC => SK ! All kinds are supported.
    use pm_kind, only: IKC => IK ! All kinds are supported.
    use pm_kind, only: LKC => LK ! All kinds are supported.
    use pm_kind, only: CKC => CK ! All kinds are supported.
    use pm_kind, only: RKC => RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayMembership, only: operator(.in.)

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        character(:,SKC), allocatable :: val, Set
        Set = "ParaMonte is a Machine Learning Library."
        val = "M"
        call disp%skip()
        call disp%show("val")
        call disp%show( val, deliml = SKC_"""" )
        call disp%show("Set")
        call disp%show( Set, deliml = SKC_"""" )
        call disp%show("val .in. Set")
        call disp%show( val .in. Set )
        call disp%skip()
    end block

    block
        character(:,SKC), allocatable :: val, Set
        Set = "ParaMonte is a Machine Learning Library."
        val = "paramonte"
        call disp%skip()
        call disp%show("val")
        call disp%show( val, deliml = SKC_"""" )
        call disp%show("Set")
        call disp%show( Set, deliml = SKC_"""" )
        call disp%show("val .in. Set")
        call disp%show( val .in. Set )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        character(9,SKC) :: val
        character(10,SKC), allocatable :: Set(:)
        Set = [character(10,SKC) :: "ParaMonte", "is       ", "a        ", "Monte    ", "Carlo    ", "Library. "]
        val = "paramonte"
        call disp%skip()
        call disp%show("val")
        call disp%show( val, deliml = SKC_"""" )
        call disp%show("Set")
        call disp%show( Set, deliml = SKC_"""" )
        call disp%show("val .in. Set")
        call disp%show( val .in. Set )
        call disp%skip()
    end block

    block
        character(9,SKC), allocatable :: val(:)
        character(10,SKC), allocatable :: Set(:)
        Set = [character(10,SKC) :: "ParaMonte", "is       ", "a        ", "Monte    ", "Carlo    ", "Library. "]
        val = [character( 9,SKC) :: "paramonte", "ParaMonte", "Carlo", "carlo"]
        call disp%skip()
        call disp%show("val")
        call disp%show( val, deliml = SKC_"""" )
        call disp%show("Set")
        call disp%show( Set, deliml = SKC_"""" )
        call disp%show("val .in. Set")
        call disp%show( val .in. Set )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        integer(IKC), allocatable :: val, Set(:)
        call disp%skip()
        call disp%show("val = 1_IKC")
                        val = 1_IKC
        call disp%show("Set = [integer(IKC) :: 0, 1, 2, 3, 4]")
                        Set = [integer(IKC) :: 0, 1, 2, 3, 4]
        call disp%show("val .in. Set")
        call disp%show( val .in. Set )
        call disp%skip()
    end block

    block
        integer(IKC), allocatable :: val(:), Set(:)
        call disp%skip()
        call disp%show("val = [integer(IKC) :: -1, 3, 5]")
                        val = [integer(IKC) :: -1, 3, 5]
        call disp%show("Set = [integer(IKC) :: 0, 1, 2, 3, 4]")
                        Set = [integer(IKC) :: 0, 1, 2, 3, 4]
        call disp%show("val .in. Set")
        call disp%show( val .in. Set )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find logical array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        logical(LKC), allocatable :: val, Set(:)
        call disp%skip()
        call disp%show("val = .false.")
                        val = .false.
        call disp%show("Set = [.true., .true., .true., .true.]")
                        Set = [.true., .true., .true., .true.]
        call disp%show("val .in. Set")
        call disp%show( val .in. Set )
        call disp%skip()
    end block

    block
        logical(LKC), allocatable :: val(:), Set(:)
        call disp%skip()
        call disp%show("val = [.false., .true., .false.]")
                        val = [.false., .true., .false.]
        call disp%show("Set = [.true., .true., .true., .true.]")
                        Set = [.true., .true., .true., .true.]
        call disp%show("val .in. Set")
        call disp%show( val .in. Set )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find complex array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        integer(IKC), allocatable :: val, Set(:)
        call disp%skip()
        call disp%show("val = (1., -1.)")
                        val = (1., -1.)
        call disp%show("Set = [complex(CKC) :: (1., 0.), -(1., 0.), (1., -1.)]")
                        Set = [complex(CKC) :: (1., 0.), -(1., 0.), (1., -1.)]
        call disp%show("val .in. Set")
        call disp%show( val .in. Set )
        call disp%skip()
    end block

    block
        integer(IKC), allocatable :: val(:), Set(:)
        call disp%skip()
        call disp%show("val = [(1., -1.), (0., -1.)]")
                        val = [(1., -1.), (0., -1.)]
        call disp%show("Set = [complex(CKC) :: (1., 0.), -(1., 0.), (1., -1.)]")
                        Set = [complex(CKC) :: (1., 0.), -(1., 0.), (1., -1.)]
        call disp%show("val .in. Set")
        call disp%show( val .in. Set )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%show("! Find real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(RKC), allocatable :: val, Set(:)
        call disp%skip()
        call disp%show("val = 1._RKC")
                        val = 1._RKC
        call disp%show("Set = [real(RKC) :: 0, 1, 2, 3, 4]")
                        Set = [real(RKC) :: 0, 1, 2, 3, 4]
        call disp%show("val .in. Set")
        call disp%show( val .in. Set )
        call disp%skip()
    end block

    block
        real(RKC), allocatable :: val(:), Set(:)
        call disp%skip()
        call disp%show("val = [real(RKC) :: -1, 1, 5]")
                        val = [real(RKC) :: -1, 1, 5]
        call disp%show("Set = [real(RKC) :: 0, 1, 2, 3, 4]")
                        Set = [real(RKC) :: 0, 1, 2, 3, 4]
        call disp%show("val .in. Set")
        call disp%show( val .in. Set )
        call disp%skip()
    end block

end program example