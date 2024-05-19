program example

    use pm_kind, only: LK
    use pm_kind, only: SKG => SK ! All kinds are supported.
    use pm_kind, only: IKG => IK ! All kinds are supported.
    use pm_kind, only: LKG => LK ! All kinds are supported.
    use pm_kind, only: CKG => CK ! All kinds are supported.
    use pm_kind, only: RKG => RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayMembership, only: operator(.allin.)

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Check character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        character(:,SKG), allocatable :: val, Set
        Set = "ParaMonte is a Machine Learning Library."
        val = "M"
        call disp%skip()
        call disp%show("val")
        call disp%show( val(:), deliml = SKG_"""" )
        call disp%show("Set")
        call disp%show( Set, deliml = SKG_"""" )
        call disp%show("val .allin. Set")
        call disp%show( val .allin. Set )
        call disp%skip()
    end block

    block
        character(:,SKG), allocatable :: val, Set
        Set = "ParaMonte is a Machine Learning Library."
        val = "paramonte"
        call disp%skip()
        call disp%show("val")
        call disp%show( val(:), deliml = SKG_"""" )
        call disp%show("Set")
        call disp%show( Set, deliml = SKG_"""" )
        call disp%show("val .allin. Set")
        call disp%show( val .allin. Set )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Check character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        character(9,SKG), allocatable :: val(:)
        character(10,SKG), allocatable :: Set(:)
        Set = [character(10,SKG) :: "ParaMonte", "is       ", "a        ", "Monte    ", "Carlo    ", "Library. "]
        val = ["paramonte"]
        call disp%skip()
        call disp%show("val")
        call disp%show( val(:), deliml = SKG_"""" )
        call disp%show("Set")
        call disp%show( Set, deliml = SKG_"""" )
        call disp%show("val .allin. Set")
        call disp%show( val .allin. Set )
        call disp%skip()
    end block

    block
        character(9,SKG), allocatable :: val(:)
        character(10,SKG), allocatable :: Set(:)
        Set = [character(10,SKG) :: "ParaMonte", "is       ", "a        ", "Monte    ", "Carlo    ", "Library. "]
        val = [character( 9,SKG) :: "paramonte", "ParaMonte", "Carlo", "carlo"]
        call disp%skip()
        call disp%show("val")
        call disp%show( val(:), deliml = SKG_"""" )
        call disp%show("Set")
        call disp%show( Set, deliml = SKG_"""" )
        call disp%show("val .allin. Set")
        call disp%show( val .allin. Set )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Check integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        integer(IKG), allocatable :: val(:), Set(:)
        call disp%skip()
        call disp%show("val = [1_IKG]")
                        val = [1_IKG]
        call disp%show("Set = [integer(IKG) :: 0, 1, 2, 3, 4]")
                        Set = [integer(IKG) :: 0, 1, 2, 3, 4]
        call disp%show("val .allin. Set")
        call disp%show( val .allin. Set )
        call disp%skip()
    end block

    block
        integer(IKG), allocatable :: val(:), Set(:)
        call disp%skip()
        call disp%show("val = [integer(IKG) :: -1, 3, 5]")
                        val = [integer(IKG) :: -1, 3, 5]
        call disp%show("Set = [integer(IKG) :: 0, 1, 2, 3, 4]")
                        Set = [integer(IKG) :: 0, 1, 2, 3, 4]
        call disp%show("val .allin. Set")
        call disp%show( val .allin. Set )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Check logical array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        logical(LKG), allocatable :: val(:), Set(:)
        call disp%skip()
        call disp%show("val = [.false.]")
                        val = [.false.]
        call disp%show("Set = [.true., .true., .true., .true.]")
                        Set = [.true., .true., .true., .true.]
        call disp%show("val .allin. Set")
        call disp%show( val .allin. Set )
        call disp%skip()
    end block

    block
        logical(LKG), allocatable :: val(:), Set(:)
        call disp%skip()
        call disp%show("val = [.false., .true., .false.]")
                        val = [.false., .true., .false.]
        call disp%show("Set = [.true., .true., .true., .true.]")
                        Set = [.true., .true., .true., .true.]
        call disp%show("val .allin. Set")
        call disp%show( val .allin. Set )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Check complex array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        integer(IKG), allocatable :: val(:), Set(:)
        call disp%skip()
        call disp%show("val = [(1., -1.)]")
                        val = [(1., -1.)]
        call disp%show("Set = [complex(CKG) :: (1., 0.), -(1., 0.), (1., -1.)]")
                        Set = [complex(CKG) :: (1., 0.), -(1., 0.), (1., -1.)]
        call disp%show("val .allin. Set")
        call disp%show( val .allin. Set )
        call disp%skip()
    end block

    block
        integer(IKG), allocatable :: val(:), Set(:)
        call disp%skip()
        call disp%show("val = [(1., -1.), (0., -1.)]")
                        val = [(1., -1.), (0., -1.)]
        call disp%show("Set = [complex(CKG) :: (1., 0.), -(1., 0.), (1., -1.)]")
                        Set = [complex(CKG) :: (1., 0.), -(1., 0.), (1., -1.)]
        call disp%show("val .allin. Set")
        call disp%show( val .allin. Set )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%show("!Check real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(RKG), allocatable :: val(:), Set(:)
        call disp%skip()
        call disp%show("val = [1._RKG]")
                        val = [1._RKG]
        call disp%show("Set = [real(RKG) :: 0, 1, 2, 3, 4]")
                        Set = [real(RKG) :: 0, 1, 2, 3, 4]
        call disp%show("val .allin. Set")
        call disp%show( val .allin. Set )
        call disp%skip()
    end block

    block
        real(RKG), allocatable :: val(:), Set(:)
        call disp%skip()
        call disp%show("val = [real(RKG) :: -1, 1, 5]")
                        val = [real(RKG) :: -1, 1, 5]
        call disp%show("Set = [real(RKG) :: 0, 1, 2, 3, 4]")
                        Set = [real(RKG) :: 0, 1, 2, 3, 4]
        call disp%show("val .allin. Set")
        call disp%show( val .allin. Set )
        call disp%skip()
    end block

end program example