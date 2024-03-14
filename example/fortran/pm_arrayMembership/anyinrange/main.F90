program example

    use pm_kind, only: LK
    use pm_kind, only: SKC => SK ! All kinds are supported.
    use pm_kind, only: IKC => IK ! All kinds are supported.
    use pm_kind, only: LKC => LK ! All kinds are supported.
    use pm_kind, only: CKC => CK ! All kinds are supported.
    use pm_kind, only: RKC => RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayMembership, only: operator(.anyinrange.)

    implicit none

    type(display_type) :: disp
    logical(LK), allocatable :: Member(:)
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Check character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        character(:,SKC), allocatable :: val, Set
        Set = "AZ"
        val = "M"
        call disp%skip()
        call disp%show("val")
        call disp%show( val, deliml = SKC_"""" )
        call disp%show("Set")
        call disp%show( Set, deliml = SKC_"""" )
        call disp%show("val .anyinrange. Set")
        call disp%show( val .anyinrange. Set )
        call disp%skip()
    end block

    block
        character(:,SKC), allocatable :: val, Set
        Set = "Az"
        val = "ParaMonte."
        call disp%skip()
        call disp%show("val")
        call disp%show( val, deliml = SKC_"""" )
        call disp%show("Set")
        call disp%show( Set, deliml = SKC_"""" )
        call disp%show("val .anyinrange. Set")
        call disp%show( val .anyinrange. Set )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Check character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        character(:,SKC), allocatable :: val(:)
        character(2,SKC), allocatable :: Set(:)
        Set = [character(2,SKC) :: "aa", "zz"]
        val = ["paramonte"]
        call disp%skip()
        call disp%show("val")
        call disp%show( val, deliml = SKC_"""" )
        call disp%show("Set")
        call disp%show( Set, deliml = SKC_"""" )
        call disp%show("val .anyinrange. Set")
        call disp%show( val .anyinrange. Set )
        call disp%skip()
    end block

    block
        character( 9,SKC), allocatable :: val(:)
        character(10,SKC), allocatable :: Set(:)
        Set = [character(10,SKC) :: "aa", "zz"]
        val = [character( 9,SKC) :: "paramonte", "ParaMonte", "Carlo", "carlo"]
        call disp%skip()
        call disp%show("val")
        call disp%show( val, deliml = SKC_"""" )
        call disp%show("Set")
        call disp%show( Set, deliml = SKC_"""" )
        call disp%show("val .anyinrange. Set")
        call disp%show( val .anyinrange. Set )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Check integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        integer(IKC), allocatable :: val(:), Set(:)
        call disp%skip()
        call disp%show("val = [1_IKC]")
                        val = [1_IKC]
        call disp%show("Set = [integer(IKC) :: 1, 3]")
                        Set = [integer(IKC) :: 1, 3]
        call disp%show("val .anyinrange. Set")
        call disp%show( val .anyinrange. Set )
        call disp%skip()
    end block

    block
        integer(IKC), allocatable :: val(:), Set(:)
        call disp%skip()
        call disp%show("val = [integer(IKC) :: -1, 3, 5]")
                        val = [integer(IKC) :: -1, 3, 5]
        call disp%show("Set = [integer(IKC) :: 1, 3]")
                        Set = [integer(IKC) :: 1, 3]
        call disp%show("val .anyinrange. Set")
        call disp%show( val .anyinrange. Set )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Check logical array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        logical(LKC), allocatable :: val(:), Set(:)
        call disp%skip()
        call disp%show("val = [.false.]")
                        val = [.false.]
        call disp%show("Set = [.false., .true.]")
                        Set = [.false., .true.]
        call disp%show("val .anyinrange. Set")
        call disp%show( val .anyinrange. Set )
        call disp%skip()
    end block

    block
        logical(LKC), allocatable :: val(:), Set(:)
        call disp%skip()
        call disp%show("val = [.false., .true., .false.]")
                        val = [.false., .true., .false.]
        call disp%show("Set = [.true., .true.]")
                        Set = [.true., .true.]
        call disp%show("val .anyinrange. Set")
        call disp%show( val .anyinrange. Set )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Check complex array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        integer(IKC), allocatable :: val(:), Set(:)
        call disp%skip()
        call disp%show("val = [(1., -1.)]")
                        val = [(1., -1.)]
        call disp%show("Set = [complex(CKC) :: (1., 0.), (1., -1.)]")
                        Set = [complex(CKC) :: (1., 0.), (1., -1.)]
        call disp%show("val .anyinrange. Set")
        call disp%show( val .anyinrange. Set )
        call disp%skip()
    end block

    block
        integer(IKC), allocatable :: val(:), Set(:)
        call disp%skip()
        call disp%show("val = [(1., -1.), (0., -1.), (2., -1.)]")
                        val = [(1., -1.), (0., -1.), (2., -1.)]
        call disp%show("Set = [complex(CKC) :: (1., 0.), (1., -1.)]")
                        Set = [complex(CKC) :: (1., 0.), (1., -1.)]
        call disp%show("val .anyinrange. Set")
        call disp%show( val .anyinrange. Set )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%show("!Check real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(RKC), allocatable :: val(:), Set(:)
        call disp%skip()
        call disp%show("val = [1._RKC]")
                        val = [1._RKC]
        call disp%show("Set = [real(RKC) :: 1, 3]")
                        Set = [real(RKC) :: 1, 3]
        call disp%show("val .anyinrange. Set")
        call disp%show( val .anyinrange. Set )
        call disp%skip()
    end block

    block
        real(RKC), allocatable :: val(:), Set(:)
        call disp%skip()
        call disp%show("val = [real(RKC) :: -1, 1, 5]")
                        val = [real(RKC) :: -1, 1, 5]
        call disp%show("Set = [real(RKC) :: 1, 3]")
                        Set = [real(RKC) :: 1, 3]
        call disp%show("val .anyinrange. Set")
        call disp%show( val .anyinrange. Set )
        call disp%skip()
    end block

end program example