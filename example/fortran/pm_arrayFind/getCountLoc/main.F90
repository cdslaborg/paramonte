program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: IKG => IK ! All kinds are supported.
    use pm_kind, only: CKG => CK ! All kinds are supported.
    use pm_kind, only: RKG => RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayFind, only: getCountLoc, discrete

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    !array = [1., 0., 2., 0., 3., 0., 4.]
    !array = [(1., -1.), (0., -0.), (2., -2.), (0., -0.), (3., -3.), (0., -0.), (4., -4.)]
    !array = [.false._LK, .true._LK, .false._LK, .true._LK, .false._LK, .true._LK, .false._LK]

    !pattern = ["XXX"]
    !pattern = [0_IK]
    !pattern = [0._RK]
    !pattern = [(0._CK, -0._CK)]
    !pattern = [.true._LK]

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!count character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: SKG => SK ! All kinds are supported.
        character(:,SKG), allocatable :: array, pattern
        call disp%show("array = 'Paramonte is a Machine Learning   Library '")
                        array = 'Paramonte is a Machine Learning   Library '
        call disp%show("pattern = ' '")
                        pattern = ' '
        call disp%show("getCountLoc(array, pattern)")
        call disp%show( getCountLoc(array, pattern) )
        call disp%show("getCountLoc(array, pattern, blindness = 3_IK)")
        call disp%show( getCountLoc(array, pattern, blindness = 3_IK) )
        call disp%show("getCountLoc(array, pattern, border = discrete)")
        call disp%show( getCountLoc(array, pattern, border = discrete) )
        call disp%show("getCountLoc(array, pattern, border = discrete, blindness = 3_IK)")
        call disp%show( getCountLoc(array, pattern, border = discrete, blindness = 3_IK) )
        call disp%show("getCountLoc(array, SKG_'m', iseq_SK) ! find with custom case-insensitive search.")
        call disp%show( getCountLoc(array, SKG_'m', iseq_SK) )
        call disp%show("getCountLoc(array, SKG_'m')")
        call disp%show( getCountLoc(array, SKG_'m') )
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!count character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: SKG => SK ! All kinds are supported.
        character(9,SKG), allocatable :: array(:), pattern
        call disp%show("array = [character(9,SKG) :: 'xxx', 'Paramonte', 'XXX', 'is', 'XXX', 'a', 'XXX', 'Monte', 'XXX', 'Carlo', 'XXX', 'XXX', 'XXX', 'Library.', 'XXX']")
                        array = [character(9,SKG) :: 'xxx', 'Paramonte', 'XXX', 'is', 'XXX', 'a', 'XXX', 'Monte', 'XXX', 'Carlo', 'XXX', 'XXX', 'XXX', 'Library.', 'XXX']
        call disp%show("pattern = 'XXX'")
                        pattern = 'XXX'
        call disp%show("getCountLoc(array, pattern)")
        call disp%show( getCountLoc(array, pattern) )
        call disp%show("getCountLoc(array, pattern, blindness = 3_IK)")
        call disp%show( getCountLoc(array, pattern, blindness = 3_IK) )
        call disp%show("getCountLoc(array, pattern, border = discrete)")
        call disp%show( getCountLoc(array, pattern, border = discrete) )
        call disp%show("getCountLoc(array, pattern, border = discrete, blindness = 3_IK)")
        call disp%show( getCountLoc(array, pattern, border = discrete, blindness = 3_IK) )
        call disp%show("getCountLoc(array, SKG_'xxx', iseq_SK) ! find with custom case-insensitive search.")
        call disp%show( getCountLoc(array, SKG_'xxx', iseq_SK) )
        call disp%show("getCountLoc(array, SKG_'xxx')")
        call disp%show( getCountLoc(array, SKG_'xxx') )
        call disp%show("getCountLoc(array, [character(3,SKG) :: 'XXX', 'XXX']) ! vector pattern")
        call disp%show( getCountLoc(array, [character(3,SKG) :: 'XXX', 'XXX']) )
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!count integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: IKG => IK ! All kinds are supported.
        integer(IKG), allocatable :: array(:), pattern
        call disp%show("array = [-1, 1, -2, 0, 0, 0, 2, 0, 3, 0, 4]")
                        array = [-1, 1, -2, 0, 0, 0, 2, 0, 3, 0, 4]
        call disp%show("pattern = 0")
                        pattern = 0
        call disp%show("getCountLoc(array, pattern)")
        call disp%show( getCountLoc(array, pattern) )
        call disp%show("getCountLoc(array, pattern, blindness = 3_IK)")
        call disp%show( getCountLoc(array, pattern, blindness = 3_IK) )
        call disp%show("getCountLoc(array, pattern, border = discrete)")
        call disp%show( getCountLoc(array, pattern, border = discrete) )
        call disp%show("getCountLoc(array, pattern, border = discrete, blindness = 3_IK)")
        call disp%show( getCountLoc(array, pattern, border = discrete, blindness = 3_IK) )
        call disp%show("pattern = 2")
                        pattern = 2
        call disp%show("getCountLoc(array, pattern, iseq_IK) ! find any pattern+-1 with custom search.")
        call disp%show( getCountLoc(array, pattern, iseq_IK) )
        call disp%show("getCountLoc(array, pattern)")
        call disp%show( getCountLoc(array, pattern) )
        call disp%show("pattern = 0")
                        pattern = 0
        call disp%show("getCountLoc(array, [pattern, pattern]) ! vector pattern")
        call disp%show( getCountLoc(array, [pattern, pattern]) )
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!count logical array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: LKG => LK ! All kinds are supported.
        logical(LKG), allocatable :: array(:), pattern
        call disp%show("array = [.false., .true., .false., .true., .false., .false., .false., .true., .false.]")
                        array = [.false., .true., .false., .true., .false., .false., .false., .true., .false.]
        call disp%show("pattern = .false.")
                        pattern = .false.
        call disp%show("getCountLoc(array, pattern)")
        call disp%show( getCountLoc(array, pattern) )
        call disp%show("getCountLoc(array, pattern, blindness = 3_IK)")
        call disp%show( getCountLoc(array, pattern, blindness = 3_IK) )
        call disp%show("getCountLoc(array, pattern, border = discrete)")
        call disp%show( getCountLoc(array, pattern, border = discrete) )
        call disp%show("getCountLoc(array, pattern, border = discrete, blindness = 3_IK)")
        call disp%show( getCountLoc(array, pattern, border = discrete, blindness = 3_IK) )
        call disp%show("pattern = .true.")
                        pattern = .true.
        call disp%show("getCountLoc(array, [pattern, .not. pattern], iseqall_LK) ! find any non-equivalent logical pair with custom search: [.true., .false.] or [.false., .true.].")
        call disp%show( getCountLoc(array, [pattern, .not. pattern], iseqall_LK) )
        call disp%show("getCountLoc(array, [pattern, .not. pattern])")
        call disp%show( getCountLoc(array, [pattern, .not. pattern]) )
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!count complex array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: CKG => RKS ! All kinds are supported.
        complex(CKG), allocatable :: array(:), pattern
        call disp%show("array = [(-1., +1.), (1., 1.), (-1., 0.), (0., 0.), (0., 0.), (0., 0.), (0., 1.), (0., 0.), (3., -3.), (0., 0.), (4., -4)]")
                        array = [(-1., +1.), (1., 1.), (-1., 0.), (0., 0.), (0., 0.), (0., 0.), (0., 1.), (0., 0.), (3., -3.), (0., 0.), (4., -4)]
        call disp%show("pattern = 0")
                        pattern = 0
        call disp%show("getCountLoc(array, pattern)")
        call disp%show( getCountLoc(array, pattern) )
        call disp%show("getCountLoc(array, pattern, blindness = 3_IK)")
        call disp%show( getCountLoc(array, pattern, blindness = 3_IK) )
        call disp%show("getCountLoc(array, pattern, border = discrete)")
        call disp%show( getCountLoc(array, pattern, border = discrete) )
        call disp%show("getCountLoc(array, pattern, border = discrete, blindness = 3_IK)")
        call disp%show( getCountLoc(array, pattern, border = discrete, blindness = 3_IK) )
        call disp%show("pattern = (0., 0.) ! dummy search value.")
                        pattern = (0., 0.)
        call disp%show("getCountLoc(array, pattern, iseq_CK) ! find any complex value whose components have opposite signs.")
        call disp%show( getCountLoc(array, pattern, iseq_CK) )
        call disp%show("getCountLoc(array, pattern)")
        call disp%show( getCountLoc(array, pattern) )
        call disp%show("pattern = 0")
                        pattern = 0
        call disp%show("getCountLoc(array, [pattern, pattern]) ! vector pattern")
        call disp%show( getCountLoc(array, [pattern, pattern]) )
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%show("!count real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: RKG => RKS ! All kinds are supported.
        real(RKG), allocatable :: array(:), pattern
        call disp%show("array = [-1, 1, -2, 0, 0, 0, 2, 0, 3, 0, 4]")
                        array = [-1, 1, -2, 0, 0, 0, 2, 0, 3, 0, 4]
        call disp%show("pattern = 0")
                        pattern = 0
        call disp%show("getCountLoc(array, pattern)")
        call disp%show( getCountLoc(array, pattern) )
        call disp%show("getCountLoc(array, pattern, blindness = 3_IK)")
        call disp%show( getCountLoc(array, pattern, blindness = 3_IK) )
        call disp%show("getCountLoc(array, pattern, border = discrete)")
        call disp%show( getCountLoc(array, pattern, border = discrete) )
        call disp%show("getCountLoc(array, pattern, border = discrete, blindness = 3_IK)")
        call disp%show( getCountLoc(array, pattern, border = discrete, blindness = 3_IK) )
        call disp%show("pattern = 2")
                        pattern = 2
        call disp%show("getCountLoc(array, pattern, iseq_RK) ! find any whose absolute value is within 0.5 of `pattern`.")
        call disp%show( getCountLoc(array, pattern, iseq_RK) )
        call disp%show("getCountLoc(array, pattern)")
        call disp%show( getCountLoc(array, pattern) )
        call disp%show("pattern = 0")
                        pattern = 0
        call disp%show("getCountLoc(array, [pattern, pattern]) ! vector pattern")
        call disp%show( getCountLoc(array, [pattern, pattern]) )
    end block

contains

    pure function iseq_SK(segment, pattern) result(equivalent)
        use pm_strASCII, only: getStrLower
        character(*, SK), intent(in)    :: segment, pattern
        logical(LK)                     :: equivalent
        equivalent = pattern == getStrLower(segment)
    end function

    function iseq_IK(segment, pattern) result(equivalent)
        use pm_kind, only: IKG => IK
        integer(IKG)    , intent(in)    :: segment, pattern
        logical(LK)                     :: equivalent
        equivalent = pattern - 2 < segment .and. segment < pattern + 2
    end function

    function iseqall_LK(segment, pattern, lenPattern) result(equivalent)
        use pm_kind, only: LKG => LK
        integer(IK)     , intent(in)    :: lenPattern
        logical(LKG)    , intent(in)    :: pattern(lenPattern), segment(lenPattern)
        logical(LK)                     :: equivalent
        equivalent = segment(1) .neqv. segment(2)
    end function

    function iseq_CK(segment, pattern) result(equivalent)
        use pm_kind, only: CKG => RKS
        complex(CKG)    , intent(in)    :: segment, pattern
        logical(LK)                     :: equivalent
        equivalent = segment%re * segment%im < 0._CKG
    end function

    function iseq_RK(segment, pattern) result(equivalent)
        use pm_kind, only: RKG => RKS
        real(RKG)       , intent(in)    :: segment, pattern
        logical(LK)                     :: equivalent
        equivalent = abs(abs(segment) - pattern) <= 0.5_RKG
    end function

end program example