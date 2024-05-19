program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! Any kinds are supported.
    use pm_kind, only: IK ! Any kinds are supported.
    use pm_kind, only: CK ! Any kinds are supported.
    use pm_kind, only: RKG => RK ! Any kinds are supported.
    use pm_io, only: display_type
    use pm_arrayUnique, only: isUnique

    implicit none

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("isUnique('')")
    call disp%show( isUnique('') )

    call disp%show("isUnique('a')")
    call disp%show( isUnique('a') )

    call disp%show("isUnique('abcd')")
    call disp%show( isUnique('abcd') )

    call disp%show("isUnique('abbd')")
    call disp%show( isUnique('abbd') )

    call disp%show("isUnique('aaa ')")
    call disp%show( isUnique('aaa ') )

    call disp%show("isUnique('abc ')")
    call disp%show( isUnique('abc ') )

    call disp%show("isUnique('aaa')")
    call disp%show( isUnique('aaa') )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("isUnique([character(1)::])")
    call disp%show( isUnique([character(1)::]) )

    call disp%show("isUnique(['a'])")
    call disp%show( isUnique(['a']) )

    call disp%show("isUnique(['a', 'b', 'c', 'd', ' '])")
    call disp%show( isUnique(['a', 'b', 'c', 'd', ' ']) )

    call disp%show("isUnique(['a', 'b', 'b', 'd', ' '])")
    call disp%show( isUnique(['a', 'b', 'b', 'd', ' ']) )

    call disp%show("isUnique(['a', 'a', 'a', 'a', ' '])")
    call disp%show( isUnique(['a', 'a', 'a', 'a', ' ']) )

    call disp%show("isUnique([character(2) :: 'a', 'b', 'c', 'd', 'a '])")
    call disp%show( isUnique([character(2) :: 'a', 'b', 'c', 'd', 'a ']) )

    call disp%show("isUnique(['a', 'a', 'a', 'a'])")
    call disp%show( isUnique(['a', 'a', 'a', 'a']) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("isUnique([integer::])")
    call disp%show( isUnique([integer::]) )

    call disp%show("isUnique([1])")
    call disp%show( isUnique([1]) )

    call disp%show("isUnique([1, 2, 3, 4, 0])")
    call disp%show( isUnique([1, 2, 3, 4, 0]) )

    call disp%show("isUnique([1, 2, 2, 4, 0])")
    call disp%show( isUnique([1, 2, 2, 4, 0]) )

    call disp%show("isUnique([1, 1, 1, 1, 0])")
    call disp%show( isUnique([1, 1, 1, 1, 0]) )

    call disp%show("isUnique([1, 1, 1, 1])")
    call disp%show( isUnique([1, 1, 1, 1]) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in logical array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("isUnique([logical::])")
    call disp%show( isUnique([logical::]) )

    call disp%show("isUnique([.true.])")
    call disp%show( isUnique([.true.]) )

    call disp%show("isUnique([.true., .false., .true., .false.])")
    call disp%show( isUnique([.true., .false., .true., .false.]) )

    call disp%show("isUnique([.false., .false., .false., .true.])")
    call disp%show( isUnique([.false., .false., .false., .true.]) )

    call disp%show("isUnique([.true., .true., .true., .false.])")
    call disp%show( isUnique([.true., .true., .true., .false.]) )

    call disp%show("isUnique([.false., .false., .false.])")
    call disp%show( isUnique([.false., .false., .false.]) )

    call disp%show("isUnique([.true., .true., .true.])")
    call disp%show( isUnique([.true., .true., .true.]) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in complex array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("isUnique([complex::])")
    call disp%show( isUnique([complex::]) )

    call disp%show("isUnique([(1., -1)])")
    call disp%show( isUnique([(1., -1)]) )

    call disp%show("isUnique([(1., -1.), (2., -2.), (3., -3.), (4., -4.), (0., 0.)])")
    call disp%show( isUnique([(1., -1.), (2., -2.), (3., -3.), (4., -4.), (0., 0.)]) )

    call disp%show("isUnique([(1., -1.), (2., -2.), (2., -2.), (4., -4.), (0., 0.)])")
    call disp%show( isUnique([(1., -1.), (2., -2.), (2., -2.), (4., -4.), (0., 0.)]) )

    call disp%show("isUnique([(1., -1.), (1., -1.), (1., -1.), (1., -1.), (0., 0.)])")
    call disp%show( isUnique([(1., -1.), (1., -1.), (1., -1.), (1., -1.), (0., 0.)]) )

    call disp%show("isUnique([(1., -1.), (1., -1.), (1., -1.), (1., -1.)])")
    call disp%show( isUnique([(1., -1.), (1., -1.), (1., -1.), (1., -1.)]) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("isUnique([real::])")
    call disp%show( isUnique([real::]) )

    call disp%show("isUnique([1.])")
    call disp%show( isUnique([1.]) )

    call disp%show("isUnique([1., 2., 3., 4., 0.])")
    call disp%show( isUnique([1., 2., 3., 4., 0.]) )

    call disp%show("isUnique([1., 2., 2., 4., 0.])")
    call disp%show( isUnique([1., 2., 2., 4., 0.]) )

    call disp%show("isUnique([1., 1., 1., 1., 0.])")
    call disp%show( isUnique([1., 1., 1., 1., 0.]) )

    call disp%show("isUnique([1., 1., 1., 1.])")
    call disp%show( isUnique([1., 1., 1., 1.]) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("isUnique([real(RKG) :: 1.01, 1.04, 0.98, 1.0, 1.02, 2.], iseq = iseq_RK)")
    call disp%show( isUnique([real(RKG) :: 1.01, 1.04, 0.98, 1.0, 1.02, 2.], iseq = iseq_RK) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique case-insensitive instances within the character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("isUnique('ABBAbbA', iseq = iseq_SK)")
    call disp%show( isUnique('ABBAbbA', iseq = iseq_SK), deliml = SK_"""" )

contains

    pure function iseq_RK(element1, element2) result(iseq)
        real(RKG)   , intent(in)    :: element1, element2
        logical(LK)                 :: iseq
        iseq = abs(abs(element1) - abs(element2)) < 0.05_RKG
    end function

    pure function iseq_SK(element1, element2) result(iseq)
        use pm_strASCII, only: getStrLower
        character(*, SK)    , intent(in)    :: element1, element2
        logical(LK)                         :: iseq
        iseq = getStrLower(element1) == getStrLower(element2)
    end function

end program example