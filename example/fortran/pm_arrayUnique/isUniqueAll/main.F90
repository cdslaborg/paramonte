program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! Any kinds are supported.
    use pm_kind, only: IK ! Any kinds are supported.
    use pm_kind, only: CK ! Any kinds are supported.
    use pm_kind, only: RKG => RK ! Any kinds are supported.
    use pm_io, only: display_type
    use pm_arrayUnique, only: isUniqueAll

    implicit none

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("isUniqueAll('')")
    call disp%show( isUniqueAll('') )

    call disp%show("isUniqueAll('a')")
    call disp%show( isUniqueAll('a') )

    call disp%show("isUniqueAll('abcd')")
    call disp%show( isUniqueAll('abcd') )

    call disp%show("isUniqueAll('abbd')")
    call disp%show( isUniqueAll('abbd') )

    call disp%show("isUniqueAll('aaa ')")
    call disp%show( isUniqueAll('aaa ') )

    call disp%show("isUniqueAll('abc ')")
    call disp%show( isUniqueAll('abc ') )

    call disp%show("isUniqueAll('aaa')")
    call disp%show( isUniqueAll('aaa') )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("isUniqueAll([character(1)::])")
    call disp%show( isUniqueAll([character(1)::]) )

    call disp%show("isUniqueAll(['a'])")
    call disp%show( isUniqueAll(['a']) )

    call disp%show("isUniqueAll(['a', 'b', 'c', 'd', ' '])")
    call disp%show( isUniqueAll(['a', 'b', 'c', 'd', ' ']) )

    call disp%show("isUniqueAll(['a', 'b', 'b', 'd', ' '])")
    call disp%show( isUniqueAll(['a', 'b', 'b', 'd', ' ']) )

    call disp%show("isUniqueAll(['a', 'a', 'a', 'a', ' '])")
    call disp%show( isUniqueAll(['a', 'a', 'a', 'a', ' ']) )

    call disp%show("isUniqueAll([character(2) :: 'a', 'b', 'c', 'd', 'a '])")
    call disp%show( isUniqueAll([character(2) :: 'a', 'b', 'c', 'd', 'a ']) )

    call disp%show("isUniqueAll(['a', 'a', 'a', 'a'])")
    call disp%show( isUniqueAll(['a', 'a', 'a', 'a']) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("isUniqueAll([integer::])")
    call disp%show( isUniqueAll([integer::]) )

    call disp%show("isUniqueAll([1])")
    call disp%show( isUniqueAll([1]) )

    call disp%show("isUniqueAll([1, 2, 3, 4, 0])")
    call disp%show( isUniqueAll([1, 2, 3, 4, 0]) )

    call disp%show("isUniqueAll([1, 2, 2, 4, 0])")
    call disp%show( isUniqueAll([1, 2, 2, 4, 0]) )

    call disp%show("isUniqueAll([1, 1, 1, 1, 0])")
    call disp%show( isUniqueAll([1, 1, 1, 1, 0]) )

    call disp%show("isUniqueAll([1, 1, 1, 1])")
    call disp%show( isUniqueAll([1, 1, 1, 1]) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in logical array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("isUniqueAll([logical::])")
    call disp%show( isUniqueAll([logical::]) )

    call disp%show("isUniqueAll([.true.])")
    call disp%show( isUniqueAll([.true.]) )

    call disp%show("isUniqueAll([.true., .false., .true., .false.])")
    call disp%show( isUniqueAll([.true., .false., .true., .false.]) )

    call disp%show("isUniqueAll([.false., .false., .false., .true.])")
    call disp%show( isUniqueAll([.false., .false., .false., .true.]) )

    call disp%show("isUniqueAll([.true., .true., .true., .false.])")
    call disp%show( isUniqueAll([.true., .true., .true., .false.]) )

    call disp%show("isUniqueAll([.false., .false., .false.])")
    call disp%show( isUniqueAll([.false., .false., .false.]) )

    call disp%show("isUniqueAll([.true., .true., .true.])")
    call disp%show( isUniqueAll([.true., .true., .true.]) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in complex array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("isUniqueAll([complex::])")
    call disp%show( isUniqueAll([complex::]) )

    call disp%show("isUniqueAll([(1., -1)])")
    call disp%show( isUniqueAll([(1., -1)]) )

    call disp%show("isUniqueAll([(1., -1.), (2., -2.), (3., -3.), (4., -4.), (0., 0.)])")
    call disp%show( isUniqueAll([(1., -1.), (2., -2.), (3., -3.), (4., -4.), (0., 0.)]) )

    call disp%show("isUniqueAll([(1., -1.), (2., -2.), (2., -2.), (4., -4.), (0., 0.)])")
    call disp%show( isUniqueAll([(1., -1.), (2., -2.), (2., -2.), (4., -4.), (0., 0.)]) )

    call disp%show("isUniqueAll([(1., -1.), (1., -1.), (1., -1.), (1., -1.), (0., 0.)])")
    call disp%show( isUniqueAll([(1., -1.), (1., -1.), (1., -1.), (1., -1.), (0., 0.)]) )

    call disp%show("isUniqueAll([(1., -1.), (1., -1.), (1., -1.), (1., -1.)])")
    call disp%show( isUniqueAll([(1., -1.), (1., -1.), (1., -1.), (1., -1.)]) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("isUniqueAll([real::])")
    call disp%show( isUniqueAll([real::]) )

    call disp%show("isUniqueAll([1.])")
    call disp%show( isUniqueAll([1.]) )

    call disp%show("isUniqueAll([1., 2., 3., 4., 0.])")
    call disp%show( isUniqueAll([1., 2., 3., 4., 0.]) )

    call disp%show("isUniqueAll([1., 2., 2., 4., 0.])")
    call disp%show( isUniqueAll([1., 2., 2., 4., 0.]) )

    call disp%show("isUniqueAll([1., 1., 1., 1., 0.])")
    call disp%show( isUniqueAll([1., 1., 1., 1., 0.]) )

    call disp%show("isUniqueAll([1., 1., 1., 1.])")
    call disp%show( isUniqueAll([1., 1., 1., 1.]) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("isUniqueAll([real(RKG) :: 1.01, 1.04, 0.98, 1.0, 1.02, 2.], iseq = iseq_RK)")
    call disp%show( isUniqueAll([real(RKG) :: 1.01, 1.04, 0.98, 1.0, 1.02, 2.], iseq = iseq_RK) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique case-insensitive instances within the character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("isUniqueAll('ABBAbbA', iseq = iseq_SK)")
    call disp%show( isUniqueAll('ABBAbbA', iseq = iseq_SK), deliml = SK_"""" )

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