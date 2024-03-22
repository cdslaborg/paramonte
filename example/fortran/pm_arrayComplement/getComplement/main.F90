program example

    use pm_kind, only: LK
    use pm_kind, only: SK, IK, LK, CK => CKS, RK => RKS ! all processor types and kinds are supported.
    use pm_arrayComplement, only: getComplement
    use pm_arrayUnique, only: getUnique
    use pm_io, only: display_type
    use pm_arraySort, only: setSorted

    implicit none

    character(:, SK), allocatable   :: StrA_SK      , StrB_SK
    character(2, SK), allocatable   :: SetA_SK(:)   , SetB_SK(:)
    integer(IK)     , allocatable   :: SetA_IK(:)   , SetB_IK(:)
    complex(CK)     , allocatable   :: SetA_CK(:)   , SetB_CK(:)
    real(RK)        , allocatable   :: SetA_RK(:)   , SetB_RK(:)

    type(display_type)              :: disp

    disp = display_type(file = "main.out.F90")

    SetA_SK = ["AA", "BB", "CC", "AA", "FF", "CC", "DD"]
    SetB_SK = ["aa", "CC", "CC", "CC", "DD", "EE", "EE"]

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find the complement of string scalar A in B: B \ A.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    StrA_SK = "AADBCCFEEZ"
    StrB_SK = "EADGBBCGCE"

    call disp%skip()
    call disp%show("StrA_SK")
    call disp%show( StrA_SK, deliml = SK_"""" )
    call disp%show("StrB_SK")
    call disp%show( StrB_SK, deliml = SK_"""" )
    call disp%show("getComplement(StrA_SK, StrB_SK)")
    call disp%show( getComplement(StrA_SK, StrB_SK), deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("call setSorted(StrA_SK)")
                    call setSorted(StrA_SK)
    call disp%show("StrA_SK")
    call disp%show( StrA_SK, deliml = SK_"""" )
    call disp%show("call setSorted(StrB_SK)")
                    call setSorted(StrB_SK)
    call disp%show("StrB_SK")
    call disp%show( StrB_SK, deliml = SK_"""" )
    call disp%show("getComplement(StrA_SK, StrB_SK, sorted = .true._LK, unique = .false._LK)")
    call disp%show( getComplement(StrA_SK, StrB_SK, sorted = .true._LK, unique = .false._LK), deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("StrA_SK = getUnique(StrA_SK)")
                    StrA_SK = getUnique(StrA_SK)
    call disp%show("StrA_SK")
    call disp%show( StrA_SK, deliml = SK_"""" )
    call disp%show("StrB_SK = getUnique(StrB_SK)")
                    StrB_SK = getUnique(StrB_SK)
    call disp%show("StrB_SK")
    call disp%show( StrB_SK, deliml = SK_"""" )
    call disp%show("getComplement(StrA_SK, StrB_SK, sorted = .true._LK, unique = .true._LK)")
    call disp%show( getComplement(StrA_SK, StrB_SK, sorted = .true._LK, unique = .true._LK), deliml = SK_"""" )
    call disp%skip()

    StrA_SK = "AADBCCFEEZ"
    StrB_SK = "AAAAAAAAAA"

    call disp%skip()
    call disp%show("StrA_SK")
    call disp%show( StrA_SK, deliml = SK_"""" )
    call disp%show("StrB_SK")
    call disp%show( StrB_SK, deliml = SK_"""" )
    call disp%show("getComplement(StrA_SK, StrB_SK, sorted = .true._LK, unique = .true._LK) ! The result is wrong because both sets must contain unique elements when `unique = .true.`.")
    call disp%show( getComplement(StrA_SK, StrB_SK, sorted = .true._LK, unique = .true._LK), deliml = SK_"""" )
    call disp%skip()

    StrA_SK = "AADBCCFEEZ"
    StrB_SK = "BAAAAAAAAA"

    call disp%skip()
    call disp%show("StrA_SK")
    call disp%show( StrA_SK, deliml = SK_"""" )
    call disp%show("StrB_SK")
    call disp%show( StrB_SK, deliml = SK_"""" )
    call disp%show("getComplement(StrA_SK, StrB_SK, sorted = .true._LK, unique = .false._LK) ! The result is wrong because both sets must be similarly sorted when `sorted = .true.`.")
    call disp%show( getComplement(StrA_SK, StrB_SK, sorted = .true._LK, unique = .false._LK), deliml = SK_"""" )
    call disp%skip()

    StrA_SK = "AABCCDEEFZ"
    StrB_SK = "ZFEEDCCBAA"

    call disp%skip()
    call disp%show("StrA_SK")
    call disp%show( StrA_SK, deliml = SK_"""" )
    call disp%show("StrB_SK")
    call disp%show( StrB_SK, deliml = SK_"""" )
    call disp%show("getComplement(StrA_SK, StrB_SK, sorted = .true._LK, unique = .false._LK) ! The result is wrong because both sets must be similarly sorted when `sorted = .true.`.")
    call disp%show( getComplement(StrA_SK, StrB_SK, sorted = .true._LK, unique = .false._LK), deliml = SK_"""" )
    call disp%show("getComplement(StrA_SK, StrB_SK, sorted = .false._LK, unique = .false._LK)")
    call disp%show( getComplement(StrA_SK, StrB_SK, sorted = .false._LK, unique = .false._LK), deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find the complement of string vector A in B: B \ A.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("SetA_SK")
    call disp%show( SetA_SK, deliml = SK_"""" )
    call disp%show("SetB_SK")
    call disp%show( SetB_SK, deliml = SK_"""" )
    call disp%show("getComplement(SetA_SK, SetB_SK)")
    call disp%show( getComplement(SetA_SK, SetB_SK), deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("call setSorted(SetA_SK)")
                    call setSorted(SetA_SK)
    call disp%show("SetA_SK")
    call disp%show( SetA_SK, deliml = SK_"""" )
    call disp%show("call setSorted(SetB_SK)")
                    call setSorted(SetB_SK)
    call disp%show("SetB_SK")
    call disp%show( SetB_SK, deliml = SK_"""" )
    call disp%show("getComplement(SetA_SK, SetB_SK, sorted = .true._LK, unique = .false._LK)")
    call disp%show( getComplement(SetA_SK, SetB_SK, sorted = .true._LK, unique = .false._LK), deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("SetA_SK = getUnique(SetA_SK)")
                    SetA_SK = getUnique(SetA_SK)
    call disp%show("SetA_SK")
    call disp%show( SetA_SK, deliml = SK_"""" )
    call disp%show("SetB_SK = getUnique(SetB_SK)")
                    SetB_SK = getUnique(SetB_SK)
    call disp%show("SetB_SK")
    call disp%show( SetB_SK, deliml = SK_"""" )
    call disp%show("getComplement(SetA_SK, SetB_SK, sorted = .true._LK, unique = .true._LK)")
    call disp%show( getComplement(SetA_SK, SetB_SK, sorted = .true._LK, unique = .true._LK), deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find the complement of integer vector A in B: B \ A.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    SetA_IK = [integer(IK) :: 4, 5, 4, 1, 2, 1, 5]
    SetB_IK = [integer(IK) :: 6, 3, 4, 2, 1, 2, 6]

    call disp%skip()
    call disp%show("SetA_IK")
    call disp%show( SetA_IK )
    call disp%show("SetB_IK")
    call disp%show( SetB_IK )
    call disp%show("getComplement(SetA_IK, SetB_IK)")
    call disp%show( getComplement(SetA_IK, SetB_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setSorted(SetA_IK)")
                    call setSorted(SetA_IK)
    call disp%show("SetA_IK")
    call disp%show( SetA_IK )
    call disp%show("call setSorted(SetB_IK)")
                    call setSorted(SetB_IK)
    call disp%show("SetB_IK")
    call disp%show( SetB_IK )
    call disp%show("getComplement(SetA_IK, SetB_IK, sorted = .true._LK, unique = .false._LK)")
    call disp%show( getComplement(SetA_IK, SetB_IK, sorted = .true._LK, unique = .false._LK) )
    call disp%skip()

    call disp%skip()
    call disp%show("SetA_IK = getUnique(SetA_IK)")
                    SetA_IK = getUnique(SetA_IK)
    call disp%show("SetA_IK")
    call disp%show( SetA_IK )
    call disp%show("SetB_IK = getUnique(SetB_IK)")
                    SetB_IK = getUnique(SetB_IK)
    call disp%show("SetB_IK")
    call disp%show( SetB_IK )
    call disp%show("getComplement(SetA_IK, SetB_IK, sorted = .true._LK, unique = .true._LK)")
    call disp%show( getComplement(SetA_IK, SetB_IK, sorted = .true._LK, unique = .true._LK) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find the complement of complex vector A in B: B \ A.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    SetA_CK = [complex(CK) :: (4, -4), (5, -5), (4, -4), (1, -1), (2, -2), (8, -8), (1, -1), (5, -5)]
    SetB_CK = [complex(CK) :: (6, -6), (3, -3), (4, -4), (2, -2), (1, -1), (7, -7), (2, -2), (6, -6)]

    call disp%skip()
    call disp%show("SetA_CK")
    call disp%show( SetA_CK )
    call disp%show("SetB_CK")
    call disp%show( SetB_CK )
    call disp%show("getComplement(SetA_CK, SetB_CK)")
    call disp%show( getComplement(SetA_CK, SetB_CK) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setSorted(SetA_CK)")
                    call setSorted(SetA_CK)
    call disp%show("SetA_CK")
    call disp%show( SetA_CK )
    call disp%show("call setSorted(SetB_CK)")
                    call setSorted(SetB_CK)
    call disp%show("SetB_CK")
    call disp%show( SetB_CK )
    call disp%show("getComplement(SetA_CK, SetB_CK, sorted = .true._LK, unique = .false._LK)")
    call disp%show( getComplement(SetA_CK, SetB_CK, sorted = .true._LK, unique = .false._LK) )
    call disp%skip()

    call disp%skip()
    call disp%show("SetA_CK = getUnique(SetA_CK)")
                    SetA_CK = getUnique(SetA_CK)
    call disp%show("SetA_CK")
    call disp%show( SetA_CK )
    call disp%show("SetB_CK = getUnique(SetB_CK)")
                    SetB_CK = getUnique(SetB_CK)
    call disp%show("SetB_CK")
    call disp%show( SetB_CK )
    call disp%show("getComplement(SetA_CK, SetB_CK, sorted = .true._LK, unique = .true._LK)")
    call disp%show( getComplement(SetA_CK, SetB_CK, sorted = .true._LK, unique = .true._LK) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find the complement of real    vector A in B: B \ A.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    SetA_RK = [real(RK) :: 4, 5, 4, 1, 2, 8, 1, 5]
    SetB_RK = [real(RK) :: 6, 3, 4, 2, 1, 7, 2, 6]

    call disp%skip()
    call disp%show("SetA_RK")
    call disp%show( SetA_RK )
    call disp%show("SetB_RK")
    call disp%show( SetB_RK )
    call disp%show("getComplement(SetA_RK, SetB_RK)")
    call disp%show( getComplement(SetA_RK, SetB_RK) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setSorted(SetA_RK)")
                    call setSorted(SetA_RK)
    call disp%show("SetA_RK")
    call disp%show( SetA_RK )
    call disp%show("call setSorted(SetB_RK)")
                    call setSorted(SetB_RK)
    call disp%show("SetB_RK")
    call disp%show( SetB_RK )
    call disp%show("getComplement(SetA_RK, SetB_RK, sorted = .true._LK, unique = .false._LK)")
    call disp%show( getComplement(SetA_RK, SetB_RK, sorted = .true._LK, unique = .false._LK) )
    call disp%skip()

    call disp%skip()
    call disp%show("SetA_RK = getUnique(SetA_RK)")
                    SetA_RK = getUnique(SetA_RK)
    call disp%show("SetA_RK")
    call disp%show( SetA_RK )
    call disp%show("SetB_RK = getUnique(SetB_RK)")
                    SetB_RK = getUnique(SetB_RK)
    call disp%show("SetB_RK")
    call disp%show( SetB_RK )
    call disp%show("getComplement(SetA_RK, SetB_RK, sorted = .true._LK, unique = .true._LK)")
    call disp%show( getComplement(SetA_RK, SetB_RK, sorted = .true._LK, unique = .true._LK) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find the complement using a custom equivalence check.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    StrA_SK = "AADBCCFEEZ"
    StrB_SK = "eadgbbcgce"

    call disp%skip()
    call disp%show("StrA_SK")
    call disp%show( StrA_SK, deliml = SK_"""" )
    call disp%show("StrB_SK")
    call disp%show( StrB_SK, deliml = SK_"""" )
    call disp%show("getComplement(StrA_SK, StrB_SK)")
    call disp%show( getComplement(StrA_SK, StrB_SK), deliml = SK_"""" )
    call disp%show("getComplement(StrA_SK, StrB_SK, iseq = iseq_SK)")
    call disp%show( getComplement(StrA_SK, StrB_SK, iseq = iseq_SK), deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("call setSorted(StrA_SK)")
                    call setSorted(StrA_SK)
    call disp%show("StrA_SK")
    call disp%show( StrA_SK, deliml = SK_"""" )
    call disp%show("call setSorted(StrB_SK)")
                    call setSorted(StrB_SK)
    call disp%show("StrB_SK")
    call disp%show( StrB_SK, deliml = SK_"""" )
    call disp%show("getComplement(StrA_SK, StrB_SK, sorted = .true._LK, unique = .false._LK)")
    call disp%show( getComplement(StrA_SK, StrB_SK, sorted = .true._LK, unique = .false._LK), deliml = SK_"""" )
    call disp%show("getComplement(StrA_SK, StrB_SK, sorted = .true._LK, unique = .false._LK, iseq = iseq_SK)")
    call disp%show( getComplement(StrA_SK, StrB_SK, sorted = .true._LK, unique = .false._LK, iseq = iseq_SK), deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("StrA_SK = getUnique(StrA_SK)")
                    StrA_SK = getUnique(StrA_SK)
    call disp%show("StrA_SK")
    call disp%show( StrA_SK, deliml = SK_"""" )
    call disp%show("StrB_SK = getUnique(StrB_SK)")
                    StrB_SK = getUnique(StrB_SK)
    call disp%show("StrB_SK")
    call disp%show( StrB_SK, deliml = SK_"""" )
    call disp%show("getComplement(StrA_SK, StrB_SK, sorted = .true._LK, unique = .true._LK)")
    call disp%show( getComplement(StrA_SK, StrB_SK, sorted = .true._LK, unique = .true._LK), deliml = SK_"""" )
    call disp%show("getComplement(StrA_SK, StrB_SK, sorted = .true._LK, unique = .true._LK, iseq = iseq_SK)")
    call disp%show( getComplement(StrA_SK, StrB_SK, sorted = .true._LK, unique = .true._LK, iseq = iseq_SK), deliml = SK_"""" )
    call disp%skip()

contains

    pure function iseq_SK(element1, element2) result(iseq)
        use pm_strASCII, only: getCharLower
        character(1, SK)    , intent(in)    :: element1, element2
        logical(LK)                         :: iseq
        iseq = getCharLower(element1) == getCharLower(element2)
    end function

end program example