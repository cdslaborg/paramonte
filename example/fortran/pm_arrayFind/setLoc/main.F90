program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! All kinds are supported.
    use pm_kind, only: IK ! All kinds are supported.
    use pm_kind, only: CK ! All kinds are supported.
    use pm_kind, only: RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayResize, only: setResized
    use pm_arrayFind, only: setLoc

    implicit none

    integer(IK)     , allocatable   :: instance(:)  , loc(:)                ! Must be of default kind IK

    character(:, SK), allocatable   :: string_SK    , stringPattern_SK
    character(9, SK), allocatable   :: array_SK(:)  , patvec_SK(:)    ! Can be any processor-supported kind.
    integer(IK)     , allocatable   :: array_IK(:)  , patvec_IK(:)    ! Can be any processor-supported kind.
    complex(CK)     , allocatable   :: array_CK(:)  , patvec_CK(:)    ! Can be any processor-supported kind.
    real(RK)        , allocatable   :: array_RK(:)  , patvec_RK(:)    ! Can be any processor-supported kind.
    logical(LK)     , allocatable   :: array_LK(:)  , patvec_LK(:)

    integer(IK) :: nloc
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    string_SK = "ParaMonte is a Machine Learning Library "
    array_SK = ["ParaMonte", "XXXXXXXXX", "is       ", "XXXXXXXXX", "a        ", "XXXXXXXXX", "Monte    ", "XXXXXXXXX", "Carlo    ", "XXXXXXXXX", "Library. ", "XXXXXXXXX"]
    array_IK = [1_IK, 0_IK, 2_IK, 0_IK, 3_IK, 0_IK, 4_IK]
    array_RK = [1._RK, 0._RK, 2._RK, 0._RK, 3._RK, 0._RK, 4._RK]
    array_CK = [(1._CK, -1._CK), (0._CK, -0._CK), (2._CK, -2._CK), (0._CK, -0._CK), (3._CK, -3._CK), (0._CK, -0._CK), (4._CK, -4._CK)]
    array_LK = [.false._LK, .true._LK, .false._LK, .true._LK, .false._LK, .true._LK, .false._LK]

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find all instances of pattern in array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    stringPattern_SK = " "
    patvec_SK = ["XXXXXXXXX"]
    patvec_IK = [0_IK]
    patvec_RK = [0._RK]
    patvec_CK = [(0._CK, -0._CK)]
    patvec_LK = [.true._LK]

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringPattern_SK")
    call disp%show( stringPattern_SK, deliml = SK_"""" )
    call disp%show("call setResized(loc, 0_IK)")
                    call setResized(loc, 0_IK)
    call disp%show("call setLoc(loc, nloc, string_SK, stringPattern_SK, blindness = 1_IK)")
                    call setLoc(loc, nloc, string_SK, stringPattern_SK, blindness = 1_IK)
    call disp%show("loc(1 : nloc)")
    call disp%show( loc(1 : nloc) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_SK")
    call disp%show( array_SK, deliml = SK_"""" )
    call disp%show("patvec_SK")
    call disp%show( patvec_SK, deliml = SK_"""" )
    call disp%show("call setLoc(loc, nloc, array_SK, patvec_SK, blindness = 1_IK)")
                    call setLoc(loc, nloc, array_SK, patvec_SK, blindness = 1_IK)
    call disp%show("loc(1 : nloc)")
    call disp%show( loc(1 : nloc) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find logical array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_LK")
    call disp%show( array_LK )
    call disp%show("patvec_LK")
    call disp%show( patvec_LK )
    call disp%show("call setLoc(loc, nloc, array_LK, patvec_LK, blindness = 1_IK)")
                    call setLoc(loc, nloc, array_LK, patvec_LK, blindness = 1_IK)
    call disp%show("loc(1 : nloc)")
    call disp%show( loc(1 : nloc) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_IK")
    call disp%show( array_IK )
    call disp%show("patvec_IK")
    call disp%show( patvec_IK )
    call disp%show("call setLoc(loc, nloc, array_IK, patvec_IK, blindness = 1_IK)")
                    call setLoc(loc, nloc, array_IK, patvec_IK, blindness = 1_IK)
    call disp%show("loc(1 : nloc)")
    call disp%show( loc(1 : nloc) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find complex array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_CK")
    call disp%show( array_CK )
    call disp%show("patvec_CK")
    call disp%show( patvec_CK )
    call disp%show("call setLoc(loc, nloc, array_CK, patvec_CK, blindness = 1_IK)")
                    call setLoc(loc, nloc, array_CK, patvec_CK, blindness = 1_IK)
    call disp%show("loc(1 : nloc)")
    call disp%show( loc(1 : nloc) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%show("! Find real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_RK")
    call disp%show( array_RK )
    call disp%show("patvec_RK")
    call disp%show( patvec_RK )
    call disp%show("call setLoc(loc, nloc, array_RK, patvec_RK, blindness = 1_IK)")
                    call setLoc(loc, nloc, array_RK, patvec_RK, blindness = 1_IK)
    call disp%show("loc(1 : nloc)")
    call disp%show( loc(1 : nloc) )


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find only particular instances of pattern in array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    string_SK = "A_A_A_A_A_A_A_A_A"
    array_SK = ["A", "_", "A", "_", "A", "_", "A", "_", "A", "_", "A", "_", "A", "_", "A", "_", "A"]
    array_IK = [0_IK, 1_IK, 0_IK, 2_IK, 3_IK, 0_IK, 4_IK, 5_IK, 0_IK, 0_IK]
    array_RK = [0._RK, 1._RK, 0._RK, 2._RK, 3._RK, 0._RK, 4._RK, 5._RK, 0._RK, 0._RK]
    array_CK = [(0._CK, -0._CK), (1._CK, -1._CK), (0._CK, -0._CK), (2._CK, -2._CK), (3._CK, -3._CK), (0._CK, -0._CK), (4._CK, -4._CK), (5._CK, -5._CK), (0._CK, -0._CK), (0._CK, -0._CK)]
    array_LK = [.false._LK, .true._LK, .false._LK, .true._LK, .true._LK, .false._LK, .true._LK, .true._LK, .false._LK, .false._LK]

    stringPattern_SK = "_"
    patvec_SK = ["_"]
    patvec_IK = [0_IK]
    patvec_RK = [0._RK]
    patvec_CK = [(0._CK, -0._CK)]
    patvec_LK = [.false._LK]

    instance = [-3_IK, 2_IK, -4_IK, 20_IK] ! Find at the second occurrence from the beginning and the third occurrence from the end. Identical instances yield identical indices. Out-of-bound instance indices are ignored.

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringPattern_SK")
    call disp%show( stringPattern_SK, deliml = SK_"""" )
    call disp%show("instance ! Identical instances yield identical indices. Out-of-bound instance indices are ignored.")
    call disp%show( instance )
    call disp%show("call setLoc(loc, nloc, string_SK, stringPattern_SK, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)")
                    call setLoc(loc, nloc, string_SK, stringPattern_SK, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
    call disp%show("loc(1 : nloc)")
    call disp%show( loc(1 : nloc) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find vector `pattern` from character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_SK")
    call disp%show( array_SK, deliml = SK_"""" )
    call disp%show("patvec_SK")
    call disp%show( patvec_SK, deliml = SK_"""" )
    call disp%show("instance ! Identical instances yield identical indices. Out-of-bound instance indices are ignored.")
    call disp%show( instance )
    call disp%show("call setLoc(loc, nloc, array_SK, patvec_SK, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)")
                    call setLoc(loc, nloc, array_SK, patvec_SK, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
    call disp%show("loc(1 : nloc)")
    call disp%show( loc(1 : nloc) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find character array with scalar `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_SK")
    call disp%show( array_SK, deliml = SK_"""" )
    call disp%show("patvec_SK(1)")
    call disp%show( patvec_SK(1), deliml = SK_"""" )
    call disp%show("instance ! Identical instances yield identical indices. Out-of-bound instance indices are ignored.")
    call disp%show( instance )
    call disp%show("call setLoc(loc, nloc, array_SK, patvec_SK(1), instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)")
                    call setLoc(loc, nloc, array_SK, patvec_SK(1), instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
    call disp%show("loc(1 : nloc)")
    call disp%show( loc(1 : nloc) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find logical array with vector `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_LK")
    call disp%show( array_LK )
    call disp%show("patvec_LK")
    call disp%show( patvec_LK )
    call disp%show("instance ! Identical instances yield identical indices. Out-of-bound instance indices are ignored.")
    call disp%show( instance )
    call disp%show("call setLoc(loc, nloc, array_LK, patvec_LK, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)")
                    call setLoc(loc, nloc, array_LK, patvec_LK, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
    call disp%show("loc(1 : nloc)")
    call disp%show( loc(1 : nloc) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find logical array with scalar `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_LK")
    call disp%show( array_LK )
    call disp%show("patvec_LK(1)")
    call disp%show( patvec_LK(1) )
    call disp%show("instance ! Identical instances yield identical indices. Out-of-bound instance indices are ignored.")
    call disp%show( instance )
    call disp%show("call setLoc(loc, nloc, array_LK, patvec_LK(1), instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)")
                    call setLoc(loc, nloc, array_LK, patvec_LK(1), instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
    call disp%show("loc(1 : nloc)")
    call disp%show( loc(1 : nloc) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find integer array with vector `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_IK")
    call disp%show( array_IK )
    call disp%show("patvec_IK")
    call disp%show( patvec_IK )
    call disp%show("instance ! Identical instances yield identical indices. Out-of-bound instance indices are ignored.")
    call disp%show( instance )
    call disp%show("call setLoc(loc, nloc, array_IK, patvec_IK, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)")
                    call setLoc(loc, nloc, array_IK, patvec_IK, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
    call disp%show("loc(1 : nloc)")
    call disp%show( loc(1 : nloc) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find integer array with scalar `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_IK")
    call disp%show( array_IK )
    call disp%show("patvec_IK(1)")
    call disp%show( patvec_IK(1) )
    call disp%show("instance ! Identical instances yield identical indices. Out-of-bound instance indices are ignored.")
    call disp%show( instance )
    call disp%show("call setLoc(loc, nloc, array_IK, patvec_IK(1), instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)")
                    call setLoc(loc, nloc, array_IK, patvec_IK(1), instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
    call disp%show("loc(1 : nloc)")
    call disp%show( loc(1 : nloc) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find complex array with vector `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_CK")
    call disp%show( array_CK )
    call disp%show("patvec_CK")
    call disp%show( patvec_CK )
    call disp%show("instance ! Identical instances yield identical indices. Out-of-bound instance indices are ignored.")
    call disp%show( instance )
    call disp%show("call setLoc(loc, nloc, array_CK, patvec_CK, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)")
                    call setLoc(loc, nloc, array_CK, patvec_CK, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
    call disp%show("loc(1 : nloc)")
    call disp%show( loc(1 : nloc) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find complex array with scalar `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_CK")
    call disp%show( array_CK )
    call disp%show("patvec_CK(1)")
    call disp%show( patvec_CK(1) )
    call disp%show("instance ! Identical instances yield identical indices. Out-of-bound instance indices are ignored.")
    call disp%show( instance )
    call disp%show("call setLoc(loc, nloc, array_CK, patvec_CK(1), instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)")
                    call setLoc(loc, nloc, array_CK, patvec_CK(1), instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
    call disp%show("loc(1 : nloc)")
    call disp%show( loc(1 : nloc) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find real array with vector `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_RK")
    call disp%show( array_RK )
    call disp%show("patvec_RK")
    call disp%show( patvec_RK )
    call disp%show("instance ! Identical instances yield identical indices. Out-of-bound instance indices are ignored.")
    call disp%show( instance )
    call disp%show("call setLoc(loc, nloc, array_RK, patvec_RK, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)")
                    call setLoc(loc, nloc, array_RK, patvec_RK, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
    call disp%show("loc(1 : nloc)")
    call disp%show( loc(1 : nloc) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find real array with scalar `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_RK")
    call disp%show( array_RK )
    call disp%show("patvec_RK(1)")
    call disp%show( patvec_RK(1) )
    call disp%show("instance ! Identical instances yield identical indices. Out-of-bound instance indices are ignored.")
    call disp%show( instance )
    call disp%show("call setLoc(loc, nloc, array_RK, patvec_RK(1), instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)")
                    call setLoc(loc, nloc, array_RK, patvec_RK(1), instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
    call disp%show("loc(1 : nloc)")
    call disp%show( loc(1 : nloc) )


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Adjust the blindness to find exclusively non-overlapping instances of `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()


    string_SK = "AAAAAAAA"
    stringPattern_SK = "AAA"

    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringPattern_SK")
    call disp%show( stringPattern_SK, deliml = SK_"""" )
    call disp%show("call setLoc(loc, nloc, string_SK, stringPattern_SK, blindness = 1_IK) ! The default blindness episode after each detection is `blindness = 1_IK`")
                    call setLoc(loc, nloc, string_SK, stringPattern_SK, blindness = 1_IK)
    call disp%show("loc(1 : nloc)")
    call disp%show( loc(1 : nloc) )

    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringPattern_SK")
    call disp%show( stringPattern_SK, deliml = SK_"""" )
    call disp%show("call setLoc(loc, nloc, string_SK, stringPattern_SK, blindness = size(stringPattern_SK, kind = IK)) ! Find only non-overlapping patterns.")
                    call setLoc(loc, nloc, string_SK, stringPattern_SK, blindness = len(stringPattern_SK, kind = IK))
    call disp%show("loc(1 : nloc)")
    call disp%show( loc(1 : nloc) )

    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringPattern_SK")
    call disp%show( stringPattern_SK, deliml = SK_"""" )
    call disp%show("call setLoc(loc, nloc, string_SK, stringPattern_SK, blindness = 2_IK) ! Find instances with jumps of size 2.")
                    call setLoc(loc, nloc, string_SK, stringPattern_SK, blindness = 2_IK)
    call disp%show("loc(1 : nloc)")
    call disp%show( loc(1 : nloc) )


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find specific instances with a user-defined equivalence test.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find case-insensitive instances of vector `pattern` within the character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    string_SK = "ABBAbbA"
    stringPattern_SK = "bb"

    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringPattern_SK")
    call disp%show( stringPattern_SK, deliml = SK_"""" )
    call disp%show("call setLoc(loc, nloc, string_SK, stringPattern_SK, iseq = iseq_SK, blindness = 1_IK)")
                    call setLoc(loc, nloc, string_SK, stringPattern_SK, iseq = iseq_SK, blindness = 1_IK)
    call disp%show("loc(1 : nloc)")
    call disp%show( loc(1 : nloc) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find specific instances of vector `pattern` within the real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    array_RK = [0._RK, 1.01_RK, 1.04_RK, 0.98_RK, 1.0_RK, 1.02_RK]
    patvec_RK = [-1._RK, 1._RK]
    instance = [-2_IK]

    call disp%show("array_RK")
    call disp%show( array_RK )
    call disp%show("patvec_RK")
    call disp%show( patvec_RK )
    call disp%show("instance ! Identical instances yield identical indices. Out-of-bound instance indices are ignored.")
    call disp%show( instance )
    call disp%show("call setLoc(loc, nloc, array_RK, patvec_RK, iseq = iseq_vec_RK, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)")
                    call setLoc(loc, nloc, array_RK, patvec_RK, iseq = iseq_vec_RK, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
    call disp%show("loc(1 : nloc)")
    call disp%show( loc(1 : nloc) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find specific instances of scalar `pattern` within the real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_RK")
    call disp%show( array_RK )
    call disp%show("patvec_RK(1)")
    call disp%show( patvec_RK(1) )
    call disp%show("instance ! Identical instances yield identical indices. Out-of-bound instance indices are ignored.")
    call disp%show( instance )
    call disp%show("call setLoc(loc, nloc, array_RK, patvec_RK(1), iseq = iseq_RK, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)")
                    call setLoc(loc, nloc, array_RK, patvec_RK(1), iseq = iseq_RK, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
    call disp%show("loc(1 : nloc)")
    call disp%show( loc(1 : nloc) )

contains

    pure function iseq_SK(ArraySegment, pattern) result(equivalent)
        use pm_strASCII, only: getStrLower
        character(*, SK), intent(in)    :: pattern, ArraySegment
        logical(LK)                     :: equivalent
        equivalent = pattern == getStrLower(ArraySegment)
    end function

    function iseq_RK(arraysegment, pattern) result(equivalent)
        real(RK)        , intent(in)    :: pattern, arraySegment
        logical(LK)                     :: equivalent
        equivalent = abs(abs(pattern) - abs(arraySegment)) < 0.05_RK
    end function

    function iseq_vec_RK(ArraySegment, pattern, lenPattern) result(equivalent)
        integer(IK)     , intent(in)    :: lenPattern
        real(RK)        , intent(in)    :: pattern(lenPattern), ArraySegment(lenPattern)
        logical(LK)                     :: equivalent
        equivalent = all(abs(abs(pattern) - abs(ArraySegment)) < 0.05_RK)
    end function

end program example