program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! all other character kinds are also supported.
    use pm_kind, only: IK ! All kinds are supported.
    use pm_kind, only: CK ! All kinds are supported.
    use pm_kind, only: RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arraySplit, only: setSplit

    ! Bypass gfortran 11 pdt bug

#if PDT_ENABLED
    use pm_container, only: cvi_pdt ! Containers of all kinds are supported.
    use pm_container, only: cvc_pdt ! Containers of all kinds are supported.
    use pm_container, only: cvr_pdt ! Containers of all kinds are supported.
    use pm_container, only: cvs_pdt ! Containers of all kinds are supported.
    use pm_container, only: cvl_pdt ! Containers of all kinds are supported.
    use pm_container, only: css_pdt ! Containers of all kinds are supported.
#endif

    implicit none

    integer(IK)     , allocatable   :: instance(:)      ! Must be of default kind IK
    integer(IK)     , allocatable   :: splitIndex(:,:)  ! Array of shape (2,*) containing the start and end indices of the input array corresponding to the split segments.

    character(:, SK), allocatable   :: string_SK    , stringSep_SK    ! Can be any processor-supported kind.
    character(9, SK), allocatable   :: array_SK(:)  , sepArray_SK(:)  ! Can be any processor-supported kind.
    integer(IK)     , allocatable   :: array_IK(:)  , sepArray_IK(:)  ! Can be any processor-supported kind.
    complex(CK)     , allocatable   :: array_CK(:)  , sepArray_CK(:)  ! Can be any processor-supported kind.
    real(RK)        , allocatable   :: array_RK(:)  , sepArray_RK(:)  ! Can be any processor-supported kind.
    logical(LK)     , allocatable   :: array_LK(:)  , sepArray_LK(:)  ! Can be any processor-supported kind.

#if PDT_ENABLED
    type(css_pdt), allocatable :: parts_SSK(:)
    type(cvs_pdt), allocatable :: parts_VSK(:)
    type(cvi_pdt), allocatable :: parts_VIK(:)
    type(cvc_pdt), allocatable :: parts_VCK(:)
    type(cvr_pdt), allocatable :: parts_VRK(:)
    type(cvl_pdt), allocatable :: parts_VLK(:)
#endif

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Split all instances of Sep in array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    string_SK = "ParaMonte is a Machine Learning Library "
    array_SK = ["ParaMonte", "XXXXXXXXX", "is       ", "XXXXXXXXX", "a        ", "XXXXXXXXX", "Monte    ", "XXXXXXXXX", "Carlo    ", "XXXXXXXXX", "Library. ", "XXXXXXXXX"]
    array_IK = [integer(IK) :: 1, 0, 2, 1, 0, 4, 5, 1, 0, 7, 8, 9, 10]
    array_LK = [logical(LK) :: .false., .true., .true., .false., .true., .true., .true., .false., .true.]
    array_CK = [complex(CK) :: (1., -1.), (0., -0.), (2., -2.), (1., -1.), (0., -0.), (4., -4.), (5., -5.), (6., -6.), (0., -0.), (4., -4.)]
    array_RK = [real(RK)    :: 1, 0, 2, 1, 0, 4, 5, 1, 0, 7, 8, 9, 10]

    stringSep_SK = " "
    sepArray_SK = ["XXXXXXXXX"]
    sepArray_IK = [integer(IK) :: 1, 0]
    sepArray_LK = [logical(LK) :: .false., .true.]
    sepArray_CK = [complex(CK) :: (1., -1.), (0., -0.)]
    sepArray_RK = [real(RK)    :: 1, 0]

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Split character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringSep_SK")
    call disp%show( stringSep_SK, deliml = SK_"""" )
    call disp%show("call setSplit(parts_SSK, string_SK, stringSep_SK)")
                    call setSplit(parts_SSK, string_SK, stringSep_SK)
    call disp%show("parts_SSK")
    call disp%show( parts_SSK, deliml = SK_"""" )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringSep_SK")
    call disp%show( stringSep_SK, deliml = SK_"""" )
    call disp%show("call setSplit(splitIndex, string_SK, stringSep_SK)")
                    call setSplit(splitIndex, string_SK, stringSep_SK)
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringSep_SK")
    call disp%show( stringSep_SK, deliml = SK_"""" )
    call disp%show("call setSplit(parts_SSK, string_SK, stringSep_SK, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(parts_SSK, string_SK, stringSep_SK, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("parts_SSK")
    call disp%show( parts_SSK, deliml = SK_"""" )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringSep_SK")
    call disp%show( stringSep_SK, deliml = SK_"""" )
    call disp%show("call setSplit(splitIndex, string_SK, stringSep_SK, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(splitIndex, string_SK, stringSep_SK, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Split character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_SK")
    call disp%show( array_SK, deliml = SK_"""" )
    call disp%show("sepArray_SK")
    call disp%show( sepArray_SK, deliml = SK_"""" )
    call disp%show("call setSplit(parts_VSK, array_SK, sepArray_SK)")
                    call setSplit(parts_VSK, array_SK, sepArray_SK)
    call disp%show("parts_VSK")
    call disp%show( parts_VSK, deliml = SK_"""" )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_SK")
    call disp%show( array_SK, deliml = SK_"""" )
    call disp%show("sepArray_SK")
    call disp%show( sepArray_SK, deliml = SK_"""" )
    call disp%show("call setSplit(splitIndex, array_SK, sepArray_SK)")
                    call setSplit(splitIndex, array_SK, sepArray_SK)
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_SK")
    call disp%show( array_SK, deliml = SK_"""" )
    call disp%show("sepArray_SK")
    call disp%show( sepArray_SK, deliml = SK_"""" )
    call disp%show("call setSplit(parts_VSK, array_SK, sepArray_SK, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(parts_VSK, array_SK, sepArray_SK, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("parts_VSK")
    call disp%show( parts_VSK, deliml = SK_"""" )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_SK")
    call disp%show( array_SK, deliml = SK_"""" )
    call disp%show("sepArray_SK")
    call disp%show( sepArray_SK, deliml = SK_"""" )
    call disp%show("call setSplit(splitIndex, array_SK, sepArray_SK, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(splitIndex, array_SK, sepArray_SK, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Split logical array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_LK")
    call disp%show( array_LK )
    call disp%show("sepArray_LK")
    call disp%show( sepArray_LK )
    call disp%show("call setSplit(parts_VLK, array_LK, sepArray_LK)")
                    call setSplit(parts_VLK, array_LK, sepArray_LK)
    call disp%show("parts_VLK")
    call disp%show( parts_VLK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_LK")
    call disp%show( array_LK )
    call disp%show("sepArray_LK")
    call disp%show( sepArray_LK )
    call disp%show("call setSplit(splitIndex, array_LK, sepArray_LK)")
                    call setSplit(splitIndex, array_LK, sepArray_LK)
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_LK")
    call disp%show( array_LK )
    call disp%show("sepArray_LK")
    call disp%show( sepArray_LK )
    call disp%show("call setSplit(parts_VLK, array_LK, sepArray_LK, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(parts_VLK, array_LK, sepArray_LK, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("parts_VLK")
    call disp%show( parts_VLK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_LK")
    call disp%show( array_LK )
    call disp%show("sepArray_LK")
    call disp%show( sepArray_LK )
    call disp%show("call setSplit(splitIndex, array_LK, sepArray_LK, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(splitIndex, array_LK, sepArray_LK, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Split integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_IK")
    call disp%show( array_IK )
    call disp%show("sepArray_IK")
    call disp%show( sepArray_IK )
    call disp%show("call setSplit(parts_VIK, array_IK, sepArray_IK)")
                    call setSplit(parts_VIK, array_IK, sepArray_IK)
    call disp%show("parts_VIK")
    call disp%show( parts_VIK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_IK")
    call disp%show( array_IK )
    call disp%show("sepArray_IK")
    call disp%show( sepArray_IK )
    call disp%show("call setSplit(splitIndex, array_IK, sepArray_IK)")
                    call setSplit(splitIndex, array_IK, sepArray_IK)
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_IK")
    call disp%show( array_IK )
    call disp%show("sepArray_IK")
    call disp%show( sepArray_IK )
    call disp%show("call setSplit(parts_VIK, array_IK, sepArray_IK, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(parts_VIK, array_IK, sepArray_IK, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("parts_VIK")
    call disp%show( parts_VIK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_IK")
    call disp%show( array_IK )
    call disp%show("sepArray_IK")
    call disp%show( sepArray_IK )
    call disp%show("call setSplit(splitIndex, array_IK, sepArray_IK, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(splitIndex, array_IK, sepArray_IK, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Split complex array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_CK")
    call disp%show( array_CK )
    call disp%show("sepArray_CK")
    call disp%show( sepArray_CK )
    call disp%show("call setSplit(parts_VCK, array_CK, sepArray_CK)")
                    call setSplit(parts_VCK, array_CK, sepArray_CK)
    call disp%show("parts_VCK")
    call disp%show( parts_VCK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_CK")
    call disp%show( array_CK )
    call disp%show("sepArray_CK")
    call disp%show( sepArray_CK )
    call disp%show("call setSplit(splitIndex, array_CK, sepArray_CK)")
                    call setSplit(splitIndex, array_CK, sepArray_CK)
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_CK")
    call disp%show( array_CK )
    call disp%show("sepArray_CK")
    call disp%show( sepArray_CK )
    call disp%show("call setSplit(parts_VCK, array_CK, sepArray_CK, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(parts_VCK, array_CK, sepArray_CK, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("parts_VCK")
    call disp%show( parts_VCK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_CK")
    call disp%show( array_CK )
    call disp%show("sepArray_CK")
    call disp%show( sepArray_CK )
    call disp%show("call setSplit(splitIndex, array_CK, sepArray_CK, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(splitIndex, array_CK, sepArray_CK, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%")
    call disp%show("! Split real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_RK")
    call disp%show( array_RK )
    call disp%show("sepArray_RK")
    call disp%show( sepArray_RK )
    call disp%show("call setSplit(parts_VRK, array_RK, sepArray_RK)")
                    call setSplit(parts_VRK, array_RK, sepArray_RK)
    call disp%show("parts_VRK")
    call disp%show( parts_VRK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_RK")
    call disp%show( array_RK )
    call disp%show("sepArray_RK")
    call disp%show( sepArray_RK )
    call disp%show("call setSplit(splitIndex, array_RK, sepArray_RK)")
                    call setSplit(splitIndex, array_RK, sepArray_RK)
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_RK")
    call disp%show( array_RK )
    call disp%show("sepArray_RK")
    call disp%show( sepArray_RK )
    call disp%show("call setSplit(parts_VRK, array_RK, sepArray_RK, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(parts_VRK, array_RK, sepArray_RK, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("parts_VRK")
    call disp%show( parts_VRK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_RK")
    call disp%show( array_RK )
    call disp%show("sepArray_RK")
    call disp%show( sepArray_RK )
    call disp%show("call setSplit(splitIndex, array_RK, sepArray_RK, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(splitIndex, array_RK, sepArray_RK, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Split only particular instances of Sep in array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    string_SK = "AAAAAAAAA"
    array_SK = ["AAAAAAAAA", "AAAAAAAAA", "AAAAAAAAA", "AAAAAAAAA", "AAAAAAAAA", "AAAAAAAAA", "AAAAAAAAA", "AAAAAAAAA", "AAAAAAAAA"]
    array_IK = [integer(IK) :: 0, 1, 0, 2, 3, 0, 4, 5, 0, 0]
    array_RK = [real(RK) :: 0., 1., 0., 2., 3., 0., 4., 5., 0., 0.]
    array_CK = [complex(CK) :: (0., -0.), (1., -1.), (0., -0.), (2., -2.), (3., -3.), (0., -0.), (4., -4.), (5., -5.), (0., -0.), (0., -0.)]
    array_LK = [.false._LK, .true._LK, .false._LK, .true._LK, .true._LK, .false._LK, .true._LK, .true._LK, .false._LK, .false._LK]

    stringSep_SK = "A"
    sepArray_SK = ["AAAAAAAAA"]
    sepArray_IK = [0_IK]
    sepArray_RK = [0._RK]
    sepArray_CK = [(0._CK, -0._CK)]
    sepArray_LK = [.false._LK]

    instance = [-3, 2, -4] ! split at the second occurrence from the beginning and the third occurrence from the end. Duplicate indices are ignored.

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Split character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringSep_SK")
    call disp%show( stringSep_SK, deliml = SK_"""" )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_SSK, string_SK, stringSep_SK, instance = instance)")
                    call setSplit(parts_SSK, string_SK, stringSep_SK, instance = instance)
    call disp%show("parts_SSK")
    call disp%show( parts_SSK, deliml = SK_"""" )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringSep_SK")
    call disp%show( stringSep_SK, deliml = SK_"""" )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, string_SK, stringSep_SK, instance = instance)")
                    call setSplit(splitIndex, string_SK, stringSep_SK, instance = instance)
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringSep_SK")
    call disp%show( stringSep_SK, deliml = SK_"""" )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_SSK, string_SK, stringSep_SK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(parts_SSK, string_SK, stringSep_SK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("parts_SSK")
    call disp%show( parts_SSK, deliml = SK_"""" )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringSep_SK")
    call disp%show( stringSep_SK, deliml = SK_"""" )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, string_SK, stringSep_SK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(splitIndex, string_SK, stringSep_SK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Split character array with vector `sep`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_SK")
    call disp%show( array_SK, deliml = SK_"""" )
    call disp%show("sepArray_SK")
    call disp%show( sepArray_SK, deliml = SK_"""" )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VSK, array_SK, sepArray_SK, instance = instance)")
                    call setSplit(parts_VSK, array_SK, sepArray_SK, instance = instance)
    call disp%show("parts_VSK")
    call disp%show( parts_VSK, deliml = SK_"""" )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_SK")
    call disp%show( array_SK, deliml = SK_"""" )
    call disp%show("sepArray_SK")
    call disp%show( sepArray_SK, deliml = SK_"""" )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_SK, sepArray_SK, instance = instance)")
                    call setSplit(splitIndex, array_SK, sepArray_SK, instance = instance)
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_SK")
    call disp%show( array_SK, deliml = SK_"""" )
    call disp%show("sepArray_SK")
    call disp%show( sepArray_SK, deliml = SK_"""" )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VSK, array_SK, sepArray_SK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(parts_VSK, array_SK, sepArray_SK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("parts_VSK")
    call disp%show( parts_VSK, deliml = SK_"""" )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_SK")
    call disp%show( array_SK, deliml = SK_"""" )
    call disp%show("sepArray_SK")
    call disp%show( sepArray_SK, deliml = SK_"""" )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_SK, sepArray_SK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(splitIndex, array_SK, sepArray_SK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Split character array with scalar `sep`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_SK")
    call disp%show( array_SK, deliml = SK_"""" )
    call disp%show("sepArray_SK(1)")
    call disp%show( sepArray_SK(1), deliml = SK_"""" )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VSK, array_SK, sepArray_SK(1), instance = instance)")
                    call setSplit(parts_VSK, array_SK, sepArray_SK(1), instance = instance)
    call disp%show("parts_VSK")
    call disp%show( parts_VSK, deliml = SK_"""" )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_SK")
    call disp%show( array_SK, deliml = SK_"""" )
    call disp%show("sepArray_SK(1)")
    call disp%show( sepArray_SK(1), deliml = SK_"""" )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_SK, sepArray_SK(1), instance = instance)")
                    call setSplit(splitIndex, array_SK, sepArray_SK(1), instance = instance)
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_SK")
    call disp%show( array_SK, deliml = SK_"""" )
    call disp%show("sepArray_SK(1)")
    call disp%show( sepArray_SK(1), deliml = SK_"""" )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VSK, array_SK, sepArray_SK(1), instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(parts_VSK, array_SK, sepArray_SK(1), instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("parts_VSK")
    call disp%show( parts_VSK, deliml = SK_"""" )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_SK")
    call disp%show( array_SK, deliml = SK_"""" )
    call disp%show("sepArray_SK(1)")
    call disp%show( sepArray_SK(1), deliml = SK_"""" )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_SK, sepArray_SK(1), instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(splitIndex, array_SK, sepArray_SK(1), instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Split logical array with vector `sep`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_LK")
    call disp%show( array_LK )
    call disp%show("sepArray_LK")
    call disp%show( sepArray_LK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VLK, array_LK, sepArray_LK, instance = instance)")
                    call setSplit(parts_VLK, array_LK, sepArray_LK, instance = instance)
    call disp%show("parts_VLK")
    call disp%show( parts_VLK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_LK")
    call disp%show( array_LK )
    call disp%show("sepArray_LK")
    call disp%show( sepArray_LK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_LK, sepArray_LK, instance = instance)")
                    call setSplit(splitIndex, array_LK, sepArray_LK, instance = instance)
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_LK")
    call disp%show( array_LK )
    call disp%show("sepArray_LK")
    call disp%show( sepArray_LK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VLK, array_LK, sepArray_LK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(parts_VLK, array_LK, sepArray_LK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("parts_VLK")
    call disp%show( parts_VLK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_LK")
    call disp%show( array_LK )
    call disp%show("sepArray_LK")
    call disp%show( sepArray_LK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_LK, sepArray_LK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(splitIndex, array_LK, sepArray_LK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Split logical array with scalar `sep`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_LK")
    call disp%show( array_LK )
    call disp%show("sepArray_LK(1)")
    call disp%show( sepArray_LK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VLK, array_LK, sepArray_LK(1), instance = instance)")
                    call setSplit(parts_VLK, array_LK, sepArray_LK(1), instance = instance)
    call disp%show("parts_VLK")
    call disp%show( parts_VLK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_LK")
    call disp%show( array_LK )
    call disp%show("sepArray_LK(1)")
    call disp%show( sepArray_LK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_LK, sepArray_LK(1), instance = instance)")
                    call setSplit(splitIndex, array_LK, sepArray_LK(1), instance = instance)
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_LK")
    call disp%show( array_LK )
    call disp%show("sepArray_LK(1)")
    call disp%show( sepArray_LK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VLK, array_LK, sepArray_LK(1), instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(parts_VLK, array_LK, sepArray_LK(1), instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("parts_VLK")
    call disp%show( parts_VLK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_LK")
    call disp%show( array_LK )
    call disp%show("sepArray_LK(1)")
    call disp%show( sepArray_LK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_LK, sepArray_LK(1), instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(splitIndex, array_LK, sepArray_LK(1), instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Split integer array with vector `sep`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_IK")
    call disp%show( array_IK )
    call disp%show("sepArray_IK")
    call disp%show( sepArray_IK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VIK, array_IK, sepArray_IK, instance = instance)")
                    call setSplit(parts_VIK, array_IK, sepArray_IK, instance = instance)
    call disp%show("parts_VIK")
    call disp%show( parts_VIK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_IK")
    call disp%show( array_IK )
    call disp%show("sepArray_IK")
    call disp%show( sepArray_IK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_IK, sepArray_IK, instance = instance)")
                    call setSplit(splitIndex, array_IK, sepArray_IK, instance = instance)
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_IK")
    call disp%show( array_IK )
    call disp%show("sepArray_IK")
    call disp%show( sepArray_IK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VIK, array_IK, sepArray_IK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(parts_VIK, array_IK, sepArray_IK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("parts_VIK")
    call disp%show( parts_VIK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_IK")
    call disp%show( array_IK )
    call disp%show("sepArray_IK")
    call disp%show( sepArray_IK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_IK, sepArray_IK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(splitIndex, array_IK, sepArray_IK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Split integer array with scalar `sep`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_IK")
    call disp%show( array_IK )
    call disp%show("sepArray_IK(1)")
    call disp%show( sepArray_IK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VIK, array_IK, sepArray_IK(1), instance = instance)")
                    call setSplit(parts_VIK, array_IK, sepArray_IK(1), instance = instance)
    call disp%show("parts_VIK")
    call disp%show( parts_VIK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_IK")
    call disp%show( array_IK )
    call disp%show("sepArray_IK(1)")
    call disp%show( sepArray_IK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_IK, sepArray_IK(1), instance = instance)")
                    call setSplit(splitIndex, array_IK, sepArray_IK(1), instance = instance)
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_IK")
    call disp%show( array_IK )
    call disp%show("sepArray_IK(1)")
    call disp%show( sepArray_IK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VIK, array_IK, sepArray_IK(1), instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(parts_VIK, array_IK, sepArray_IK(1), instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("parts_VIK")
    call disp%show( parts_VIK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_IK")
    call disp%show( array_IK )
    call disp%show("sepArray_IK(1)")
    call disp%show( sepArray_IK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_IK, sepArray_IK(1), instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(splitIndex, array_IK, sepArray_IK(1), instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Split complex array with vector `sep`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_CK")
    call disp%show( array_CK )
    call disp%show("sepArray_CK")
    call disp%show( sepArray_CK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VCK, array_CK, sepArray_CK, instance = instance)")
                    call setSplit(parts_VCK, array_CK, sepArray_CK, instance = instance)
    call disp%show("parts_VCK")
    call disp%show( parts_VCK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_CK")
    call disp%show( array_CK )
    call disp%show("sepArray_CK")
    call disp%show( sepArray_CK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_CK, sepArray_CK, instance = instance)")
                    call setSplit(splitIndex, array_CK, sepArray_CK, instance = instance)
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_CK")
    call disp%show( array_CK )
    call disp%show("sepArray_CK")
    call disp%show( sepArray_CK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VCK, array_CK, sepArray_CK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(parts_VCK, array_CK, sepArray_CK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("parts_VCK")
    call disp%show( parts_VCK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_CK")
    call disp%show( array_CK )
    call disp%show("sepArray_CK")
    call disp%show( sepArray_CK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_CK, sepArray_CK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(splitIndex, array_CK, sepArray_CK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Split complex array with scalar `sep`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_CK")
    call disp%show( array_CK )
    call disp%show("sepArray_CK(1)")
    call disp%show( sepArray_CK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VCK, array_CK, sepArray_CK(1), instance = instance)")
                    call setSplit(parts_VCK, array_CK, sepArray_CK(1), instance = instance)
    call disp%show("parts_VCK")
    call disp%show( parts_VCK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_CK")
    call disp%show( array_CK )
    call disp%show("sepArray_CK(1)")
    call disp%show( sepArray_CK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_CK, sepArray_CK(1), instance = instance)")
                    call setSplit(splitIndex, array_CK, sepArray_CK(1), instance = instance)
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_CK")
    call disp%show( array_CK )
    call disp%show("sepArray_CK(1)")
    call disp%show( sepArray_CK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VCK, array_CK, sepArray_CK(1), instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(parts_VCK, array_CK, sepArray_CK(1), instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("parts_VCK")
    call disp%show( parts_VCK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_CK")
    call disp%show( array_CK )
    call disp%show("sepArray_CK(1)")
    call disp%show( sepArray_CK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_CK, sepArray_CK(1), instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(splitIndex, array_CK, sepArray_CK(1), instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Split real array with vector `sep`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_RK")
    call disp%show( array_RK )
    call disp%show("sepArray_RK")
    call disp%show( sepArray_RK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VRK, array_RK, sepArray_RK, instance = instance)")
                    call setSplit(parts_VRK, array_RK, sepArray_RK, instance = instance)
    call disp%show("parts_VRK")
    call disp%show( parts_VRK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_RK")
    call disp%show( array_RK )
    call disp%show("sepArray_RK")
    call disp%show( sepArray_RK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_RK, sepArray_RK, instance = instance)")
                    call setSplit(splitIndex, array_RK, sepArray_RK, instance = instance)
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_RK")
    call disp%show( array_RK )
    call disp%show("sepArray_RK")
    call disp%show( sepArray_RK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VRK, array_RK, sepArray_RK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(parts_VRK, array_RK, sepArray_RK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("parts_VRK")
    call disp%show( parts_VRK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_RK")
    call disp%show( array_RK )
    call disp%show("sepArray_RK")
    call disp%show( sepArray_RK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_RK, sepArray_RK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(splitIndex, array_RK, sepArray_RK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Split real array with scalar `sep`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_RK")
    call disp%show( array_RK )
    call disp%show("sepArray_RK(1)")
    call disp%show( sepArray_RK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VRK, array_RK, sepArray_RK(1), instance = instance)")
                    call setSplit(parts_VRK, array_RK, sepArray_RK(1), instance = instance)
    call disp%show("parts_VRK")
    call disp%show( parts_VRK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_RK")
    call disp%show( array_RK )
    call disp%show("sepArray_RK(1)")
    call disp%show( sepArray_RK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_RK, sepArray_RK(1), instance = instance)")
                    call setSplit(splitIndex, array_RK, sepArray_RK(1), instance = instance)
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_RK")
    call disp%show( array_RK )
    call disp%show("sepArray_RK(1)")
    call disp%show( sepArray_RK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VRK, array_RK, sepArray_RK(1), instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(parts_VRK, array_RK, sepArray_RK(1), instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("parts_VRK")
    call disp%show( parts_VRK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_RK")
    call disp%show( array_RK )
    call disp%show("sepArray_RK(1)")
    call disp%show( sepArray_RK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_RK, sepArray_RK(1), instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(splitIndex, array_RK, sepArray_RK(1), instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Split at specific instances with a user-defined equivalence test.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Split at case-insensitive instances of vector `sep` within the character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    string_SK = "ABBAbbA"
    stringSep_SK = "bb"

#if PDT_ENABLED
    call disp%skip()
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""")
    call disp%show("stringSep_SK")
    call disp%show( stringSep_SK, deliml = SK_"""")
    call disp%show("call setSplit(parts_SSK, string_SK, stringSep_SK, iseq = iseq_SK)")
                    call setSplit(parts_SSK, string_SK, stringSep_SK, iseq = iseq_SK)
    call disp%show("parts_SSK")
    call disp%show( parts_SSK, deliml = SK_"""")
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""")
    call disp%show("stringSep_SK")
    call disp%show( stringSep_SK, deliml = SK_"""")
    call disp%show("call setSplit(splitIndex, string_SK, stringSep_SK, iseq = iseq_SK)")
                    call setSplit(splitIndex, string_SK, stringSep_SK, iseq = iseq_SK)
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""")
    call disp%show("stringSep_SK")
    call disp%show( stringSep_SK, deliml = SK_"""")
    call disp%show("call setSplit(parts_SSK, string_SK, stringSep_SK, iseq = iseq_SK, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(parts_SSK, string_SK, stringSep_SK, iseq = iseq_SK, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("parts_SSK")
    call disp%show( parts_SSK, deliml = SK_"""")
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""")
    call disp%show("stringSep_SK")
    call disp%show( stringSep_SK, deliml = SK_"""")
    call disp%show("call setSplit(splitIndex, string_SK, stringSep_SK, iseq = iseq_SK, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(splitIndex, string_SK, stringSep_SK, iseq = iseq_SK, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Split at specific instances of vector `sep` within the real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    array_RK = [0._RK, 1.01_RK, 1.04_RK, 0.98_RK, 1.0_RK, 1.02_RK]
    sepArray_RK = [-1._RK, 1._RK]
    instance = [-2_IK]

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_RK")
    call disp%show( array_RK )
    call disp%show("sepArray_RK")
    call disp%show( sepArray_RK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VRK, array_RK, sepArray_RK, iseq = iseq_vec_RK, instance = instance)")
                    call setSplit(parts_VRK, array_RK, sepArray_RK, iseq = iseq_vec_RK, instance = instance)
    call disp%show("parts_VRK")
    call disp%show( parts_VRK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_RK")
    call disp%show( array_RK )
    call disp%show("sepArray_RK")
    call disp%show( sepArray_RK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_RK, sepArray_RK, iseq = iseq_vec_RK, instance = instance)")
                    call setSplit(splitIndex, array_RK, sepArray_RK, iseq = iseq_vec_RK, instance = instance)
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_RK")
    call disp%show( array_RK )
    call disp%show("sepArray_RK")
    call disp%show( sepArray_RK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VRK, array_RK, sepArray_RK, iseq = iseq_vec_RK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(parts_VRK, array_RK, sepArray_RK, iseq = iseq_vec_RK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("parts_VRK")
    call disp%show( parts_VRK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_RK")
    call disp%show( array_RK )
    call disp%show("sepArray_RK")
    call disp%show( sepArray_RK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_RK, sepArray_RK, iseq = iseq_vec_RK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(splitIndex, array_RK, sepArray_RK, iseq = iseq_vec_RK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Split at specific instances of scalar `sep` within the real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_RK")
    call disp%show( array_RK )
    call disp%show("sepArray_RK(1)")
    call disp%show( sepArray_RK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VRK, array_RK, sepArray_RK(1), iseq = iseq_RK, instance = instance)")
                    call setSplit(parts_VRK, array_RK, sepArray_RK(1), iseq = iseq_RK, instance = instance)
    call disp%show("parts_VRK")
    call disp%show( parts_VRK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_RK")
    call disp%show( array_RK )
    call disp%show("sepArray_RK(1)")
    call disp%show( sepArray_RK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_RK, sepArray_RK(1), iseq = iseq_RK, instance = instance)")
                    call setSplit(splitIndex, array_RK, sepArray_RK(1), iseq = iseq_RK, instance = instance)
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()

#if PDT_ENABLED
    call disp%skip()
    call disp%show("Array_RK")
    call disp%show( array_RK )
    call disp%show("sepArray_RK(1)")
    call disp%show( sepArray_RK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(parts_VRK, array_RK, sepArray_RK(1), iseq = iseq_RK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(parts_VRK, array_RK, sepArray_RK(1), iseq = iseq_RK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("parts_VRK")
    call disp%show( parts_VRK )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("Array_RK")
    call disp%show( array_RK )
    call disp%show("sepArray_RK(1)")
    call disp%show( sepArray_RK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setSplit(splitIndex, array_RK, sepArray_RK(1), iseq = iseq_RK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.")
                    call setSplit(splitIndex, array_RK, sepArray_RK(1), iseq = iseq_RK, instance = instance, keep = .true._LK) ! Keep the `sep` cases in the output.
    call disp%show("splitIndex")
    call disp%show( splitIndex )
    call disp%skip()

contains

    pure function iseq_SK(ArraySegment, Sep) result(equivalent)
        use pm_strASCII, only: getStrLower
        character(*, SK), intent(in)    :: Sep, ArraySegment
        logical(LK)                     :: equivalent
        equivalent = Sep == getStrLower(ArraySegment)
    end function

    function iseq_RK(arraySegment, sep) result(equivalent)
        real(RK)    , intent(in)    :: sep, arraySegment
        logical(LK)                 :: equivalent
        equivalent = abs(abs(sep) - abs(arraySegment)) < 0.05_RK
    end function

    function iseq_vec_RK(ArraySegment, Sep, lenSep) result(equivalent)
        integer(IK) , intent(in)    :: lenSep
        real(RK)    , intent(in)    :: Sep(lenSep), ArraySegment(lenSep)
        logical(LK)                 :: equivalent
        equivalent = all(abs(abs(Sep) - abs(ArraySegment)) < 0.05_RK)
    end function

end program example