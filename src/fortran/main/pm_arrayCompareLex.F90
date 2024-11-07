!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief
!>  This module contains procedures and generic interfaces for performing
!>  lexicographic comparisons of two arrays of similar type, kind, and rank.<br>
!>
!>  \details
!>  A **lexicographical comparison** is the kind of comparison generally used to sort words alphabetically in dictionaries.<br>
!>  It involves comparing sequentially the elements that have the same position in
!>  both ranges against each other until one element is not equivalent to the other.<br>
!>  The result of comparing these first non-matching elements is the result of the lexicographical comparison.<br>
!>  If the two arrays are of different lengths but are equivalent up to the last common element comparison,
!>  then the result of the comparison of the sizes of the two arrays will be returned.
!>
!>  \note
!>  The **lexicographical order** is also known as the **Dictionary order**
!>  because the words in a dictionary are commonly sorted lexicographically.<br>
!>
!>  \note
!>  Note that the Fortran language intrinsically defines the lexicographic comparison of scalar character objects.<br>
!>  However, the intrinsic Fortran definition pads the smaller array with blanks before comparison.<br>
!>  Therefore, as long as the extra elements of the longer array correspond to characters that appear after
!>  the blank character in the **collating sequence** of the processor, the intrinsic Fortran scalar character
!>  comparison matches the above definition for lexicographic comparison that is used in this module.<br>
!>
!>  \note
!>  Note that the result of a lexicographic comparison is always a scalar regardless of the ranks of the objects being compared.<br>
!>
!>  \note
!>  By definition, a lexicographic comparison can only be defined for objects of the same rank.<br>
!>
!>  \note
!>  By default, a `.true.` value is assumed to be larger than the `.false.` value in this module.<br>
!>  Therefore, `.false. < .true.` evaluates to `.true.`.<br>
!>
!>  \note
!>  By default, the two components of complex numbers are also compared lexically in this module.<br>
!>  For example, `(1.,2.) < (1.,1.)` evaluates to `.false.` while `(1.,2.) < (2.,1.)` evaluates to `.true.`.<br>
!>
!>  \test
!>  [test_pm_arrayCompareLex](@ref test_pm_arrayCompareLex)<br>
!>
!>  \todo
!>  \pmed A generic interface for supplying a user-defined arbitrary comparison operator should be added.<br>
!>
!>  \todo
!>  \pmed Ideally, a procedure with unlimited polymorphic arguments could be also added for lexical
!>  comparison of all types of input arrays, albeit with a user-defined operator supplied to the procedure.<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayCompareLex

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_arrayCompareLex"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_arrayCompareLex_isllt
    !>  Generate and return the result of the lexicographic comparison of two
    !>  input objects of the same type, kind, and rank using the `<` operator.
    !>
    !>  \details
    !>  For more information about the lexicographic comparisons defined in this generic procedure, in particular,
    !>  in the context of `complex` and `logical` arrays, see the detailed description of [pm_arrayCompareLex](@ref pm_arrayCompareLex).
    !>
    !>  \param[in]  array1  :   The input `contiguous` **array of rank `1`** of either<br>
    !>                          <ul>
    !>                              <li>    type `character` of kind \SKALL or, <br>
    !>                              <li>    type `integer` of kind \IKALL or, <br>
    !>                              <li>    type `logical` of kind \LKALL or, <br>
    !>                              <li>    type `complex` of kind \CKALL or, <br>
    !>                              <li>    type `real` of kind \RKALL or, <br>
    !>                          </ul>
    !>                          or,
    !>                          <ul>
    !>                              <li>    a **scalar assumed-length** `character` of kind \SKALL, <br>
    !>                          </ul>
    !>                          whose contents will be lexicographically compared with the contents of `array2`.
    !>  \param[in]  array2  :   The input of the same type, kind, and rank as `array1`.
    !>
    !>  \return
    !>  `compares`          :   The output scalar of type `logical` of default kind \LK containing the
    !>                          result of the lexical comparison of the values of the two input objects.
    !>
    !>  \interface{llt}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_arrayCompareLex, only: operator(.llt.)
    !>      logical(LK) :: compares
    !>
    !>      compares = array1 .llt. array2
    !>
    !>  \endcode
    !>
    !>  \remark
    !>  The operator `.llt.` stands for **lexically less than**.
    !>
    !>  \pure
    !>
    !>  \see
    !>  [operator(.llt.)](@ref pm_arrayCompareLex_isllt)<br>
    !>  [operator(.lle.)](@ref pm_arrayCompareLex_islle)<br>
    !>  [operator(.lge.)](@ref pm_arrayCompareLex_islge)<br>
    !>  [operator(.lgt.)](@ref pm_arrayCompareLex_islgt)<br>
    !>
    !>  \example{llt}
    !>  \include{lineno} example/pm_arrayCompareLex/llt/main.F90
    !>  \compilef{llt}
    !>  \output{llt}
    !>  \include{lineno} example/pm_arrayCompareLex/llt/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayCompareLex](@ref test_pm_arrayCompareLex)
    !>
    !>  \todo
    !>  \pvlow The functionality of this generic interface should be extended to input container arguments, also to arrays of higher ranks.
    !>
    !>  \final{llt}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(.llt.)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function isllt_D0_D0_SK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D0_D0_SK5
#endif
        use pm_kind, only: LK, SKG => SK5
        character(*,SKG), intent(in)                :: array1
        character(*,SKG), intent(in)                :: array2
        logical(LK)                                 :: compares
    end function
#endif

#if SK4_ENABLED
    pure module function isllt_D0_D0_SK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D0_D0_SK4
#endif
        use pm_kind, only: LK, SKG => SK4
        character(*,SKG), intent(in)                :: array1
        character(*,SKG), intent(in)                :: array2
        logical(LK)                                 :: compares
    end function
#endif

#if SK3_ENABLED
    pure module function isllt_D0_D0_SK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D0_D0_SK3
#endif
        use pm_kind, only: LK, SKG => SK3
        character(*,SKG), intent(in)                :: array1
        character(*,SKG), intent(in)                :: array2
        logical(LK)                                 :: compares
    end function
#endif

#if SK2_ENABLED
    pure module function isllt_D0_D0_SK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D0_D0_SK2
#endif
        use pm_kind, only: LK, SKG => SK2
        character(*,SKG), intent(in)                :: array1
        character(*,SKG), intent(in)                :: array2
        logical(LK)                                 :: compares
    end function
#endif

#if SK1_ENABLED
    pure module function isllt_D0_D0_SK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D0_D0_SK1
#endif
        use pm_kind, only: LK, SKG => SK1
        character(*,SKG), intent(in)                :: array1
        character(*,SKG), intent(in)                :: array2
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function isllt_D1_D1_SK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_SK5
#endif
        use pm_kind, only: LK, SKG => SK5
        character(*,SKG), intent(in), contiguous    :: array1(:)
        character(*,SKG), intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if SK4_ENABLED
    pure module function isllt_D1_D1_SK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_SK4
#endif
        use pm_kind, only: LK, SKG => SK4
        character(*,SKG), intent(in), contiguous    :: array1(:)
        character(*,SKG), intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if SK3_ENABLED
    pure module function isllt_D1_D1_SK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_SK3
#endif
        use pm_kind, only: LK, SKG => SK3
        character(*,SKG), intent(in), contiguous    :: array1(:)
        character(*,SKG), intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if SK2_ENABLED
    pure module function isllt_D1_D1_SK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_SK2
#endif
        use pm_kind, only: LK, SKG => SK2
        character(*,SKG), intent(in), contiguous    :: array1(:)
        character(*,SKG), intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if SK1_ENABLED
    pure module function isllt_D1_D1_SK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_SK1
#endif
        use pm_kind, only: LK, SKG => SK1
        character(*,SKG), intent(in), contiguous    :: array1(:)
        character(*,SKG), intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module function isllt_D1_D1_IK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_IK5
#endif
        use pm_kind, only: LK, IKG => IK5
        integer(IKG)    , intent(in), contiguous    :: array1(:)
        integer(IKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if IK4_ENABLED
    pure module function isllt_D1_D1_IK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_IK4
#endif
        use pm_kind, only: LK, IKG => IK4
        integer(IKG)    , intent(in), contiguous    :: array1(:)
        integer(IKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if IK3_ENABLED
    pure module function isllt_D1_D1_IK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_IK3
#endif
        use pm_kind, only: LK, IKG => IK3
        integer(IKG)    , intent(in), contiguous    :: array1(:)
        integer(IKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if IK2_ENABLED
    pure module function isllt_D1_D1_IK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_IK2
#endif
        use pm_kind, only: LK, IKG => IK2
        integer(IKG)    , intent(in), contiguous    :: array1(:)
        integer(IKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if IK1_ENABLED
    pure module function isllt_D1_D1_IK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_IK1
#endif
        use pm_kind, only: LK, IKG => IK1
        integer(IKG)    , intent(in), contiguous    :: array1(:)
        integer(IKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module function isllt_D1_D1_LK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_LK5
#endif
        use pm_kind, only: LK, LKG => LK5
        logical(LKG)    , intent(in), contiguous    :: array1(:)
        logical(LKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if LK4_ENABLED
    pure module function isllt_D1_D1_LK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_LK4
#endif
        use pm_kind, only: LK, LKG => LK4
        logical(LKG)    , intent(in), contiguous    :: array1(:)
        logical(LKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if LK3_ENABLED
    pure module function isllt_D1_D1_LK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_LK3
#endif
        use pm_kind, only: LK, LKG => LK3
        logical(LKG)    , intent(in), contiguous    :: array1(:)
        logical(LKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if LK2_ENABLED
    pure module function isllt_D1_D1_LK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_LK2
#endif
        use pm_kind, only: LK, LKG => LK2
        logical(LKG)    , intent(in), contiguous    :: array1(:)
        logical(LKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if LK1_ENABLED
    pure module function isllt_D1_D1_LK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_LK1
#endif
        use pm_kind, only: LK, LKG => LK1
        logical(LKG)    , intent(in), contiguous    :: array1(:)
        logical(LKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function isllt_D1_D1_CK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_CK5
#endif
        use pm_kind, only: LK, CKG => CK5
        complex(CKG)    , intent(in), contiguous    :: array1(:)
        complex(CKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if CK4_ENABLED
    pure module function isllt_D1_D1_CK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_CK4
#endif
        use pm_kind, only: LK, CKG => CK4
        complex(CKG)    , intent(in), contiguous    :: array1(:)
        complex(CKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if CK3_ENABLED
    pure module function isllt_D1_D1_CK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_CK3
#endif
        use pm_kind, only: LK, CKG => CK3
        complex(CKG)    , intent(in), contiguous    :: array1(:)
        complex(CKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if CK2_ENABLED
    pure module function isllt_D1_D1_CK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_CK2
#endif
        use pm_kind, only: LK, CKG => CK2
        complex(CKG)    , intent(in), contiguous    :: array1(:)
        complex(CKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if CK1_ENABLED
    pure module function isllt_D1_D1_CK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_CK1
#endif
        use pm_kind, only: LK, CKG => CK1
        complex(CKG)    , intent(in), contiguous    :: array1(:)
        complex(CKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module function isllt_D1_D1_RK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_RK5
#endif
        use pm_kind, only: LK, RKG => RK5
        real(RKG)       , intent(in), contiguous    :: array1(:)
        real(RKG)       , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if RK4_ENABLED
    pure module function isllt_D1_D1_RK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_RK4
#endif
        use pm_kind, only: LK, RKG => RK4
        real(RKG)       , intent(in), contiguous    :: array1(:)
        real(RKG)       , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if RK3_ENABLED
    pure module function isllt_D1_D1_RK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_RK3
#endif
        use pm_kind, only: LK, RKG => RK3
        real(RKG)       , intent(in), contiguous    :: array1(:)
        real(RKG)       , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if RK2_ENABLED
    pure module function isllt_D1_D1_RK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_RK2
#endif
        use pm_kind, only: LK, RKG => RK2
        real(RKG)       , intent(in), contiguous    :: array1(:)
        real(RKG)       , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if RK1_ENABLED
    pure module function isllt_D1_D1_RK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isllt_D1_D1_RK1
#endif
        use pm_kind, only: LK, RKG => RK1
        real(RKG)       , intent(in), contiguous    :: array1(:)
        real(RKG)       , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_arrayCompareLex_islle
    !>  Generate and return the result of the lexicographic comparison of two
    !>  input objects of the same type, kind, and rank using the `<=` operator.
    !>
    !>  \details
    !>  For more information about the lexicographic comparisons defined in this generic procedure, in particular,
    !>  in the context of `complex` and `logical` arrays, see the detailed description of [pm_arrayCompareLex](@ref pm_arrayCompareLex).
    !>
    !>  \param[in]  array1  :   The input `contiguous` **array of rank `1`** of either<br>
    !>                          <ul>
    !>                              <li>    type `character` of kind \SKALL or, <br>
    !>                              <li>    type `integer` of kind \IKALL or, <br>
    !>                              <li>    type `logical` of kind \LKALL or, <br>
    !>                              <li>    type `complex` of kind \CKALL or, <br>
    !>                              <li>    type `real` of kind \RKALL or, <br>
    !>                          </ul>
    !>                          or,
    !>                          <ul>
    !>                              <li>    a **scalar assumed-length** `character` of kind \SKALL, <br>
    !>                          </ul>
    !>                          whose contents will be lexicographically compared with the contents of `array2`.
    !>  \param[in]  array2  :   The input of the same type, kind, and rank as `array1`.
    !>
    !>  \return
    !>  `compares`          :   The output scalar of type `logical` of default kind \LK containing the
    !>                          result of the lexical comparison of the values of the two input objects.
    !>
    !>  \interface{lle}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_arrayCompareLex, only: operator(.lle.)
    !>      logical(LK) :: compares
    !>
    !>      compares = array1 .lle. array2
    !>
    !>  \endcode
    !>
    !>  \remark
    !>  The operator `.lle.` stands for **lexically less than or equal to**.
    !>
    !>  \pure
    !>
    !>  \see
    !>  [operator(.llt.)](@ref pm_arrayCompareLex_isllt)<br>
    !>  [operator(.lle.)](@ref pm_arrayCompareLex_islle)<br>
    !>  [operator(.lge.)](@ref pm_arrayCompareLex_islge)<br>
    !>  [operator(.lgt.)](@ref pm_arrayCompareLex_islgt)<br>
    !>
    !>  \example{lle}
    !>  \include{lineno} example/pm_arrayCompareLex/lle/main.F90
    !>  \compilef{lle}
    !>  \output{lle}
    !>  \include{lineno} example/pm_arrayCompareLex/lle/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayCompareLex](@ref test_pm_arrayCompareLex)
    !>
    !>  \todo
    !>  \pvlow The functionality of this generic interface should be extended to input container arguments, also to arrays of higher ranks.
    !>
    !>  \final{lle}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(.lle.)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function islle_D0_D0_SK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D0_D0_SK5
#endif
        use pm_kind, only: LK, SKG => SK5
        character(*,SKG), intent(in)                :: array1
        character(*,SKG), intent(in)                :: array2
        logical(LK)                                 :: compares
    end function
#endif

#if SK4_ENABLED
    pure module function islle_D0_D0_SK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D0_D0_SK4
#endif
        use pm_kind, only: LK, SKG => SK4
        character(*,SKG), intent(in)                :: array1
        character(*,SKG), intent(in)                :: array2
        logical(LK)                                 :: compares
    end function
#endif

#if SK3_ENABLED
    pure module function islle_D0_D0_SK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D0_D0_SK3
#endif
        use pm_kind, only: LK, SKG => SK3
        character(*,SKG), intent(in)                :: array1
        character(*,SKG), intent(in)                :: array2
        logical(LK)                                 :: compares
    end function
#endif

#if SK2_ENABLED
    pure module function islle_D0_D0_SK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D0_D0_SK2
#endif
        use pm_kind, only: LK, SKG => SK2
        character(*,SKG), intent(in)                :: array1
        character(*,SKG), intent(in)                :: array2
        logical(LK)                                 :: compares
    end function
#endif

#if SK1_ENABLED
    pure module function islle_D0_D0_SK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D0_D0_SK1
#endif
        use pm_kind, only: LK, SKG => SK1
        character(*,SKG), intent(in)                :: array1
        character(*,SKG), intent(in)                :: array2
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function islle_D1_D1_SK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_SK5
#endif
        use pm_kind, only: LK, SKG => SK5
        character(*,SKG), intent(in), contiguous    :: array1(:)
        character(*,SKG), intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if SK4_ENABLED
    pure module function islle_D1_D1_SK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_SK4
#endif
        use pm_kind, only: LK, SKG => SK4
        character(*,SKG), intent(in), contiguous    :: array1(:)
        character(*,SKG), intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if SK3_ENABLED
    pure module function islle_D1_D1_SK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_SK3
#endif
        use pm_kind, only: LK, SKG => SK3
        character(*,SKG), intent(in), contiguous    :: array1(:)
        character(*,SKG), intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if SK2_ENABLED
    pure module function islle_D1_D1_SK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_SK2
#endif
        use pm_kind, only: LK, SKG => SK2
        character(*,SKG), intent(in), contiguous    :: array1(:)
        character(*,SKG), intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if SK1_ENABLED
    pure module function islle_D1_D1_SK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_SK1
#endif
        use pm_kind, only: LK, SKG => SK1
        character(*,SKG), intent(in), contiguous    :: array1(:)
        character(*,SKG), intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module function islle_D1_D1_IK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_IK5
#endif
        use pm_kind, only: LK, IKG => IK5
        integer(IKG)    , intent(in), contiguous    :: array1(:)
        integer(IKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if IK4_ENABLED
    pure module function islle_D1_D1_IK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_IK4
#endif
        use pm_kind, only: LK, IKG => IK4
        integer(IKG)    , intent(in), contiguous    :: array1(:)
        integer(IKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if IK3_ENABLED
    pure module function islle_D1_D1_IK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_IK3
#endif
        use pm_kind, only: LK, IKG => IK3
        integer(IKG)    , intent(in), contiguous    :: array1(:)
        integer(IKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if IK2_ENABLED
    pure module function islle_D1_D1_IK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_IK2
#endif
        use pm_kind, only: LK, IKG => IK2
        integer(IKG)    , intent(in), contiguous    :: array1(:)
        integer(IKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if IK1_ENABLED
    pure module function islle_D1_D1_IK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_IK1
#endif
        use pm_kind, only: LK, IKG => IK1
        integer(IKG)    , intent(in), contiguous    :: array1(:)
        integer(IKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module function islle_D1_D1_LK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_LK5
#endif
        use pm_kind, only: LK, LKG => LK5
        logical(LKG)    , intent(in), contiguous    :: array1(:)
        logical(LKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if LK4_ENABLED
    pure module function islle_D1_D1_LK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_LK4
#endif
        use pm_kind, only: LK, LKG => LK4
        logical(LKG)    , intent(in), contiguous    :: array1(:)
        logical(LKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if LK3_ENABLED
    pure module function islle_D1_D1_LK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_LK3
#endif
        use pm_kind, only: LK, LKG => LK3
        logical(LKG)    , intent(in), contiguous    :: array1(:)
        logical(LKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if LK2_ENABLED
    pure module function islle_D1_D1_LK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_LK2
#endif
        use pm_kind, only: LK, LKG => LK2
        logical(LKG)    , intent(in), contiguous    :: array1(:)
        logical(LKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if LK1_ENABLED
    pure module function islle_D1_D1_LK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_LK1
#endif
        use pm_kind, only: LK, LKG => LK1
        logical(LKG)    , intent(in), contiguous    :: array1(:)
        logical(LKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function islle_D1_D1_CK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_CK5
#endif
        use pm_kind, only: LK, CKG => CK5
        complex(CKG)    , intent(in), contiguous    :: array1(:)
        complex(CKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if CK4_ENABLED
    pure module function islle_D1_D1_CK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_CK4
#endif
        use pm_kind, only: LK, CKG => CK4
        complex(CKG)    , intent(in), contiguous    :: array1(:)
        complex(CKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if CK3_ENABLED
    pure module function islle_D1_D1_CK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_CK3
#endif
        use pm_kind, only: LK, CKG => CK3
        complex(CKG)    , intent(in), contiguous    :: array1(:)
        complex(CKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if CK2_ENABLED
    pure module function islle_D1_D1_CK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_CK2
#endif
        use pm_kind, only: LK, CKG => CK2
        complex(CKG)    , intent(in), contiguous    :: array1(:)
        complex(CKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if CK1_ENABLED
    pure module function islle_D1_D1_CK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_CK1
#endif
        use pm_kind, only: LK, CKG => CK1
        complex(CKG)    , intent(in), contiguous    :: array1(:)
        complex(CKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module function islle_D1_D1_RK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_RK5
#endif
        use pm_kind, only: LK, RKG => RK5
        real(RKG)       , intent(in), contiguous    :: array1(:)
        real(RKG)       , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if RK4_ENABLED
    pure module function islle_D1_D1_RK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_RK4
#endif
        use pm_kind, only: LK, RKG => RK4
        real(RKG)       , intent(in), contiguous    :: array1(:)
        real(RKG)       , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if RK3_ENABLED
    pure module function islle_D1_D1_RK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_RK3
#endif
        use pm_kind, only: LK, RKG => RK3
        real(RKG)       , intent(in), contiguous    :: array1(:)
        real(RKG)       , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if RK2_ENABLED
    pure module function islle_D1_D1_RK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_RK2
#endif
        use pm_kind, only: LK, RKG => RK2
        real(RKG)       , intent(in), contiguous    :: array1(:)
        real(RKG)       , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if RK1_ENABLED
    pure module function islle_D1_D1_RK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islle_D1_D1_RK1
#endif
        use pm_kind, only: LK, RKG => RK1
        real(RKG)       , intent(in), contiguous    :: array1(:)
        real(RKG)       , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_arrayCompareLex_islge
    !>  Generate and return the result of the lexicographic comparison of two
    !>  input objects of the same type, kind, and rank using the `>=` operator.
    !>
    !>  \details
    !>  For more information about the lexicographic comparisons defined in this generic procedure, in particular,
    !>  in the context of `complex` and `logical` arrays, see the detailed description of [pm_arrayCompareLex](@ref pm_arrayCompareLex).
    !>
    !>  \param[in]  array1  :   The input `contiguous` **array of rank `1`** of either<br>
    !>                          <ul>
    !>                              <li>    type `character` of kind \SKALL or, <br>
    !>                              <li>    type `integer` of kind \IKALL or, <br>
    !>                              <li>    type `logical` of kind \LKALL or, <br>
    !>                              <li>    type `complex` of kind \CKALL or, <br>
    !>                              <li>    type `real` of kind \RKALL or, <br>
    !>                          </ul>
    !>                          or,
    !>                          <ul>
    !>                              <li>    a **scalar assumed-length** `character` of kind \SKALL, <br>
    !>                          </ul>
    !>                          whose contents will be lexicographically compared with the contents of `array2`.
    !>  \param[in]  array2  :   The input of the same type, kind, and rank as `array1`.
    !>
    !>  \return
    !>  `compares`          :   The output scalar of type `logical` of default kind \LK containing the
    !>                          result of the lexical comparison of the values of the two input objects.
    !>
    !>  \interface{lge}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_arrayCompareLex, only: operator(.lge.)
    !>      logical(LK) :: compares
    !>
    !>      compares = array1 .lge. array2
    !>
    !>  \endcode
    !>
    !>  \remark
    !>  The operator `.lge.` stands for **lexically greater than or equal to**.
    !>
    !>  \pure
    !>
    !>  \see
    !>  [operator(.llt.)](@ref pm_arrayCompareLex_isllt)<br>
    !>  [operator(.lle.)](@ref pm_arrayCompareLex_islle)<br>
    !>  [operator(.lge.)](@ref pm_arrayCompareLex_islge)<br>
    !>  [operator(.lgt.)](@ref pm_arrayCompareLex_islgt)<br>
    !>
    !>  \example{lge}
    !>  \include{lineno} example/pm_arrayCompareLex/lge/main.F90
    !>  \compilef{lge}
    !>  \output{lge}
    !>  \include{lineno} example/pm_arrayCompareLex/lge/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayCompareLex](@ref test_pm_arrayCompareLex)
    !>
    !>  \todo
    !>  \pvlow The functionality of this generic interface should be extended to input container arguments, also to arrays of higher ranks.
    !>
    !>  \final{lge}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(.lge.)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function islge_D0_D0_SK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D0_D0_SK5
#endif
        use pm_kind, only: LK, SKG => SK5
        character(*,SKG), intent(in)                :: array1
        character(*,SKG), intent(in)                :: array2
        logical(LK)                                 :: compares
    end function
#endif

#if SK4_ENABLED
    pure module function islge_D0_D0_SK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D0_D0_SK4
#endif
        use pm_kind, only: LK, SKG => SK4
        character(*,SKG), intent(in)                :: array1
        character(*,SKG), intent(in)                :: array2
        logical(LK)                                 :: compares
    end function
#endif

#if SK3_ENABLED
    pure module function islge_D0_D0_SK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D0_D0_SK3
#endif
        use pm_kind, only: LK, SKG => SK3
        character(*,SKG), intent(in)                :: array1
        character(*,SKG), intent(in)                :: array2
        logical(LK)                                 :: compares
    end function
#endif

#if SK2_ENABLED
    pure module function islge_D0_D0_SK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D0_D0_SK2
#endif
        use pm_kind, only: LK, SKG => SK2
        character(*,SKG), intent(in)                :: array1
        character(*,SKG), intent(in)                :: array2
        logical(LK)                                 :: compares
    end function
#endif

#if SK1_ENABLED
    pure module function islge_D0_D0_SK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D0_D0_SK1
#endif
        use pm_kind, only: LK, SKG => SK1
        character(*,SKG), intent(in)                :: array1
        character(*,SKG), intent(in)                :: array2
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function islge_D1_D1_SK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_SK5
#endif
        use pm_kind, only: LK, SKG => SK5
        character(*,SKG), intent(in), contiguous    :: array1(:)
        character(*,SKG), intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if SK4_ENABLED
    pure module function islge_D1_D1_SK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_SK4
#endif
        use pm_kind, only: LK, SKG => SK4
        character(*,SKG), intent(in), contiguous    :: array1(:)
        character(*,SKG), intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if SK3_ENABLED
    pure module function islge_D1_D1_SK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_SK3
#endif
        use pm_kind, only: LK, SKG => SK3
        character(*,SKG), intent(in), contiguous    :: array1(:)
        character(*,SKG), intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if SK2_ENABLED
    pure module function islge_D1_D1_SK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_SK2
#endif
        use pm_kind, only: LK, SKG => SK2
        character(*,SKG), intent(in), contiguous    :: array1(:)
        character(*,SKG), intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if SK1_ENABLED
    pure module function islge_D1_D1_SK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_SK1
#endif
        use pm_kind, only: LK, SKG => SK1
        character(*,SKG), intent(in), contiguous    :: array1(:)
        character(*,SKG), intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module function islge_D1_D1_IK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_IK5
#endif
        use pm_kind, only: LK, IKG => IK5
        integer(IKG)    , intent(in), contiguous    :: array1(:)
        integer(IKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if IK4_ENABLED
    pure module function islge_D1_D1_IK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_IK4
#endif
        use pm_kind, only: LK, IKG => IK4
        integer(IKG)    , intent(in), contiguous    :: array1(:)
        integer(IKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if IK3_ENABLED
    pure module function islge_D1_D1_IK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_IK3
#endif
        use pm_kind, only: LK, IKG => IK3
        integer(IKG)    , intent(in), contiguous    :: array1(:)
        integer(IKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if IK2_ENABLED
    pure module function islge_D1_D1_IK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_IK2
#endif
        use pm_kind, only: LK, IKG => IK2
        integer(IKG)    , intent(in), contiguous    :: array1(:)
        integer(IKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if IK1_ENABLED
    pure module function islge_D1_D1_IK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_IK1
#endif
        use pm_kind, only: LK, IKG => IK1
        integer(IKG)    , intent(in), contiguous    :: array1(:)
        integer(IKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module function islge_D1_D1_LK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_LK5
#endif
        use pm_kind, only: LK, LKG => LK5
        logical(LKG)    , intent(in), contiguous    :: array1(:)
        logical(LKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if LK4_ENABLED
    pure module function islge_D1_D1_LK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_LK4
#endif
        use pm_kind, only: LK, LKG => LK4
        logical(LKG)    , intent(in), contiguous    :: array1(:)
        logical(LKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if LK3_ENABLED
    pure module function islge_D1_D1_LK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_LK3
#endif
        use pm_kind, only: LK, LKG => LK3
        logical(LKG)    , intent(in), contiguous    :: array1(:)
        logical(LKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if LK2_ENABLED
    pure module function islge_D1_D1_LK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_LK2
#endif
        use pm_kind, only: LK, LKG => LK2
        logical(LKG)    , intent(in), contiguous    :: array1(:)
        logical(LKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if LK1_ENABLED
    pure module function islge_D1_D1_LK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_LK1
#endif
        use pm_kind, only: LK, LKG => LK1
        logical(LKG)    , intent(in), contiguous    :: array1(:)
        logical(LKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function islge_D1_D1_CK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_CK5
#endif
        use pm_kind, only: LK, CKG => CK5
        complex(CKG)    , intent(in), contiguous    :: array1(:)
        complex(CKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if CK4_ENABLED
    pure module function islge_D1_D1_CK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_CK4
#endif
        use pm_kind, only: LK, CKG => CK4
        complex(CKG)    , intent(in), contiguous    :: array1(:)
        complex(CKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if CK3_ENABLED
    pure module function islge_D1_D1_CK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_CK3
#endif
        use pm_kind, only: LK, CKG => CK3
        complex(CKG)    , intent(in), contiguous    :: array1(:)
        complex(CKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if CK2_ENABLED
    pure module function islge_D1_D1_CK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_CK2
#endif
        use pm_kind, only: LK, CKG => CK2
        complex(CKG)    , intent(in), contiguous    :: array1(:)
        complex(CKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if CK1_ENABLED
    pure module function islge_D1_D1_CK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_CK1
#endif
        use pm_kind, only: LK, CKG => CK1
        complex(CKG)    , intent(in), contiguous    :: array1(:)
        complex(CKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module function islge_D1_D1_RK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_RK5
#endif
        use pm_kind, only: LK, RKG => RK5
        real(RKG)       , intent(in), contiguous    :: array1(:)
        real(RKG)       , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if RK4_ENABLED
    pure module function islge_D1_D1_RK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_RK4
#endif
        use pm_kind, only: LK, RKG => RK4
        real(RKG)       , intent(in), contiguous    :: array1(:)
        real(RKG)       , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if RK3_ENABLED
    pure module function islge_D1_D1_RK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_RK3
#endif
        use pm_kind, only: LK, RKG => RK3
        real(RKG)       , intent(in), contiguous    :: array1(:)
        real(RKG)       , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if RK2_ENABLED
    pure module function islge_D1_D1_RK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_RK2
#endif
        use pm_kind, only: LK, RKG => RK2
        real(RKG)       , intent(in), contiguous    :: array1(:)
        real(RKG)       , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if RK1_ENABLED
    pure module function islge_D1_D1_RK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islge_D1_D1_RK1
#endif
        use pm_kind, only: LK, RKG => RK1
        real(RKG)       , intent(in), contiguous    :: array1(:)
        real(RKG)       , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_arrayCompareLex_islgt
    !>  Generate and return the result of the lexicographic comparison of two
    !>  input objects of the same type, kind, and rank using the `>` operator.
    !>
    !>  \details
    !>  For more information about the lexicographic comparisons defined in this generic procedure, in particular,
    !>  in the context of `complex` and `logical` arrays, see the detailed description of [pm_arrayCompareLex](@ref pm_arrayCompareLex).
    !>
    !>  \param[in]  array1  :   The input `contiguous` **array of rank `1`** of either<br>
    !>                          <ul>
    !>                              <li>    type `character` of kind \SKALL or, <br>
    !>                              <li>    type `integer` of kind \IKALL or, <br>
    !>                              <li>    type `logical` of kind \LKALL or, <br>
    !>                              <li>    type `complex` of kind \CKALL or, <br>
    !>                              <li>    type `real` of kind \RKALL or, <br>
    !>                          </ul>
    !>                          or,
    !>                          <ul>
    !>                              <li>    a **scalar assumed-length** `character` of kind \SKALL, <br>
    !>                          </ul>
    !>                          whose contents will be lexicographically compared with the contents of `array2`.
    !>  \param[in]  array2  :   The input of the same type, kind, and rank as `array1`.
    !>
    !>  \return
    !>  `compares`          :   The output scalar of type `logical` of default kind \LK containing the
    !>                          result of the lexical comparison of the values of the two input objects.
    !>
    !>  \interface{lgt}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_arrayCompareLex, only: operator(.lgt.)
    !>      logical(LK) :: compares
    !>
    !>      compares = array1 .lgt. array2
    !>
    !>  \endcode
    !>
    !>  \remark
    !>  The operator `.lgt.` stands for **lexically greater than**.
    !>
    !>  \pure
    !>
    !>  \see
    !>  [operator(.llt.)](@ref pm_arrayCompareLex_isllt)<br>
    !>  [operator(.lle.)](@ref pm_arrayCompareLex_islle)<br>
    !>  [operator(.lge.)](@ref pm_arrayCompareLex_islge)<br>
    !>  [operator(.lgt.)](@ref pm_arrayCompareLex_islgt)<br>
    !>
    !>  \example{lgt}
    !>  \include{lineno} example/pm_arrayCompareLex/lgt/main.F90
    !>  \compilef{lgt}
    !>  \output{lgt}
    !>  \include{lineno} example/pm_arrayCompareLex/lgt/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayCompareLex](@ref test_pm_arrayCompareLex)
    !>
    !>  \todo
    !>  \pvlow The functionality of this generic interface should be extended to input container arguments, also to arrays of higher ranks.
    !>
    !>  \final{lgt}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(.lgt.)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function islgt_D0_D0_SK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D0_D0_SK5
#endif
        use pm_kind, only: LK, SKG => SK5
        character(*,SKG), intent(in)                :: array1
        character(*,SKG), intent(in)                :: array2
        logical(LK)                                 :: compares
    end function
#endif

#if SK4_ENABLED
    pure module function islgt_D0_D0_SK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D0_D0_SK4
#endif
        use pm_kind, only: LK, SKG => SK4
        character(*,SKG), intent(in)                :: array1
        character(*,SKG), intent(in)                :: array2
        logical(LK)                                 :: compares
    end function
#endif

#if SK3_ENABLED
    pure module function islgt_D0_D0_SK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D0_D0_SK3
#endif
        use pm_kind, only: LK, SKG => SK3
        character(*,SKG), intent(in)                :: array1
        character(*,SKG), intent(in)                :: array2
        logical(LK)                                 :: compares
    end function
#endif

#if SK2_ENABLED
    pure module function islgt_D0_D0_SK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D0_D0_SK2
#endif
        use pm_kind, only: LK, SKG => SK2
        character(*,SKG), intent(in)                :: array1
        character(*,SKG), intent(in)                :: array2
        logical(LK)                                 :: compares
    end function
#endif

#if SK1_ENABLED
    pure module function islgt_D0_D0_SK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D0_D0_SK1
#endif
        use pm_kind, only: LK, SKG => SK1
        character(*,SKG), intent(in)                :: array1
        character(*,SKG), intent(in)                :: array2
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function islgt_D1_D1_SK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_SK5
#endif
        use pm_kind, only: LK, SKG => SK5
        character(*,SKG), intent(in), contiguous    :: array1(:)
        character(*,SKG), intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if SK4_ENABLED
    pure module function islgt_D1_D1_SK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_SK4
#endif
        use pm_kind, only: LK, SKG => SK4
        character(*,SKG), intent(in), contiguous    :: array1(:)
        character(*,SKG), intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if SK3_ENABLED
    pure module function islgt_D1_D1_SK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_SK3
#endif
        use pm_kind, only: LK, SKG => SK3
        character(*,SKG), intent(in), contiguous    :: array1(:)
        character(*,SKG), intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if SK2_ENABLED
    pure module function islgt_D1_D1_SK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_SK2
#endif
        use pm_kind, only: LK, SKG => SK2
        character(*,SKG), intent(in), contiguous    :: array1(:)
        character(*,SKG), intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if SK1_ENABLED
    pure module function islgt_D1_D1_SK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_SK1
#endif
        use pm_kind, only: LK, SKG => SK1
        character(*,SKG), intent(in), contiguous    :: array1(:)
        character(*,SKG), intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module function islgt_D1_D1_IK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_IK5
#endif
        use pm_kind, only: LK, IKG => IK5
        integer(IKG)    , intent(in), contiguous    :: array1(:)
        integer(IKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if IK4_ENABLED
    pure module function islgt_D1_D1_IK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_IK4
#endif
        use pm_kind, only: LK, IKG => IK4
        integer(IKG)    , intent(in), contiguous    :: array1(:)
        integer(IKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if IK3_ENABLED
    pure module function islgt_D1_D1_IK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_IK3
#endif
        use pm_kind, only: LK, IKG => IK3
        integer(IKG)    , intent(in), contiguous    :: array1(:)
        integer(IKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if IK2_ENABLED
    pure module function islgt_D1_D1_IK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_IK2
#endif
        use pm_kind, only: LK, IKG => IK2
        integer(IKG)    , intent(in), contiguous    :: array1(:)
        integer(IKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if IK1_ENABLED
    pure module function islgt_D1_D1_IK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_IK1
#endif
        use pm_kind, only: LK, IKG => IK1
        integer(IKG)    , intent(in), contiguous    :: array1(:)
        integer(IKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module function islgt_D1_D1_LK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_LK5
#endif
        use pm_kind, only: LK, LKG => LK5
        logical(LKG)    , intent(in), contiguous    :: array1(:)
        logical(LKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if LK4_ENABLED
    pure module function islgt_D1_D1_LK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_LK4
#endif
        use pm_kind, only: LK, LKG => LK4
        logical(LKG)    , intent(in), contiguous    :: array1(:)
        logical(LKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if LK3_ENABLED
    pure module function islgt_D1_D1_LK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_LK3
#endif
        use pm_kind, only: LK, LKG => LK3
        logical(LKG)    , intent(in), contiguous    :: array1(:)
        logical(LKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if LK2_ENABLED
    pure module function islgt_D1_D1_LK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_LK2
#endif
        use pm_kind, only: LK, LKG => LK2
        logical(LKG)    , intent(in), contiguous    :: array1(:)
        logical(LKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if LK1_ENABLED
    pure module function islgt_D1_D1_LK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_LK1
#endif
        use pm_kind, only: LK, LKG => LK1
        logical(LKG)    , intent(in), contiguous    :: array1(:)
        logical(LKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function islgt_D1_D1_CK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_CK5
#endif
        use pm_kind, only: LK, CKG => CK5
        complex(CKG)    , intent(in), contiguous    :: array1(:)
        complex(CKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if CK4_ENABLED
    pure module function islgt_D1_D1_CK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_CK4
#endif
        use pm_kind, only: LK, CKG => CK4
        complex(CKG)    , intent(in), contiguous    :: array1(:)
        complex(CKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if CK3_ENABLED
    pure module function islgt_D1_D1_CK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_CK3
#endif
        use pm_kind, only: LK, CKG => CK3
        complex(CKG)    , intent(in), contiguous    :: array1(:)
        complex(CKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if CK2_ENABLED
    pure module function islgt_D1_D1_CK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_CK2
#endif
        use pm_kind, only: LK, CKG => CK2
        complex(CKG)    , intent(in), contiguous    :: array1(:)
        complex(CKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if CK1_ENABLED
    pure module function islgt_D1_D1_CK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_CK1
#endif
        use pm_kind, only: LK, CKG => CK1
        complex(CKG)    , intent(in), contiguous    :: array1(:)
        complex(CKG)    , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module function islgt_D1_D1_RK5(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_RK5
#endif
        use pm_kind, only: LK, RKG => RK5
        real(RKG)       , intent(in), contiguous    :: array1(:)
        real(RKG)       , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if RK4_ENABLED
    pure module function islgt_D1_D1_RK4(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_RK4
#endif
        use pm_kind, only: LK, RKG => RK4
        real(RKG)       , intent(in), contiguous    :: array1(:)
        real(RKG)       , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if RK3_ENABLED
    pure module function islgt_D1_D1_RK3(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_RK3
#endif
        use pm_kind, only: LK, RKG => RK3
        real(RKG)       , intent(in), contiguous    :: array1(:)
        real(RKG)       , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if RK2_ENABLED
    pure module function islgt_D1_D1_RK2(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_RK2
#endif
        use pm_kind, only: LK, RKG => RK2
        real(RKG)       , intent(in), contiguous    :: array1(:)
        real(RKG)       , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

#if RK1_ENABLED
    pure module function islgt_D1_D1_RK1(array1, array2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D1_D1_RK1
#endif
        use pm_kind, only: LK, RKG => RK1
        real(RKG)       , intent(in), contiguous    :: array1(:)
        real(RKG)       , intent(in), contiguous    :: array2(:)
        logical(LK)                                 :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    \brief
!    \anchor pm_arrayCompareLex_islgt
!    Generate and return the result of the lexicographic comparison of two
!    input objects of the same type, kind, and rank using the `>` operator.
!
!    \details
!    For more information about the lexicographic comparisons defined in this generic procedure, in particular,
!    in the context of `complex` and `logical` arrays, see the detailed description of [pm_arrayCompareLex](@ref pm_arrayCompareLex).
!
!    \param[in]  array1  :   The input `contiguous` **array of rank `1`** of either<br>
!                            <ul>
!                                <li>    type `character` of kind \SKALL or, <br>
!                                <li>    type `integer` of kind \IKALL or, <br>
!                                <li>    type `logical` of kind \LKALL or, <br>
!                                <li>    type `complex` of kind \CKALL or, <br>
!                                <li>    type `real` of kind \RKALL or, <br>
!                            </ul>
!                            or,
!                            <ul>
!                                <li>    a **scalar assumed-length** `character` of kind \SKALL, <br>
!                            </ul>
!                            whose contents will be lexicographically compared with the contents of `array2`.
!    \param[in]  array2  :   The input of the same type, kind, and rank as `array1`.
!
!    \return
!    `compares`          :   The output scalar of type `logical` of default kind \LK containing the
!                            result of the lexical comparison of the values of the two input objects.
!
!    \interface{lgt}
!    \code{.F90}
!
!        use pm_kind, only: LK
!        use pm_arrayCompareLex, only: operator(.lgt.)
!        logical(LK) :: compares
!
!        compares = array1 .lgt. array2
!
!    \endcode
!
!    \remark
!    The generic interface name `islct` stands for **lexically greater than**.
!
!    \pure
!
!    \see
!    [operator(.llt.)](@ref pm_arrayCompareLex_isllt)<br>
!    [operator(.lle.)](@ref pm_arrayCompareLex_islle)<br>
!    [operator(.lge.)](@ref pm_arrayCompareLex_islge)<br>
!    [operator(.lgt.)](@ref pm_arrayCompareLex_islgt)<br>
!
!    \example{lgt}
!    \include{lineno} example/pm_arrayCompareLex/lgt/main.F90
!    \compile
!    \output
!    \include{lineno} example/pm_arrayCompareLex/lgt/main.out.F90
!
!    \test
!    [test_pm_arrayCompareLex](@ref test_pm_arrayCompareLex)
!
!    \todo
!    \pvlow The functionality of this generic interface should be extended to input container arguments, also to arrays of higher ranks.
!
!    \final
!
!    \author
!    Amir Shahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
!    interface islct
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    pure module function islct_UP(array1, array2, iscomparable) result(compares)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: islgt_D0_D0_SK
!#endif
!        use pm_kind, only: LK, SK
!        class(*)        , intent(in), contiguous    :: array1
!        class(*)        , intent(in), contiguous    :: array2
!        procedure(logical(LK))                      :: istrue
!        logical(LK)                                 :: compares
!    end function
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayCompareLex ! LCOV_EXCL_LINE