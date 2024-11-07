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
!>  This module contains procedures and generic interfaces for checking if a `complex` number is
!>  lexicographically comparable to another `complex` number of the same kind.<br>
!>
!>  The two components of complex numbers are compared lexically.<br>
!>  For example, `(1.,2.) < (1.,1.)` evaluates to `.false.` while `(1.,2.) < (2.,0.)` evaluates to `.true.`.<br>
!>
!>  \details
!>  A **lexicographical comparison** is the kind of comparison generally used to sort words alphabetically in dictionaries.<br>
!>  It involves comparing sequentially the elements that have the same position in both ranges against each other until one element is not equivalent to the other.<br>
!>  The result of comparing these first non-matching elements is the result of the lexicographical comparison.<br>
!>
!>  \remark
!>  An example of a lexicographical comparison is the intrinsic comparison of two scalar characters in the Fortran programming language.<br>
!>
!>  \remark
!>  The primary purpose of the procedures in this module is to provide a convenient set of relational operators for generic programming.<br>
!>  Such use cases frequently occur in various library testing scenarios.<br>
!>
!>  \test
!>  [test_pm_complexCompareLex](@ref test_pm_complexCompareLex)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_complexCompareLex

    use pm_kind, only: SK, IK, LK
    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_complexCompareLex"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_complexCompareLex_islexless
    !>  Generate and return `.true.` if the real and imaginary components of the input `complex` argument `val1`
    !>  are lexicographically less than the corresponding components of the input  `complex` argument `val2`.
    !>
    !>  \param[in]  val1  :   The input scalar, or array of the same rank and shape as `val2`, of type `complex` of kind \CKALL.
    !>  \param[in]  val2  :   The input scalar, or array of the same rank and shape as `val1`, of the same type and kind as `val1`.
    !>
    !>  \return
    !>  `compares`          :   The output object of type `logical` of default kind \LK, and the higher rank of the two input arguments `val1` and `val2`.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_complexCompareLex, only: operator(<)
    !>      use pm_kind, only: LK
    !>      complex(LK) :: compares
    !>
    !>      compares = val1 < val2
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [operator(<)](@ref pm_complexCompareLex_islexless)<br>
    !>  [operator(<=)](@ref pm_complexCompareLex_islexleq)<br>
    !>  [operator(>=)](@ref pm_complexCompareLex_islexmeq)<br>
    !>  [operator(>)](@ref pm_complexCompareLex_islexmore)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_complexCompareLex/islexless/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_complexCompareLex/islexless/main.out.F90
    !>
    !>  \test
    !>  [test_pm_complexCompareLex](@ref test_pm_complexCompareLex)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(<)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function islexless_CK5(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islexless_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK4_ENABLED
    pure elemental module function islexless_CK4(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islexless_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK3_ENABLED
    pure elemental module function islexless_CK3(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islexless_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK2_ENABLED
    pure elemental module function islexless_CK2(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islexless_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK1_ENABLED
    pure elemental module function islexless_CK1(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islexless_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_complexCompareLex_islexleq
    !>  Generate and return `.true.` if the real and imaginary components of the input `complex` argument `val1`
    !>  are lexicographically less than or equal to the corresponding components of the input  `complex` argument `val2`.
    !>
    !>  \param[in]  val1  :   The input scalar, or array of the same rank and shape as `val2`, of type `complex` of kind \CKALL.
    !>  \param[in]  val2  :   The input scalar, or array of the same rank and shape as `val1`, of the same type and kind as `val1`.
    !>
    !>  \return
    !>  `compares`          :   The output object of type `logical` of default kind \LK, and
    !>                          the higher rank of the two input arguments `val1` and `val2`.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_complexCompareLex, only: operator(<=)
    !>      use pm_kind, only: LK
    !>      complex(LK) :: compares
    !>
    !>      compares = val1 <= val2
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [operator(<)](@ref pm_complexCompareLex_islexless)<br>
    !>  [operator(<=)](@ref pm_complexCompareLex_islexleq)<br>
    !>  [operator(>=)](@ref pm_complexCompareLex_islexmeq)<br>
    !>  [operator(>)](@ref pm_complexCompareLex_islexmore)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_complexCompareLex/islexleq/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_complexCompareLex/islexleq/main.out.F90
    !>
    !>  \test
    !>  [test_pm_complexCompareLex](@ref test_pm_complexCompareLex)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(<=)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function islexleq_CK5(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islexleq_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK4_ENABLED
    pure elemental module function islexleq_CK4(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islexleq_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK3_ENABLED
    pure elemental module function islexleq_CK3(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islexleq_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK2_ENABLED
    pure elemental module function islexleq_CK2(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islexleq_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK1_ENABLED
    pure elemental module function islexleq_CK1(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islexleq_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !>  \brief
!    !>  \anchor pm_complexCompareLex_isanyeq
!    !>  Generate and return `.true.` if either components of the input `complex` argument `val1`
!    !>  is equal to the corresponding components of the input `complex` argument `val2`.
!    !>
!    !>  \param[in]  val1  :   The input scalar, or array of the same rank and shape as `val2`, of type `complex` of kind \CKALL.
!    !>  \param[in]  val2  :   The input scalar, or array of the same rank and shape as `val1`, of the same type and kind as `val1`.
!    !>
!    !>  \return
!    !>  `compares`          :   The output object of type `logical` of default kind \LK, and
!    !>                          the higher rank of the two input arguments `val1` and `val2`.
!    !>
!    !>  \interface
!    !>  \code{.F90}
!    !>
!    !>      use pm_complexCompareLex, only: operator(==)
!    !>      use pm_kind, only: LK
!    !>      complex(LK) :: compares
!    !>
!    !>      compares = val1 == val2
!    !>
!    !>  \endcode
!    !>
!    !>  \pure
!    !>
!    !>  \elemental
!    !>
!    !>  \see
!    !>  [operator(<)](@ref pm_complexCompareLex_islexless)<br>
!    !>  [operator(<=)](@ref pm_complexCompareLex_islexleq)<br>
!    !>  [operator(.anyeq.)](@ref pm_complexCompareLex_islexneq)<br>
!    !>  [operator(.lexneq.)](@ref pm_complexCompareLex_islexneq)<br>
!    !>  [operator(>=)](@ref pm_complexCompareLex_islexmeq)<br>
!    !>  [operator(>)](@ref pm_complexCompareLex_islexmore)<br>
!    !>
!    !>  \example
!    !>  \include{lineno} example/pm_complexCompareLex/isanyeq/main.F90
!    !>  \compilef
!    !>  \output
!    !>  \include{lineno} example/pm_complexCompareLex/isanyeq/main.out.F90
!    !>
!    !>  \test
!    !>  [test_pm_complexCompareLex](@ref test_pm_complexCompareLex)
!    !>
!    !>  \final
!    !>
!    !>  \author
!    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
!    interface operator(.anyeq.)
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if CK3_ENABLED
!    pure elemental module function isanyeq_CK3(val1, val2) result(compares)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: isanyeq_CK3
!#endif
!        use pm_kind, only: CKG => CK3
!        complex(CKG), intent(in)                :: val1, val2
!        logical(LK)                             :: compares
!    end function
!#endif
!
!#if CK2_ENABLED
!    pure elemental module function isanyeq_CK2(val1, val2) result(compares)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: isanyeq_CK2
!#endif
!        use pm_kind, only: CKG => CK2
!        complex(CKG), intent(in)                :: val1, val2
!        logical(LK)                             :: compares
!    end function
!#endif
!
!#if CK1_ENABLED
!    pure elemental module function isanyeq_CK1(val1, val2) result(compares)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: isanyeq_CK1
!#endif
!        use pm_kind, only: CKG => CK1
!        complex(CKG), intent(in)                :: val1, val2
!        logical(LK)                             :: compares
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !>  \brief
!    !>  \anchor pm_complexCompareLex_islexneq
!    !>  Generate and return `.true.` if either components of the input `complex` argument `val1`
!    !>  are not equal to the corresponding components of the input `complex` argument `val2`.
!    !>
!    !>  \param[in]  val1  :   The input scalar, or array of the same rank and shape as `val2`, of type `complex` of kind \CKALL.
!    !>  \param[in]  val2  :   The input scalar, or array of the same rank and shape as `val1`, of the same type and kind as `val1`.
!    !>
!    !>  \return
!    !>  `compares`          :   The output object of type `logical` of default kind \LK, and
!    !>                          the higher rank of the two input arguments `val1` and `val2`.
!    !>
!    !>  \remark
!    !>  The intrinsic operator `/=` defined for complex numbers conflicts with the inequality definition in this interface.
!    !>  As such, the interface is named `.lexneq.` to avoid the name clash.
!    !>
!    !>  \interface
!    !>  \code{.F90}
!    !>
!    !>      use pm_complexCompareLex, only: operator(.lexneq.)
!    !>      use pm_kind, only: LK
!    !>      complex(LK) :: compares
!    !>
!    !>      compares = val1 .lexneq. val2
!    !>
!    !>  \endcode
!    !>
!    !>  \pure
!    !>
!    !>  \elemental
!    !>
!    !>  \see
!    !>  [operator(<)](@ref pm_complexCompareLex_islexless)<br>
!    !>  [operator(<=)](@ref pm_complexCompareLex_islexleq)<br>
!    !>  [operator(.anyeq.)](@ref pm_complexCompareLex_islexneq)<br>
!    !>  [operator(.lexneq.)](@ref pm_complexCompareLex_islexneq)<br>
!    !>  [operator(>=)](@ref pm_complexCompareLex_islexmeq)<br>
!    !>  [operator(>)](@ref pm_complexCompareLex_islexmore)<br>
!    !>
!    !>  \example
!    !>  \include{lineno} example/pm_complexCompareLex/islexneq/main.F90
!    !>  \compilef
!    !>  \output
!    !>  \include{lineno} example/pm_complexCompareLex/islexneq/main.out.F90
!    !>
!    !>  \test
!    !>  [test_pm_complexCompareLex](@ref test_pm_complexCompareLex)
!    !>
!    !>  \final
!    !>
!    !>  \author
!    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
!    interface operator(.lexneq.)
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if CK3_ENABLED
!    pure elemental module function islexneq_CK3(val1, val2) result(compares)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: islexneq_CK3
!#endif
!        use pm_kind, only: CKG => CK3
!        complex(CKG), intent(in)                :: val1, val2
!        logical(LK)                             :: compares
!    end function
!#endif
!
!#if CK2_ENABLED
!    pure elemental module function islexneq_CK2(val1, val2) result(compares)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: islexneq_CK2
!#endif
!        use pm_kind, only: CKG => CK2
!        complex(CKG), intent(in)                :: val1, val2
!        logical(LK)                             :: compares
!    end function
!#endif
!
!#if CK1_ENABLED
!    pure elemental module function islexneq_CK1(val1, val2) result(compares)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: islexneq_CK1
!#endif
!        use pm_kind, only: CKG => CK1
!        complex(CKG), intent(in)                :: val1, val2
!        logical(LK)                             :: compares
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_complexCompareLex_islexmeq
    !>  Generate and return `.true.` if the real and imaginary components of the input `complex` argument `val1`
    !>  are lexicographically more than or equal to the corresponding components of the input  `complex` argument `val2`.
    !>
    !>  \param[in]  val1  :   The input scalar, or array of the same rank and shape as `val2`, of type `complex` of kind \CKALL.
    !>  \param[in]  val2  :   The input scalar, or array of the same rank and shape as `val1`, of the same type and kind as `val1`.
    !>
    !>  \return
    !>  `compares`          :   The output object of type `logical` of default kind \LK, and
    !>                          the higher rank of the two input arguments `val1` and `val2`.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_complexCompareLex, only: operator(>=)
    !>      use pm_kind, only: LK
    !>      complex(LK) :: compares
    !>
    !>      compares = val1 >= val2
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [operator(<)](@ref pm_complexCompareLex_islexless)<br>
    !>  [operator(<=)](@ref pm_complexCompareLex_islexleq)<br>
    !>  [operator(>=)](@ref pm_complexCompareLex_islexmeq)<br>
    !>  [operator(>)](@ref pm_complexCompareLex_islexmore)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_complexCompareLex/islexmeq/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_complexCompareLex/islexmeq/main.out.F90
    !>
    !>  \test
    !>  [test_pm_complexCompareLex](@ref test_pm_complexCompareLex)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(>=)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function islexmeq_CK5(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islexmeq_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK4_ENABLED
    pure elemental module function islexmeq_CK4(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islexmeq_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK3_ENABLED
    pure elemental module function islexmeq_CK3(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islexmeq_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK2_ENABLED
    pure elemental module function islexmeq_CK2(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islexmeq_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK1_ENABLED
    pure elemental module function islexmeq_CK1(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islexmeq_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_complexCompareLex_islexmore
    !>  Generate and return `.true.` if the real and imaginary components of the input `complex` argument `val1`
    !>  are lexicographically more than the corresponding components of the input  `complex` argument `val2`.
    !>
    !>  \param[in]  val1  :   The input scalar, or array of the same rank and shape as `val2`, of type `complex` of kind \CKALL.
    !>  \param[in]  val2  :   The input scalar, or array of the same rank and shape as `val1`, of the same type and kind as `val1`.
    !>
    !>  \return
    !>  `compares`          :   The output object of type `logical` of default kind \LK, and
    !>                          the higher rank of the two input arguments `val1` and `val2`.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_complexCompareLex, only: operator(>)
    !>      use pm_kind, only: LK
    !>      complex(LK) :: compares
    !>
    !>      compares = val1 > val2
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [operator(<)](@ref pm_complexCompareLex_islexless)<br>
    !>  [operator(<=)](@ref pm_complexCompareLex_islexleq)<br>
    !>  [operator(>=)](@ref pm_complexCompareLex_islexmeq)<br>
    !>  [operator(>)](@ref pm_complexCompareLex_islexmore)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_complexCompareLex/islexmore/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_complexCompareLex/islexmore/main.out.F90
    !>
    !>  \test
    !>  [test_pm_complexCompareLex](@ref test_pm_complexCompareLex)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(>)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function islexmore_CK5(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islexmore_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK4_ENABLED
    pure elemental module function islexmore_CK4(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islexmore_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK3_ENABLED
    pure elemental module function islexmore_CK3(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islexmore_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK2_ENABLED
    pure elemental module function islexmore_CK2(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islexmore_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK1_ENABLED
    pure elemental module function islexmore_CK1(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: islexmore_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_complexCompareLex ! LCOV_EXCL_LINE