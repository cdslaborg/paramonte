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
!>  This module contains procedures and generic interfaces for checking if either of the corresponding
!>  real and imaginary components of two complex numbers satisfy a relational operator.
!>
!>  \details
!>  The primary purpose of the procedures in this module is to provide a convenient set of relational operators for generic programming.<br>
!>  Such use cases frequently occur in various library testing scenarios.<br>
!>
!>  \test
!>  [test_pm_complexCompareAny](@ref test_pm_complexCompareAny)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_complexCompareAny

    use pm_kind, only: LK, SK
    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_complexCompareAny"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_complexCompareAny_isanyless
    !>  Generate and return `.true.` if either components of the input `complex` argument `val1`
    !>  are less than the corresponding components of the input  `complex` argument `val2`.
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
    !>      use pm_complexCompareAny, only: operator(<)
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
    !>  [operator(<)](@ref pm_complexCompareAny_isanyless)<br>
    !>  [operator(<=)](@ref pm_complexCompareAny_isanyleq)<br>
    !>  [operator(.anyeq.)](@ref pm_complexCompareAny_isanyneq)<br>
    !>  [operator(.anyneq.)](@ref pm_complexCompareAny_isanyneq)<br>
    !>  [operator(>=)](@ref pm_complexCompareAny_isanymeq)<br>
    !>  [operator(>)](@ref pm_complexCompareAny_isanymore)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_complexCompareAny/isanyless/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_complexCompareAny/isanyless/main.out.F90
    !>
    !>  \test
    !>  [test_pm_complexCompareAny](@ref test_pm_complexCompareAny)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(<)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function isanyless_CK5(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanyless_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK4_ENABLED
    pure elemental module function isanyless_CK4(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanyless_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK3_ENABLED
    pure elemental module function isanyless_CK3(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanyless_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK2_ENABLED
    pure elemental module function isanyless_CK2(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanyless_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK1_ENABLED
    pure elemental module function isanyless_CK1(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanyless_CK1
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
    !>  \anchor pm_complexCompareAny_isanyleq
    !>  Generate and return `.true.` if either components of the input `complex` argument `val1` is
    !>  less than or equal to the corresponding components of the input `complex` argument `val2`.
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
    !>      use pm_complexCompareAny, only: operator(<=)
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
    !>  [operator(<)](@ref pm_complexCompareAny_isanyless)<br>
    !>  [operator(<=)](@ref pm_complexCompareAny_isanyleq)<br>
    !>  [operator(.anyeq.)](@ref pm_complexCompareAny_isanyneq)<br>
    !>  [operator(.anyneq.)](@ref pm_complexCompareAny_isanyneq)<br>
    !>  [operator(>=)](@ref pm_complexCompareAny_isanymeq)<br>
    !>  [operator(>)](@ref pm_complexCompareAny_isanymore)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_complexCompareAny/isanyleq/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_complexCompareAny/isanyleq/main.out.F90
    !>
    !>  \test
    !>  [test_pm_complexCompareAny](@ref test_pm_complexCompareAny)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(<=)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function isanyleq_CK5(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanyleq_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK4_ENABLED
    pure elemental module function isanyleq_CK4(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanyleq_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK3_ENABLED
    pure elemental module function isanyleq_CK3(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanyleq_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK2_ENABLED
    pure elemental module function isanyleq_CK2(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanyleq_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK1_ENABLED
    pure elemental module function isanyleq_CK1(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanyleq_CK1
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
    !>  \anchor pm_complexCompareAny_isanyeq
    !>  Generate and return `.true.` if either components of the input `complex` argument `val1`
    !>  is equal to the corresponding components of the input `complex` argument `val2`.
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
    !>      use pm_complexCompareAny, only: operator(==)
    !>      use pm_kind, only: LK
    !>      complex(LK) :: compares
    !>
    !>      compares = val1 == val2
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [operator(<)](@ref pm_complexCompareAny_isanyless)<br>
    !>  [operator(<=)](@ref pm_complexCompareAny_isanyleq)<br>
    !>  [operator(.anyeq.)](@ref pm_complexCompareAny_isanyneq)<br>
    !>  [operator(.anyneq.)](@ref pm_complexCompareAny_isanyneq)<br>
    !>  [operator(>=)](@ref pm_complexCompareAny_isanymeq)<br>
    !>  [operator(>)](@ref pm_complexCompareAny_isanymore)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_complexCompareAny/isanyeq/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_complexCompareAny/isanyeq/main.out.F90
    !>
    !>  \test
    !>  [test_pm_complexCompareAny](@ref test_pm_complexCompareAny)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(.anyeq.)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function isanyeq_CK5(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanyeq_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK4_ENABLED
    pure elemental module function isanyeq_CK4(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanyeq_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK3_ENABLED
    pure elemental module function isanyeq_CK3(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanyeq_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK2_ENABLED
    pure elemental module function isanyeq_CK2(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanyeq_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK1_ENABLED
    pure elemental module function isanyeq_CK1(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanyeq_CK1
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
    !>  \anchor pm_complexCompareAny_isanyneq
    !>  Generate and return `.true.` if either components of the input `complex` argument `val1`
    !>  are not equal to the corresponding components of the input `complex` argument `val2`.
    !>
    !>  \param[in]  val1  :   The input scalar, or array of the same rank and shape as `val2`, of type `complex` of kind \CKALL.
    !>  \param[in]  val2  :   The input scalar, or array of the same rank and shape as `val1`, of the same type and kind as `val1`.
    !>
    !>  \return
    !>  `compares`          :   The output object of type `logical` of default kind \LK, and
    !>                          the higher rank of the two input arguments `val1` and `val2`.
    !>
    !>  \remark
    !>  The intrinsic operator `/=` defined for complex numbers conflicts with the inequality definition in this interface.
    !>  As such, the interface is named `.anyneq.` to avoid the name clash.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_complexCompareAny, only: operator(.anyneq.)
    !>      use pm_kind, only: LK
    !>      complex(LK) :: compares
    !>
    !>      compares = val1 .anyneq. val2
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [operator(<)](@ref pm_complexCompareAny_isanyless)<br>
    !>  [operator(<=)](@ref pm_complexCompareAny_isanyleq)<br>
    !>  [operator(.anyeq.)](@ref pm_complexCompareAny_isanyneq)<br>
    !>  [operator(.anyneq.)](@ref pm_complexCompareAny_isanyneq)<br>
    !>  [operator(>=)](@ref pm_complexCompareAny_isanymeq)<br>
    !>  [operator(>)](@ref pm_complexCompareAny_isanymore)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_complexCompareAny/isanyneq/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_complexCompareAny/isanyneq/main.out.F90
    !>
    !>  \test
    !>  [test_pm_complexCompareAny](@ref test_pm_complexCompareAny)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(.anyneq.)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function isanyneq_CK5(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanyneq_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK4_ENABLED
    pure elemental module function isanyneq_CK4(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanyneq_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK3_ENABLED
    pure elemental module function isanyneq_CK3(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanyneq_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK2_ENABLED
    pure elemental module function isanyneq_CK2(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanyneq_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK1_ENABLED
    pure elemental module function isanyneq_CK1(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanyneq_CK1
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
    !>  \anchor pm_complexCompareAny_isanymeq
    !>  Generate and return `.true.` if either components of the input `complex` argument `val1`
    !>  is more than or equal to the corresponding components of the input `complex`   argument `val2`.
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
    !>      use pm_complexCompareAny, only: operator(>=)
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
    !>  [operator(<)](@ref pm_complexCompareAny_isanyless)<br>
    !>  [operator(<=)](@ref pm_complexCompareAny_isanyleq)<br>
    !>  [operator(.anyeq.)](@ref pm_complexCompareAny_isanyneq)<br>
    !>  [operator(.anyneq.)](@ref pm_complexCompareAny_isanyneq)<br>
    !>  [operator(>=)](@ref pm_complexCompareAny_isanymeq)<br>
    !>  [operator(>)](@ref pm_complexCompareAny_isanymore)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_complexCompareAny/isanymeq/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_complexCompareAny/isanymeq/main.out.F90
    !>
    !>  \test
    !>  [test_pm_complexCompareAny](@ref test_pm_complexCompareAny)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(>=)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function isanymeq_CK5(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanymeq_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK4_ENABLED
    pure elemental module function isanymeq_CK4(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanymeq_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK3_ENABLED
    pure elemental module function isanymeq_CK3(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanymeq_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK2_ENABLED
    pure elemental module function isanymeq_CK2(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanymeq_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK1_ENABLED
    pure elemental module function isanymeq_CK1(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanymeq_CK1
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
    !>  \anchor pm_complexCompareAny_isanymore
    !>  Generate and return `.true.` if either components of the input `complex` argument `val1`
    !>  is more than the corresponding components of the input `complex` argument `val2`.
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
    !>      use pm_complexCompareAny, only: operator(>)
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
    !>  [operator(<)](@ref pm_complexCompareAny_isanyless)<br>
    !>  [operator(<=)](@ref pm_complexCompareAny_isanyleq)<br>
    !>  [operator(.anyeq.)](@ref pm_complexCompareAny_isanyneq)<br>
    !>  [operator(.anyneq.)](@ref pm_complexCompareAny_isanyneq)<br>
    !>  [operator(>=)](@ref pm_complexCompareAny_isanymeq)<br>
    !>  [operator(>)](@ref pm_complexCompareAny_isanymore)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_complexCompareAny/isanymore/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_complexCompareAny/isanymore/main.out.F90
    !>
    !>  \test
    !>  [test_pm_complexCompareAny](@ref test_pm_complexCompareAny)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(>)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function isanymore_CK5(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanymore_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK4_ENABLED
    pure elemental module function isanymore_CK4(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanymore_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK3_ENABLED
    pure elemental module function isanymore_CK3(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanymore_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK2_ENABLED
    pure elemental module function isanymore_CK2(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanymore_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK1_ENABLED
    pure elemental module function isanymore_CK1(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isanymore_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_complexCompareAny ! LCOV_EXCL_LINE