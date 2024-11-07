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
!>  This module contains procedures and generic interfaces for checking if both of the corresponding
!>  real and imaginary components of two complex numbers satisfy a relational operator.
!>
!>  \details
!>  The primary purpose of the procedures in this module is to provide a convenient set of relational operators for generic programming.<br>
!>  Such use cases frequently occur in various library testing scenarios.
!>
!>  \test
!>  [test_pm_complexCompareAll](@ref test_pm_complexCompareAll)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_complexCompareAll

    use pm_kind, only: LK, SK
    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_complexCompareAll"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_complexCompareAll_isallless
    !>  Generate and return `.true.` if both components of the input `complex` argument `val1`
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
    !>      use pm_complexCompareAll, only: operator(<)
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
    !>  [operator(<)](@ref pm_complexCompareAll_isallless)<br>
    !>  [operator(<=)](@ref pm_complexCompareAll_isallleq)<br>
    !>  [operator(.allneq.)](@ref pm_complexCompareAll_isallneq)<br>
    !>  [operator(>=)](@ref pm_complexCompareAll_isallmeq)<br>
    !>  [operator(>)](@ref pm_complexCompareAll_isallmore)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_complexCompareAll/isallless/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_complexCompareAll/isallless/main.out.F90
    !>
    !>  \test
    !>  [test_pm_complexCompareAll](@ref test_pm_complexCompareAll)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(<)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function isallless_CK5(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallless_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK4_ENABLED
    pure elemental module function isallless_CK4(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallless_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK3_ENABLED
    pure elemental module function isallless_CK3(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallless_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK2_ENABLED
    pure elemental module function isallless_CK2(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallless_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK1_ENABLED
    pure elemental module function isallless_CK1(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallless_CK1
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
    !>  \anchor pm_complexCompareAll_isallleq
    !>  Generate and return `.true.` if both components of the input `complex` argument `val1` is
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
    !>      use pm_complexCompareAll, only: operator(<=)
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
    !>  [operator(<)](@ref pm_complexCompareAll_isallless)<br>
    !>  [operator(<=)](@ref pm_complexCompareAll_isallleq)<br>
    !>  [operator(.allneq.)](@ref pm_complexCompareAll_isallneq)<br>
    !>  [operator(>=)](@ref pm_complexCompareAll_isallmeq)<br>
    !>  [operator(>)](@ref pm_complexCompareAll_isallmore)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_complexCompareAll/isallleq/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_complexCompareAll/isallleq/main.out.F90
    !>
    !>  \test
    !>  [test_pm_complexCompareAll](@ref test_pm_complexCompareAll)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(<=)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function isallleq_CK5(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallleq_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK4_ENABLED
    pure elemental module function isallleq_CK4(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallleq_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK3_ENABLED
    pure elemental module function isallleq_CK3(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallleq_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK2_ENABLED
    pure elemental module function isallleq_CK2(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallleq_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK1_ENABLED
    pure elemental module function isallleq_CK1(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallleq_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    !>  \brief
!    !>  \anchor pm_complexCompareAll_isalleq
!    !>  Generate and return `.true.` if both components of the input `complex` argument `val1`
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
!    !>      use pm_complexCompareAll, only: operator(==)
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
!    !>  [operator(<)](@ref pm_complexCompareAll_isallless)<br>
!    !>  [operator(<=)](@ref pm_complexCompareAll_isallleq)<br>
!    !>  [operator(.allneq.)](@ref pm_complexCompareAll_isallneq)<br>
!    !>  [operator(>=)](@ref pm_complexCompareAll_isallmeq)<br>
!    !>  [operator(>)](@ref pm_complexCompareAll_isallmore)<br>
!    !>
!    !>  \example
!    !>  \include{lineno} example/pm_complexCompareAll/isalleq/main.F90
!    !>  \compilef
!    !>  \output
!    !>  \include{lineno} example/pm_complexCompareAll/isalleq/main.out.F90
!    !>
!    !>  \test
!    !>  [test_pm_complexCompareAll](@ref test_pm_complexCompareAll)
!    !>
!    !>  \final
!    !>
!    !>  \author
!    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
!    interface operator(==)
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if CK3_ENABLED
!    pure elemental module function isalleq_CK3(val1, val2) result(compares)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: isalleq_CK3
!#endif
!        use pm_kind, only: CKG => CK3
!        complex(CKG), intent(in)                :: val1, val2
!        logical(LK)                             :: compares
!    end function
!#endif
!
!#if CK2_ENABLED
!    pure elemental module function isalleq_CK2(val1, val2) result(compares)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: isalleq_CK2
!#endif
!        use pm_kind, only: CKG => CK2
!        complex(CKG), intent(in)                :: val1, val2
!        logical(LK)                             :: compares
!    end function
!#endif
!
!#if CK1_ENABLED
!    pure elemental module function isalleq_CK1(val1, val2) result(compares)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: isalleq_CK1
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
    !>  \anchor pm_complexCompareAll_isallneq
    !>  Generate and return `.true.` if both components of the input `complex` argument `val1`
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
    !>  As such, the interface is named `.allneq.` to avoid the name clash.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_complexCompareAll, only: operator(.allneq.)
    !>      use pm_kind, only: LK
    !>      complex(LK) :: compares
    !>
    !>      compares = val1 .allneq. val2
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [operator(<)](@ref pm_complexCompareAll_isallless)<br>
    !>  [operator(<=)](@ref pm_complexCompareAll_isallleq)<br>
    !>  [operator(.allneq.)](@ref pm_complexCompareAll_isallneq)<br>
    !>  [operator(>=)](@ref pm_complexCompareAll_isallmeq)<br>
    !>  [operator(>)](@ref pm_complexCompareAll_isallmore)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_complexCompareAll/isallneq/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_complexCompareAll/isallneq/main.out.F90
    !>
    !>  \test
    !>  [test_pm_complexCompareAll](@ref test_pm_complexCompareAll)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(.allneq.)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function isallneq_CK5(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallneq_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK4_ENABLED
    pure elemental module function isallneq_CK4(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallneq_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK3_ENABLED
    pure elemental module function isallneq_CK3(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallneq_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK2_ENABLED
    pure elemental module function isallneq_CK2(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallneq_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK1_ENABLED
    pure elemental module function isallneq_CK1(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallneq_CK1
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
    !>  \anchor pm_complexCompareAll_isallmeq
    !>  Generate and return `.true.` if both components of the input `complex` argument `val1`
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
    !>      use pm_complexCompareAll, only: operator(>=)
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
    !>  [operator(<)](@ref pm_complexCompareAll_isallless)<br>
    !>  [operator(<=)](@ref pm_complexCompareAll_isallleq)<br>
    !>  [operator(.allneq.)](@ref pm_complexCompareAll_isallneq)<br>
    !>  [operator(>=)](@ref pm_complexCompareAll_isallmeq)<br>
    !>  [operator(>)](@ref pm_complexCompareAll_isallmore)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_complexCompareAll/isallmeq/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_complexCompareAll/isallmeq/main.out.F90
    !>
    !>  \test
    !>  [test_pm_complexCompareAll](@ref test_pm_complexCompareAll)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(>=)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function isallmeq_CK5(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallmeq_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK4_ENABLED
    pure elemental module function isallmeq_CK4(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallmeq_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK3_ENABLED
    pure elemental module function isallmeq_CK3(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallmeq_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK2_ENABLED
    pure elemental module function isallmeq_CK2(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallmeq_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK1_ENABLED
    pure elemental module function isallmeq_CK1(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallmeq_CK1
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
    !>  \anchor pm_complexCompareAll_isallmore
    !>  Generate and return `.true.` if both components of the input `complex` argument `val1`
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
    !>      use pm_complexCompareAll, only: operator(>)
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
    !>  [operator(<)](@ref pm_complexCompareAll_isallless)<br>
    !>  [operator(<=)](@ref pm_complexCompareAll_isallleq)<br>
    !>  [operator(.allneq.)](@ref pm_complexCompareAll_isallneq)<br>
    !>  [operator(>=)](@ref pm_complexCompareAll_isallmeq)<br>
    !>  [operator(>)](@ref pm_complexCompareAll_isallmore)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_complexCompareAll/isallmore/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_complexCompareAll/isallmore/main.out.F90
    !>
    !>  \test
    !>  [test_pm_complexCompareAll](@ref test_pm_complexCompareAll)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(>)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function isallmore_CK5(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallmore_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK4_ENABLED
    pure elemental module function isallmore_CK4(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallmore_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK3_ENABLED
    pure elemental module function isallmore_CK3(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallmore_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK2_ENABLED
    pure elemental module function isallmore_CK2(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallmore_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

#if CK1_ENABLED
    pure elemental module function isallmore_CK1(val1, val2) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isallmore_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)                :: val1, val2
        logical(LK)                             :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_complexCompareAll ! LCOV_EXCL_LINE