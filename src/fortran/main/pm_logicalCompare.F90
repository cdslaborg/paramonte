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
!>  This module contains procedures and generic interfaces for performing a variety of logical comparison operations
!>  using `logical` values as if `.true.` evaluates to `1` and `.false.` evaluates to `0`.
!>
!>  \test
!>  [test_pm_logicalCompare](@ref test_pm_logicalCompare)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_logicalCompare

    use pm_kind, only: SK, IK
    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_logicalCompare"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_logicalCompare_isless
    !>  Generate and return `.true.` if the input `logical` argument `lhs` is less than the input `logical` argument `rhs`.
    !>
    !>  \param[in]  lhs :   The input scalar or array of the same rank and shape as `rhs` of type `logical` of default kind \LK.
    !>  \param[in]  rhs :   The input scalar or array of the same rank and shape as `lhs` of type `logical` of default kind \LK.
    !>
    !>  \return
    !>  `compares`      :   The output object of the same type and kind as, and the higher rank of, the two input arguments `lhs` and `rhs`.
    !>
    !>  \interface{pm_logicalCompare_isless}
    !>  \code{.F90}
    !>
    !>      use pm_logicalCompare, only: operator(<)
    !>      use pm_kind, only: LK
    !>      logical(LK) :: compares
    !>
    !>      compares = lhs < rhs
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [operator(<=)](@ref pm_logicalCompare_isless)<br>
    !>
    !>  \example{pm_logicalCompare_isless}
    !>  \include{lineno} example/pm_logicalCompare/isless/main.F90
    !>  \compilef{pm_logicalCompare_isless}
    !>  \output{pm_logicalCompare_isless}
    !>  \include{lineno} example/pm_logicalCompare/isless/main.out.F90
    !>
    !>  \test
    !>  [test_pm_logicalCompare](@ref test_pm_logicalCompare)
    !>
    !>  \final{pm_logicalCompare_isless}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(<)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure elemental module function isless_LK5(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_LK5
#endif
        use pm_kind, only: LKG => LK5 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK4_ENABLED
    pure elemental module function isless_LK4(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_LK4
#endif
        use pm_kind, only: LKG => LK4 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK3_ENABLED
    pure elemental module function isless_LK3(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_LK3
#endif
        use pm_kind, only: LKG => LK3 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK2_ENABLED
    pure elemental module function isless_LK2(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_LK2
#endif
        use pm_kind, only: LKG => LK2 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK1_ENABLED
    pure elemental module function isless_LK1(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_LK1
#endif
        use pm_kind, only: LKG => LK1 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_logicalCompare_isleq
    !>  Generate and return `.true.` if the input `logical` argument `lhs` is less than or equal to the input `logical` argument `rhs`.
    !>
    !>  \param[in]  lhs :   The input scalar or array of the same rank and shape as `rhs` of type `logical` of default kind \LK.
    !>  \param[in]  rhs :   The input scalar or array of the same rank and shape as `lhs` of type `logical` of default kind \LK.
    !>
    !>  \return
    !>  `compares`      :   The output object of the same type and kind as, and the higher rank of, the two input arguments `lhs` and `rhs`.
    !>
    !>  \interface{pm_logicalCompare_isleq}
    !>  \code{.F90}
    !>
    !>      use pm_logicalCompare, only: operator(<=)
    !>      use pm_kind, only: LK
    !>      logical(LK) :: compares
    !>
    !>      compares = lhs <= rhs
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [operator(<)](@ref pm_logicalCompare_isless)<br>
    !>
    !>  \example{pm_logicalCompare_isleq}
    !>  \include{lineno} example/pm_logicalCompare/isleq/main.F90
    !>  \compilef{pm_logicalCompare_isleq}
    !>  \output{pm_logicalCompare_isleq}
    !>  \include{lineno} example/pm_logicalCompare/isleq/main.out.F90
    !>
    !>  \test
    !>  [test_pm_logicalCompare](@ref test_pm_logicalCompare)
    !>
    !>  \final{pm_logicalCompare_isleq}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(<=)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure elemental module function isleq_LK5(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_LK5
#endif
        use pm_kind, only: LKG => LK5 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK4_ENABLED
    pure elemental module function isleq_LK4(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_LK4
#endif
        use pm_kind, only: LKG => LK4 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK3_ENABLED
    pure elemental module function isleq_LK3(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_LK3
#endif
        use pm_kind, only: LKG => LK3 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK2_ENABLED
    pure elemental module function isleq_LK2(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_LK2
#endif
        use pm_kind, only: LKG => LK2 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK1_ENABLED
    pure elemental module function isleq_LK1(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_LK1
#endif
        use pm_kind, only: LKG => LK1 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_logicalCompare_iseq
    !>  Generate and return `.true.` if the input `logical` argument `lhs` is equal to the input `logical` argument `rhs`.
    !>
    !>  \param[in]  lhs :   The input scalar or array of the same rank and shape as `rhs` of type `logical` of default kind \LK.
    !>  \param[in]  rhs :   The input scalar or array of the same rank and shape as `lhs` of type `logical` of default kind \LK.
    !>
    !>  \return
    !>  `compares`      :   The output object of the same type and kind as, and the higher rank of, the two input arguments `lhs` and `rhs`.
    !>
    !>  \interface{pm_logicalCompare_iseq}
    !>  \code{.F90}
    !>
    !>      use pm_logicalCompare, only: operator(==)
    !>      use pm_kind, only: LK
    !>      logical(LK) :: compares
    !>
    !>      compares = lhs == rhs
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [operator(<)](@ref pm_logicalCompare_isless)<br>
    !>
    !>  \example{pm_logicalCompare_iseq}
    !>  \include{lineno} example/pm_logicalCompare/iseq/main.F90
    !>  \compilef{pm_logicalCompare_iseq}
    !>  \output{pm_logicalCompare_iseq}
    !>  \include{lineno} example/pm_logicalCompare/iseq/main.out.F90
    !>
    !>  \test
    !>  [test_pm_logicalCompare](@ref test_pm_logicalCompare)
    !>
    !>  \final{pm_logicalCompare_iseq}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(==)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure elemental module function iseq_LK5(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_LK5
#endif
        use pm_kind, only: LKG => LK5 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK4_ENABLED
    pure elemental module function iseq_LK4(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_LK4
#endif
        use pm_kind, only: LKG => LK4 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK3_ENABLED
    pure elemental module function iseq_LK3(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_LK3
#endif
        use pm_kind, only: LKG => LK3 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK2_ENABLED
    pure elemental module function iseq_LK2(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_LK2
#endif
        use pm_kind, only: LKG => LK2 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK1_ENABLED
    pure elemental module function iseq_LK1(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_LK1
#endif
        use pm_kind, only: LKG => LK1 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_logicalCompare_isneq
    !>  Generate and return `.true.` if the input `logical` argument `lhs` is not equal to the input `logical` argument `rhs`.
    !>
    !>  \param[in]  lhs :   The input scalar or array of the same rank and shape as `rhs` of type `logical` of default kind \LK.
    !>  \param[in]  rhs :   The input scalar or array of the same rank and shape as `lhs` of type `logical` of default kind \LK.
    !>
    !>  \return
    !>  `compares`      :   The output object of the same type and kind as, and the higher rank of, the two input arguments `lhs` and `rhs`.
    !>
    !>  \interface{pm_logicalCompare_isneq}
    !>  \code{.F90}
    !>
    !>      use pm_logicalCompare, only: operator(/=)
    !>      use pm_kind, only: LK
    !>      logical(LK) :: compares
    !>
    !>      compares = lhs /= rhs
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [operator(<)](@ref pm_logicalCompare_isless)<br>
    !>
    !>  \example{pm_logicalCompare_isneq}
    !>  \include{lineno} example/pm_logicalCompare/isneq/main.F90
    !>  \compilef{pm_logicalCompare_isneq}
    !>  \output{pm_logicalCompare_isneq}
    !>  \include{lineno} example/pm_logicalCompare/isneq/main.out.F90
    !>
    !>  \test
    !>  [test_pm_logicalCompare](@ref test_pm_logicalCompare)
    !>
    !>  \final{pm_logicalCompare_isneq}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(/=)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure elemental module function isneq_LK5(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_LK5
#endif
        use pm_kind, only: LKG => LK5 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK4_ENABLED
    pure elemental module function isneq_LK4(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_LK4
#endif
        use pm_kind, only: LKG => LK4 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK3_ENABLED
    pure elemental module function isneq_LK3(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_LK3
#endif
        use pm_kind, only: LKG => LK3 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK2_ENABLED
    pure elemental module function isneq_LK2(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_LK2
#endif
        use pm_kind, only: LKG => LK2 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK1_ENABLED
    pure elemental module function isneq_LK1(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_LK1
#endif
        use pm_kind, only: LKG => LK1 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_logicalCompare_ismeq
    !>  Generate and return `.true.` if the input `logical` argument `lhs` is more than or equal to the input `logical`   argument `rhs`.
    !>
    !>  \param[in]  lhs :   The input scalar or array of the same rank and shape as `rhs` of type `logical` of default kind \LK.
    !>  \param[in]  rhs :   The input scalar or array of the same rank and shape as `lhs` of type `logical` of default kind \LK.
    !>
    !>  \return
    !>  `compares`      :   The output object of the same type and kind as, and the higher rank of, the two input arguments `lhs` and `rhs`.
    !>
    !>  \interface{pm_logicalCompare_ismeq}
    !>  \code{.F90}
    !>
    !>      use pm_logicalCompare, only: operator(>=)
    !>      use pm_kind, only: LK
    !>      logical(LK) :: compares
    !>
    !>      compares = lhs >= rhs
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [operator(<)](@ref pm_logicalCompare_isless)<br>
    !>
    !>  \example{pm_logicalCompare_ismeq}
    !>  \include{lineno} example/pm_logicalCompare/ismeq/main.F90
    !>  \compilef{pm_logicalCompare_ismeq}
    !>  \output{pm_logicalCompare_ismeq}
    !>  \include{lineno} example/pm_logicalCompare/ismeq/main.out.F90
    !>
    !>  \test
    !>  [test_pm_logicalCompare](@ref test_pm_logicalCompare)
    !>
    !>  \final{pm_logicalCompare_ismeq}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(>=)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure elemental module function ismeq_LK5(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_LK5
#endif
        use pm_kind, only: LKG => LK5 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK4_ENABLED
    pure elemental module function ismeq_LK4(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_LK4
#endif
        use pm_kind, only: LKG => LK4 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK3_ENABLED
    pure elemental module function ismeq_LK3(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_LK3
#endif
        use pm_kind, only: LKG => LK3 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK2_ENABLED
    pure elemental module function ismeq_LK2(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_LK2
#endif
        use pm_kind, only: LKG => LK2 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK1_ENABLED
    pure elemental module function ismeq_LK1(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_LK1
#endif
        use pm_kind, only: LKG => LK1 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_logicalCompare_ismore
    !>  Generate and return `.true.` if the input `logical` argument `lhs` is more than the input `logical` argument `rhs`.
    !>
    !>  \param[in]  lhs :   The input scalar or array of the same rank and shape as `rhs` of type `logical` of default kind \LK.
    !>  \param[in]  rhs :   The input scalar or array of the same rank and shape as `lhs` of type `logical` of default kind \LK.
    !>
    !>  \return
    !>  `compares`      :   The output object of the same type and kind as, and the higher rank of, the two input arguments `lhs` and `rhs`.
    !>
    !>  \interface{pm_logicalCompare_ismore}
    !>  \code{.F90}
    !>
    !>      use pm_logicalCompare, only: operator(>)
    !>      use pm_kind, only: LK
    !>      logical(LK) :: compares
    !>
    !>      compares = lhs > rhs
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [operator(<)](@ref pm_logicalCompare_isless)<br>
    !>
    !>  \example{pm_logicalCompare_ismore}
    !>  \include{lineno} example/pm_logicalCompare/ismore/main.F90
    !>  \compilef{pm_logicalCompare_ismore}
    !>  \output{pm_logicalCompare_ismore}
    !>  \include{lineno} example/pm_logicalCompare/ismore/main.out.F90
    !>
    !>  \test
    !>  [test_pm_logicalCompare](@ref test_pm_logicalCompare)
    !>
    !>  \final{pm_logicalCompare_ismore}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(>)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure elemental module function ismore_LK5(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_LK5
#endif
        use pm_kind, only: LKG => LK5 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK4_ENABLED
    pure elemental module function ismore_LK4(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_LK4
#endif
        use pm_kind, only: LKG => LK4 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK3_ENABLED
    pure elemental module function ismore_LK3(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_LK3
#endif
        use pm_kind, only: LKG => LK3 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK2_ENABLED
    pure elemental module function ismore_LK2(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_LK2
#endif
        use pm_kind, only: LKG => LK2 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

#if LK1_ENABLED
    pure elemental module function ismore_LK1(lhs, rhs) result(compares)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_LK1
#endif
        use pm_kind, only: LKG => LK1 
        logical(LKG)            , intent(in)                :: lhs, rhs
        logical(LKG)                                        :: compares
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_logicalCompare ! LCOV_EXCL_LINE
