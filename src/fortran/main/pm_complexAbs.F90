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
!>  This module contains procedures and generic interfaces for performing element-wise comparison of the
!>  real and imaginary components of scalars and arrays of arbitrary ranks of various types.
!>
!>  \remark
!>  The primary purpose of the procedures in this module is to provide a convenient
!>  element-wise comparison of the real and imaginary components of complex numbers.<br>
!>  Such cases frequently occur in various library testing scenarios.<br>
!>  Otherwise, these procedures have minimal mathematical utilities.<br>
!>
!>  \test
!>  [test_pm_complexAbs](@ref test_pm_complexAbs)
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_complexAbs

    use pm_kind, only: LK, SK
    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_complexAbs"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the component-wise absolute value of (both real and imaginary parts of) the input `complex`.
    !>
    !>  \param[in]  val :   The input scalar or array of arbitrary rank of type `complex` of kind \CKALL whose elemental absolute value is to be returned.
    !>
    !>  \return
    !>  `absval`        :   The output scalar or array of the same rank and shape as `val` whose real and imaginary parts
    !>                      are the absolute values of the corresponding real and imaginary parts of the input `val`.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_complexAbs, only: abs
    !>
    !>      absval = abs(val) ! = (abs(val%re), abs(val%im))
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [pm_complexCompareAll](@ref pm_complexCompareAll)<br>
    !>  [pm_complexCompareAny](@ref pm_complexCompareAny)<br>
    !>  [pm_complexCompareLex](@ref pm_complexCompareLex)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_complexAbs/abs/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_complexAbs/abs/main.out.F90
    !>
    !>  \test
    !>  [test_pm_complexAbs](@ref test_pm_complexAbs)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface abs

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function abs_CK5(val) result(absval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: abs_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                :: val
        complex(CKG)                            :: absval
    end function
#endif

#if CK4_ENABLED
    pure elemental module function abs_CK4(val) result(absval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: abs_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                :: val
        complex(CKG)                            :: absval
    end function
#endif

#if CK3_ENABLED
    pure elemental module function abs_CK3(val) result(absval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: abs_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                :: val
        complex(CKG)                            :: absval
    end function
#endif

#if CK2_ENABLED
    pure elemental module function abs_CK2(val) result(absval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: abs_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                :: val
        complex(CKG)                            :: absval
    end function
#endif

#if CK1_ENABLED
    pure elemental module function abs_CK1(val) result(absval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: abs_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)                :: val
        complex(CKG)                            :: absval
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_complexAbs ! LCOV_EXCL_LINE