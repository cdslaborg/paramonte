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
!>  This module contains procedures and generic interfaces for computing the complex division robustly without potential overflow of computations.<br>
!>
!>  \note
!>  It appears that most compilers, including Intel and GNU Fortran compilers perform complex division robustly, by default.<br>
!>  As such, the explicit usage of the routines of this module over the default complex division of the Fortran standard may not be warranted.<br>
!>
!>  \test
!>  [test_pm_complexDiv](@ref test_pm_complexDiv)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_complexDiv

    use pm_kind, only: LK, SK
    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_complexDiv"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the complex division of the two input `complex` values **robustly**
    !>  without potential overflow caused by the computing the modulus of the divisor in the denominator of the result.<br>
    !>
    !>  \param[in]  dividend    :   The input scalar or (array of the same shape as other array-like input arguments) of type `complex` of kind \CKALL,
    !>                              representing the numerator of the division.
    !>  \param[in]  divisor     :   The input scalar or (array of the same shape as other array-like input arguments) of same type and kind as `dividend`,
    !>                              representing the denominator of the division.
    !>
    !>  \return
    !>  `quotient`              :   The input scalar or (array of the same shape as other array-like input arguments) of same type and kind as `dividend`,
    !>                              representing the complex division result.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_complexDiv, only: getDiv
    !>
    !>      quotient = getDiv(dividend, divisor)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `divisor /= (0._CKG, 0._CKG)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [pm_complexCompareAll](@ref pm_complexCompareAll)<br>
    !>  [pm_complexCompareAny](@ref pm_complexCompareAny)<br>
    !>  [pm_complexCompareLex](@ref pm_complexCompareLex)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_complexDiv/getDiv/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_complexDiv/getDiv/main.out.F90
    !>
    !>  \test
    !>  [test_pm_complexDiv](@ref test_pm_complexDiv)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getDiv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module function getDiv_CK5(dividend, divisor) result(quotient)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDiv_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                :: dividend, divisor
        complex(CKG)                            :: quotient
    end function
#endif

#if CK4_ENABLED
    PURE elemental module function getDiv_CK4(dividend, divisor) result(quotient)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDiv_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                :: dividend, divisor
        complex(CKG)                            :: quotient
    end function
#endif

#if CK3_ENABLED
    PURE elemental module function getDiv_CK3(dividend, divisor) result(quotient)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDiv_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                :: dividend, divisor
        complex(CKG)                            :: quotient
    end function
#endif

#if CK2_ENABLED
    PURE elemental module function getDiv_CK2(dividend, divisor) result(quotient)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDiv_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                :: dividend, divisor
        complex(CKG)                            :: quotient
    end function
#endif

#if CK1_ENABLED
    PURE elemental module function getDiv_CK1(dividend, divisor) result(quotient)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDiv_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)                :: dividend, divisor
        complex(CKG)                            :: quotient
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_complexDiv ! LCOV_EXCL_LINE