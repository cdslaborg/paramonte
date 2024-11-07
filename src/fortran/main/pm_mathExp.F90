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
!>  This module contains procedures and generic interfaces for computing the
!>  <b>previous/next integer exponent</b> for the given base that yields a number smaller/larger than the absolute input value.
!>
!>  \details
!>  For a given number \f$x\f$ and `base` \f$b\f$, the <b>next integer exponent</b> \f$p\f$ is defined as the **smallest integer** such that,<br>
!>  \f{equation}{
!>      b^{(p - 1)} < |x| \leq b^p ~~~ \forall x > 0 ~,
!>  \f}
!>
!>  \note
!>  The functionalities of this module are useful in computing the padded lengths of arrays whose
!>  auto-correlation, cross-correlation, or Fourier-Transform are to be computed using the Fast-Fourier Transform method.<br>
!>
!>  \see
!>  [pm_fftpack](@ref pm_fftpack)<br>
!>
!>  \test
!>  [test_pm_mathExp](@ref test_pm_mathExp)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 25, 2015, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_mathExp

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathExp"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` <b>if and only if</b> the input positive `integer` is a whole-number (integer) power of the specified integer `base`.
    !>
    !>  \details
    !>  The procedures of this generic interface rely on the binary
    !>  representation of numbers to determine if the number is a power of base two.<br>
    !>  For all other cases, an iterative division approach is used to determine the output.<br>
    !>
    !>  \param[in]  absx    :   The input scalar (or array of the same rank, shape, and size as other array-like arguments) of,<br>
    !>                          <ol>
    !>                              <li>    type `integer` of kind \IKALL,<br>
    !>                          </ol>
    !>                          representing the **absolute** value of the number for which the next integer exponent in specified `base` must be computed.<br>
    !>  \param[in]  base    :   The input positive scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `absx`,
    !>                          representing the base of the exponentiation.<br>
    !>                          (**optional**, default = `2`)
    !>
    !>  \return
    !>  `powisint`          :   The output scalar (or array of the same rank, shape, and size as other array-like arguments), of the type `logical` of default kind \LK
    !>                          that is `.true.` if and only if the input `absx` is an integer power of the input `absx`.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_mathExp, only: isIntPow
    !>      logical(LK) :: powisint
    !>
    !>      powisint = isIntPow(absx)
    !>      powisint = isIntPow(absx, base)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < absx` must hold.<br>
    !>  The condition `1 <= base` must hold.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [isIntPow](@ref pm_mathExp::isIntPow)<br>
    !>  [getExpNext](@ref pm_mathExp::getExpNext)<br>
    !>  [getExpPrev](@ref pm_mathExp::getExpPrev)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_mathExp/isIntPow/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_mathExp/isIntPow/main.out.F90
    !>
    !>  \test
    !>  [test_pm_mathExp](@ref test_pm_mathExp)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, April 25, 2015, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface isIntPow

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE elemental module function isIntPowArb_IK5(absx, base) result(powisint)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isIntPowArb_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG), intent(in)            :: absx
        integer(IKG), intent(in)            :: base
        logical(LK)                         :: powisint
    end function
#endif

#if IK4_ENABLED
    PURE elemental module function isIntPowArb_IK4(absx, base) result(powisint)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isIntPowArb_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG), intent(in)            :: absx
        integer(IKG), intent(in)            :: base
        logical(LK)                         :: powisint
    end function
#endif

#if IK3_ENABLED
    PURE elemental module function isIntPowArb_IK3(absx, base) result(powisint)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isIntPowArb_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG), intent(in)            :: absx
        integer(IKG), intent(in)            :: base
        logical(LK)                         :: powisint
    end function
#endif

#if IK2_ENABLED
    PURE elemental module function isIntPowArb_IK2(absx, base) result(powisint)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isIntPowArb_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG), intent(in)            :: absx
        integer(IKG), intent(in)            :: base
        logical(LK)                         :: powisint
    end function
#endif

#if IK1_ENABLED
    PURE elemental module function isIntPowArb_IK1(absx, base) result(powisint)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isIntPowArb_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG), intent(in)            :: absx
        integer(IKG), intent(in)            :: base
        logical(LK)                         :: powisint
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE elemental module function isIntPowDef_IK5(absx) result(powisint)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isIntPowDef_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG), intent(in)            :: absx
        logical(LK)                         :: powisint
    end function
#endif

#if IK4_ENABLED
    PURE elemental module function isIntPowDef_IK4(absx) result(powisint)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isIntPowDef_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG), intent(in)            :: absx
        logical(LK)                         :: powisint
    end function
#endif

#if IK3_ENABLED
    PURE elemental module function isIntPowDef_IK3(absx) result(powisint)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isIntPowDef_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG), intent(in)            :: absx
        logical(LK)                         :: powisint
    end function
#endif

#if IK2_ENABLED
    PURE elemental module function isIntPowDef_IK2(absx) result(powisint)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isIntPowDef_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG), intent(in)            :: absx
        logical(LK)                         :: powisint
    end function
#endif

#if IK1_ENABLED
    PURE elemental module function isIntPowDef_IK1(absx) result(powisint)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isIntPowDef_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG), intent(in)            :: absx
        logical(LK)                         :: powisint
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate the <b>next integer exponent</b> `expNext` for the specified input `base` and absolute value `absx = abs(x)`
    !>  such that the condition `absx <= base**expNext` holds.
    !>
    !>  \details
    !>  See the documentation of [pm_mathExp](@ref pm_mathExp) for more information.<br>
    !>
    !>  \param[in]  absx    :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of either,<br>
    !>                          <ol>
    !>                              <li>    type `integer` of kind \IKALL, or<br>
    !>                              <li>    type `real` of kind \RKALL,<br>
    !>                          </ol>
    !>                          representing the **absolute** value of the number for which the next integer exponent in specified `base` must be computed.<br>
    !>  \param[in]  base    :   The input positive scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `absx`,
    !>                          representing the base of the exponentiation.<br>
    !>                          (**optional**, default = `2`)
    !>
    !>  \return
    !>  `expNext`           :   The output scalar (or array of the same rank, shape, and size as other array-like arguments) of
    !>                          <ol>
    !>                              <li>    type `integer` of default kind \IK if `absx` is of type `real`,
    !>                              <li>    type `integer` of the same kind as `absx` if `absx` is of type `integer`,
    !>                          </ol>
    !>                          containing the next integer exponent.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_mathExp, only: getExpNext
    !>
    !>      expNext = getExpNext(absx, base = base)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `1 < base` must hold.<br>
    !>  The condition `0 < absx` must hold.<br>
    !>  \vericons
    !>
    !>  \warning
    !>  Note that setting `base` very near `1` can lead to very large output exponents that overflow the default integer kind \IK.<br>
    !>  For example, `1.0000001**2147483647 = 1.836644819690732e+93`.<br>
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \devnote
    !>  A `real` value of kind \RK32 can represent `integer` values as large as `huge(1_int128) = 170141183460469231731687303715884105727 = 1.70141183E+38 < huge(1._RK32) = 3.40282347E+38 << huge(1._RK64)`.<br>
    !>  One can envision a distant future human society with advanced computers capable of representing higher precision integer value for which \RK32 or \RK64 would be insufficient.<br>
    !>
    !>  \see
    !>  [isIntPow](@ref pm_mathExp::isIntPow)<br>
    !>  [getExpNext](@ref pm_mathExp::getExpNext)<br>
    !>  [getExpPrev](@ref pm_mathExp::getExpPrev)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_mathExp/getExpNext/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_mathExp/getExpNext/main.out.F90
    !>
    !>  \test
    !>  [test_pm_mathExp](@ref test_pm_mathExp)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, April 25, 2015, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface getExpNext

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE elemental module function getExpNext_IK5(absx, base) result(expNext)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpNext_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG), intent(in)            :: absx
        integer(IKG), intent(in), optional  :: base
        integer(IKG)                        :: expNext
    end function
#endif

#if IK4_ENABLED
    PURE elemental module function getExpNext_IK4(absx, base) result(expNext)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpNext_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG), intent(in)            :: absx
        integer(IKG), intent(in), optional  :: base
        integer(IKG)                        :: expNext
    end function
#endif

#if IK3_ENABLED
    PURE elemental module function getExpNext_IK3(absx, base) result(expNext)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpNext_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG), intent(in)            :: absx
        integer(IKG), intent(in), optional  :: base
        integer(IKG)                        :: expNext
    end function
#endif

#if IK2_ENABLED
    PURE elemental module function getExpNext_IK2(absx, base) result(expNext)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpNext_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG), intent(in)            :: absx
        integer(IKG), intent(in), optional  :: base
        integer(IKG)                        :: expNext
    end function
#endif

#if IK1_ENABLED
    PURE elemental module function getExpNext_IK1(absx, base) result(expNext)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpNext_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG), intent(in)            :: absx
        integer(IKG), intent(in), optional  :: base
        integer(IKG)                        :: expNext
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getExpNext_RK5(absx, base) result(expNext)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpNext_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)            :: absx
        real(RKG)   , intent(in), optional  :: base
        integer(IK)                         :: expNext
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getExpNext_RK4(absx, base) result(expNext)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpNext_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)            :: absx
        real(RKG)   , intent(in), optional  :: base
        integer(IK)                         :: expNext
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getExpNext_RK3(absx, base) result(expNext)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpNext_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)            :: absx
        real(RKG)   , intent(in), optional  :: base
        integer(IK)                         :: expNext
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getExpNext_RK2(absx, base) result(expNext)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpNext_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)            :: absx
        real(RKG)   , intent(in), optional  :: base
        integer(IK)                         :: expNext
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getExpNext_RK1(absx, base) result(expNext)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpNext_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)            :: absx
        real(RKG)   , intent(in), optional  :: base
        integer(IK)                         :: expNext
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate the <b>previous integer exponent</b> `expPrev` for the specified input `base` and absolute value `absx = abs(x)`
    !>  such that the condition `base**expPrev <= absx` holds.
    !>
    !>  \details
    !>  See the documentation of [pm_mathExp](@ref pm_mathExp) for more information.<br>
    !>
    !>  \param[in]  absx    :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of either,<br>
    !>                          <ol>
    !>                              <li>    type `integer` of kind \IKALL, or<br>
    !>                              <li>    type `real` of kind \RKALL,<br>
    !>                          </ol>
    !>                          representing the **absolute** value of the number for which the previous integer exponent in specified `base` must be computed.<br>
    !>  \param[in]  base    :   The input positive scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `absx`,
    !>                          representing the base of the exponentiation.<br>
    !>                          (**optional**, default = `2`)
    !>
    !>  \return
    !>  `expPrev`           :   The output scalar (or array of the same rank, shape, and size as other array-like arguments) of
    !>                          <ol>
    !>                              <li>    type `integer` of default kind \IK if `absx` is of type `real`,
    !>                              <li>    type `integer` of the same kind as `absx` if `absx` is of type `integer`,
    !>                          </ol>
    !>                          containing the previous integer exponent.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_mathExp, only: getExpPrev
    !>
    !>      expPrev = getExpPrev(absx, base = base)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `1 < base` must hold.<br>
    !>  The condition `0 < absx` must hold.<br>
    !>  \vericons
    !>
    !>  \warning
    !>  Note that setting `base` very near `1` can lead to very large output exponents that overflow the default integer kind \IK.<br>
    !>  For example, `1.0000001**2147483647 = 1.836644819690732e+93`.<br>
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \devnote
    !>  A `real` value of kind \RK32 can represent `integer` values as large as `huge(1_int128) = 170141183460469231731687303715884105727 = 1.70141183E+38 < huge(1._RK32) = 3.40282347E+38 << huge(1._RK64)`.<br>
    !>  One can envision a distant future human society with advanced computers capable of representing higher precision integer value for which \RK32 or \RK64 would be insufficient.<br>
    !>
    !>  \see
    !>  [isIntPow](@ref pm_mathExp::isIntPow)<br>
    !>  [getExpNext](@ref pm_mathExp::getExpNext)<br>
    !>  [getExpPrev](@ref pm_mathExp::getExpPrev)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_mathExp/getExpPrev/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_mathExp/getExpPrev/main.out.F90
    !>
    !>  \test
    !>  [test_pm_mathExp](@ref test_pm_mathExp)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, April 25, 2015, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface getExpPrev

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE elemental module function getExpPrev_IK5(absx, base) result(expPrev)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpPrev_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG), intent(in)            :: absx
        integer(IKG), intent(in), optional  :: base
        integer(IKG)                        :: expPrev
    end function
#endif

#if IK4_ENABLED
    PURE elemental module function getExpPrev_IK4(absx, base) result(expPrev)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpPrev_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG), intent(in)            :: absx
        integer(IKG), intent(in), optional  :: base
        integer(IKG)                        :: expPrev
    end function
#endif

#if IK3_ENABLED
    PURE elemental module function getExpPrev_IK3(absx, base) result(expPrev)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpPrev_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG), intent(in)            :: absx
        integer(IKG), intent(in), optional  :: base
        integer(IKG)                        :: expPrev
    end function
#endif

#if IK2_ENABLED
    PURE elemental module function getExpPrev_IK2(absx, base) result(expPrev)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpPrev_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG), intent(in)            :: absx
        integer(IKG), intent(in), optional  :: base
        integer(IKG)                        :: expPrev
    end function
#endif

#if IK1_ENABLED
    PURE elemental module function getExpPrev_IK1(absx, base) result(expPrev)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpPrev_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG), intent(in)            :: absx
        integer(IKG), intent(in), optional  :: base
        integer(IKG)                        :: expPrev
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getExpPrev_RK5(absx, base) result(expPrev)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpPrev_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)            :: absx
        real(RKG)   , intent(in), optional  :: base
        integer(IK)                         :: expPrev
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getExpPrev_RK4(absx, base) result(expPrev)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpPrev_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)            :: absx
        real(RKG)   , intent(in), optional  :: base
        integer(IK)                         :: expPrev
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getExpPrev_RK3(absx, base) result(expPrev)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpPrev_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)            :: absx
        real(RKG)   , intent(in), optional  :: base
        integer(IK)                         :: expPrev
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getExpPrev_RK2(absx, base) result(expPrev)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpPrev_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)            :: absx
        real(RKG)   , intent(in), optional  :: base
        integer(IK)                         :: expPrev
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getExpPrev_RK1(absx, base) result(expPrev)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpPrev_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)            :: absx
        real(RKG)   , intent(in), optional  :: base
        integer(IK)                         :: expPrev
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathExp ! LCOV_EXCL_LINE