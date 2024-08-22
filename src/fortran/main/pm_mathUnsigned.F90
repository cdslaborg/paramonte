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
!>  This module contains procedures and generic interfaces for various operations with positive integers
!>  with results that have the same binary representation as an unsigned integer.<br>
!>  Such operations (like addition or subtraction) would normally cause runtime overflow errors within the default Fortran environment.<br>
!>
!>  \details
!>  Consider the simple case of a 4-bits unsigned integer kind.<br>
!>  The binary representations of all possible values by such an integer are shown the left plot of the figure below.<br>
!>  Now, consider the case of a 4-bits signed integer kind (as illustrated in the right plot of the figure below).<br>
!>
!>  \htmlonly
!>      <img src="pm_mathUnsigned@4bits.png" style="width:50%;">
!>  \endhtmlonly
!>
!>  Because all Fortran integer kinds are by default signed, the addition of any two (signed)
!>  integer values are guaranteed to not cause overflow if added as unsigned integers.<br>
!>  The trick to compute the binary representation of such overflowed signed additions is
!>  therefore to compute the negative number in the signed integer space that corresponds
!>  to the unsigned addition of the two positive numbers.<br>
!>  For example, consider the operation `4 + 6 = 10` where all numbers are 4-bit unsigned integers.<br>
!>  This addition would be invalid and cause runtime overflow if `4` and `6` are 4-bits *signed* integers as illustrated in the figure above.<br>
!>  However, the binary representation of the result (`10`) of the unsigned addition operation can be readily computed
!>  by first rotating `4` and `6` by 180 degrees on the circular representation of the numbers in the figure to get `-4` and `-2`.<br>
!>  Then, adding the rotated numbers yields the correct binary representation for the result of `4 + 6 = 10` as an unsigned integer addition.<br>
!>  This rotation operation is equivalent to negating the binary signs of the two addition operands before adding them.<br>
!>  Note that an addition operation can never overflow if the two operands are of opposite signs.<br>
!>
!>  \warning
!>  The procedures of this module work based on the plausible assumption that 
!>  integers are represented in [two's complement format](https://en.wikipedia.org/wiki/Two%27s_complement).<br>
!>
!>  \note
!>  The sum of two integer values with opposite signs will never cause arithmetic overflow.<br>
!>
!>  \see
!>  [pm_except](@ref pm_except)<br>
!>
!>  \test
!>  [test_pm_mathUnsigned](@ref test_pm_mathUnsigned)
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_mathUnsigned

    use pm_kind, only: SK, IK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathUnsigned"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a (possibly overflowed) signed integer that is
    !>  the result of adding the two input (non-negative) integers without runtime overflow error.
    !>
    !>  \details
    !>  If the addition operation does not lead to overflow, then the sum matches the outcome of regular addition operation.<br>
    !>  If the addition operation does lead to overflow, then bit pattern of the result will match that of the result of the
    !>  unsigned addition of the two input integers.<br>
    !>  See the documentation of [pm_mathUnsigned](@ref pm_mathUnsigned) for more information.<br>
    !>
    !>  \param[in]  a   :   The input scalar, or array of the same rank and shape as other array-like arguments, of type `integer` of kind \IKALL,
    !>                      containing a (non-negative) signed integer to be added to another input signed integer `b` without causing overflow.
    !>  \param[in]  b   :   The input scalar, or array of the same rank and shape as other array-like arguments, of the same type and kind as `a`,
    !>                      containing a (non-negative) signed integer to be added to another input signed integer `a` without causing overflow.
    !>
    !>  \return
    !>  `sum`           :   The output scalar, or array of the same rank and shape as other array-like arguments, of the same type and kind as `a`,
    !>                      containing the sum of the two input (non-negative) signed integers without causing runtime overflow error.<br>
    !>
    !>  \interface{uadd}
    !>  \code{.F90}
    !>
    !>      use pm_mathUnsigned, only: operator(.uadd.)
    !>
    !>      c = a .uadd. b
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= a` must hold for the corresponding input arguments.<br>
    !>  The condition `0 <= b` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [pm_distUnif](@ref pm_distUnif)<br>
    !>
    !>  \example{uadd}
    !>  \include{lineno} example/pm_mathUnsigned/uadd/main.F90
    !>  \compilef{uadd}
    !>  \output{uadd}
    !>  \include{lineno} example/pm_mathUnsigned/uadd/main.out.F90
    !>
    !>  \test
    !>  [test_pm_mathUnsigned](@ref test_pm_mathUnsigned)
    !>
    !>  \final{uadd}
    !>
    !>  \author
    !>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface operator(.uadd.)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE elemental module function uadd_IK5(a, b) result(sum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: uadd_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG), intent(in)            :: a, b
        integer(IKG)                        :: sum
    end function
#endif

#if IK4_ENABLED
    PURE elemental module function uadd_IK4(a, b) result(sum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: uadd_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG), intent(in)            :: a, b
        integer(IKG)                        :: sum
    end function
#endif

#if IK3_ENABLED
    PURE elemental module function uadd_IK3(a, b) result(sum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: uadd_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG), intent(in)            :: a, b
        integer(IKG)                        :: sum
    end function
#endif

#if IK2_ENABLED
    PURE elemental module function uadd_IK2(a, b) result(sum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: uadd_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG), intent(in)            :: a, b
        integer(IKG)                        :: sum
    end function
#endif

#if IK1_ENABLED
    PURE elemental module function uadd_IK1(a, b) result(sum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: uadd_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG), intent(in)            :: a, b
        integer(IKG)                        :: sum
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathUnsigned