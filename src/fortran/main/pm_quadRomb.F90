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
!>  This module contains classes and procedures to perform numerical integrations.
!>
!>  \test
!>  [test_pm_distExp](@ref test_pm_distExp)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 4:41 AM, Michigan<br>
!>  \JoshuaOsborne, May 28, 2020, 9:06 PM, Arlington, TX<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_quadRomb

    use pm_kind, only: IK, RK, SK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_quadRomb"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the indicator type for generating instances of objects that indicate the integration interval is open.
    !>
    !>  \details
    !>  This is an empty derived type that exists solely for generating unique objects that are distinguishable
    !>  as input arguments to procedures under the generic interface [getQuadRomb](@ref pm_quadRomb::getQuadRomb).
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_quadRomb, only: open_type
    !>      type(open_type) :: Open
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [lbis_type](@ref pm_quadRomb::lbis_type)<br>
    !>  [nexp_type](@ref pm_quadRomb::nexp_type)<br>
    !>  [open_type](@ref pm_quadRomb::open_type)<br>
    !>  [pexp_type](@ref pm_quadRomb::pexp_type)<br>
    !>  [pwrl_type](@ref pm_quadRomb::pwrl_type)<br>
    !>  [ubis_type](@ref pm_quadRomb::ubis_type)<br>
    !>  [getQuadRomb](@ref pm_quadRomb::getQuadRomb)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_quadRomb/open_type/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_quadRomb/open_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_quadRomb](@ref test_pm_quadRomb)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: open_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the indicator type for generating instances of objects that indicate the integration interval is open \f$(a, b)\f$ and,
    !>  the intervals should be spaced assuming an integrand that behaves like,
    !>  <ol>
    !>      <li>    a  decreasing Power-Law (PWL) on a positive support \f$(a > 0, b > 0)\f$, such that the upper limit of integration is allowed to be \f$b = +\infty\f$, or
    !>      <li>    an increasing Power-Law (PWL) on a negative support \f$(a < 0, b < 0)\f$, such that the lower limit of integration is allowed to be \f$a = -\infty\f$.
    !>  </ol>
    !>
    !>  \details
    !>  This is an empty derived type that exists solely for generating unique objects that are distinguishable
    !>  as input arguments to procedures under the generic interface [getQuadRomb](@ref pm_quadRomb::getQuadRomb).<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_quadRomb, only: pwrl_type
    !>      type(pwrl_type) :: PWRL
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [lbis_type](@ref pm_quadRomb::lbis_type)<br>
    !>  [nexp_type](@ref pm_quadRomb::nexp_type)<br>
    !>  [open_type](@ref pm_quadRomb::open_type)<br>
    !>  [pexp_type](@ref pm_quadRomb::pexp_type)<br>
    !>  [pwrl_type](@ref pm_quadRomb::pwrl_type)<br>
    !>  [ubis_type](@ref pm_quadRomb::ubis_type)<br>
    !>  [getQuadRomb](@ref pm_quadRomb::getQuadRomb)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_quadRomb/pwrl_type/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_quadRomb/pwrl_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_quadRomb](@ref test_pm_quadRomb)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: pwrl_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the indicator type for generating instances of objects that indicate the integration interval is open and,
    !>  the intervals should be spaced assuming an integrand that behaves like a Negative-Exponent Exponential (NEXP),
    !>  such that the upper limit of integration is allowed to be \f$b = +\infty\f$.
    !>
    !>  \details
    !>  This is an empty derived type that exists solely for generating unique objects that are distinguishable
    !>  as input arguments to procedures under the generic interface [getQuadRomb](@ref pm_quadRomb::getQuadRomb).<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_quadRomb, only: nexp_type
    !>      type(nexp_type) :: NEXP
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [lbis_type](@ref pm_quadRomb::lbis_type)<br>
    !>  [nexp_type](@ref pm_quadRomb::nexp_type)<br>
    !>  [open_type](@ref pm_quadRomb::open_type)<br>
    !>  [pexp_type](@ref pm_quadRomb::pexp_type)<br>
    !>  [pwrl_type](@ref pm_quadRomb::pwrl_type)<br>
    !>  [ubis_type](@ref pm_quadRomb::ubis_type)<br>
    !>  [getQuadRomb](@ref pm_quadRomb::getQuadRomb)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_quadRomb/nexp_type/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_quadRomb/nexp_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_quadRomb](@ref test_pm_quadRomb)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: nexp_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the indicator type for generating instances of objects that indicate the integration interval is open and,
    !>  the intervals should be spaced assuming an integrand that behaves like a Positive-Exponent Exponential (PEXP),
    !>  such that the lower limit of integration is allowed to be \f$a = -\infty\f$.
    !>
    !>  \details
    !>  This is an empty derived type that exists solely for generating unique objects that are distinguishable
    !>  as input arguments to procedures under the generic interface [getQuadRomb](@ref pm_quadRomb::getQuadRomb).<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_quadRomb, only: pexp_type
    !>      type(pexp_type) :: PEXP
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [lbis_type](@ref pm_quadRomb::lbis_type)<br>
    !>  [nexp_type](@ref pm_quadRomb::nexp_type)<br>
    !>  [open_type](@ref pm_quadRomb::open_type)<br>
    !>  [pexp_type](@ref pm_quadRomb::pexp_type)<br>
    !>  [pwrl_type](@ref pm_quadRomb::pwrl_type)<br>
    !>  [ubis_type](@ref pm_quadRomb::ubis_type)<br>
    !>  [getQuadRomb](@ref pm_quadRomb::getQuadRomb)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_quadRomb/pexp_type/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_quadRomb/pexp_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_quadRomb](@ref test_pm_quadRomb)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: pexp_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the indicator type for generating instances of objects that indicate the integration interval is open and,
    !>  the integrand has an Integrable square-root type of Singularity at the finite Lower Bound of integration (LBIS).
    !>
    !>  \details
    !>  This is an empty derived type that exists solely for generating unique objects that are distinguishable
    !>  as input arguments to procedures under the generic interface [getQuadRomb](@ref pm_quadRomb::getQuadRomb).<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_quadRomb, only: lbis_type
    !>      type(lbis_type) :: LBIS
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [lbis_type](@ref pm_quadRomb::lbis_type)<br>
    !>  [nexp_type](@ref pm_quadRomb::nexp_type)<br>
    !>  [open_type](@ref pm_quadRomb::open_type)<br>
    !>  [pexp_type](@ref pm_quadRomb::pexp_type)<br>
    !>  [pwrl_type](@ref pm_quadRomb::pwrl_type)<br>
    !>  [ubis_type](@ref pm_quadRomb::ubis_type)<br>
    !>  [getQuadRomb](@ref pm_quadRomb::getQuadRomb)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_quadRomb/lbis_type/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_quadRomb/lbis_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_quadRomb](@ref test_pm_quadRomb)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: lbis_type
        real :: exponent = 0.5
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the indicator type for generating instances of objects that indicate the integration interval is open and,
    !>  the integrand has an Integrable square-root type of Singularity at the finite Lower Bound of integration (LBIS).
    !>
    !>  \details
    !>  This is an empty derived type that exists solely for generating unique objects that are distinguishable
    !>  as input arguments to procedures under the generic interface [getQuadRomb](@ref pm_quadRomb::getQuadRomb).<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_quadRomb, only: ubis_type
    !>      type(ubis_type) :: UBIS
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [lbis_type](@ref pm_quadRomb::lbis_type)<br>
    !>  [nexp_type](@ref pm_quadRomb::nexp_type)<br>
    !>  [open_type](@ref pm_quadRomb::open_type)<br>
    !>  [pexp_type](@ref pm_quadRomb::pexp_type)<br>
    !>  [pwrl_type](@ref pm_quadRomb::pwrl_type)<br>
    !>  [ubis_type](@ref pm_quadRomb::ubis_type)<br>
    !>  [getQuadRomb](@ref pm_quadRomb::getQuadRomb)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_quadRomb/ubis_type/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_quadRomb/ubis_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_quadRomb](@ref test_pm_quadRomb)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: ubis_type
        real :: exponent = 0.5
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the integral of the input function `getFunc()` in the closed range `[lb, ub]` using the **Adaptive Romberg extrapolation method**.
    !>
    !>  \details
    !>  This Romberg integration method is quite powerful for sufficiently smooth (e.g., analytic) integrands `getFunc()`,
    !>  integrated over bounded intervals which contain no singularities, and where the end points are also nonsingular.
    !>
    !>  \param      getFunc :   The input function to be integrated (i.e., the integrand).
    !>                              -#  On entry, it must take an input scalar of the same type and kind as `quadRomb`.<br>
    !>                              -#  On exit, it must generate an input scalar of the same type and kind as `quadRomb`, representing the corresponding function value.<br>
    !>                          The following illustrates the general interface of `getFunc`:
    !>                          \code{.F90}
    !>
    !>                              function getFunc(x) result(func)
    !>                                  use pm_kind, only: RK => RKG
    !>                                  real(RK)    , intent(in)    :: x
    !>                                  real(RK)                    :: func
    !>                              end function
    !>
    !>                          \endcode
    !>                          where `RKG` refers to any desired `real` kind supported by the processor.<br>
    !>  \param[in]  lb      :   The input scalar of the same type and kind as the output `quadRomb`, containing the lower bound of the integration.
    !>  \param[in]  ub      :   The input scalar of the same type and kind as the output `quadRomb`, containing the upper bound of the integration.
    !>  \param[in]  tol     :   The input scalar of the same type and kind as the output `quadRomb`, containing the *relative* error the integration.<br>
    !>                          The algorithm converges if \f$\ms{relerr} ~(~\equiv |\Delta\ms{quadRomb}|~) ~\leq~ \ms{tol} \times | \ms{quadRomb} |\f$.<br>
    !>                          Note that `tol > epsilon(0._RKG)` must hold at all times for integration to converge. Here `RKG` is the desired `real` kind of the output.<br>
    !>                          Ideally, set `tol` to a value such that `tol < epsilon(0._RKG) * 100` holds to ensure convergence.<br>
    !>  \param[in]  nref    :   The input scalar `integer` of default kind \IK, representing the number of refinements to be used in the Romberg method.<br>
    !>                          Think of `nref` as the maximum possible degree of the polynomial extrapolation used for approximating the integral at any stage.
    !>                          <ul>
    !>                              <li>    The smaller values of `nref` can delay an accurate estimation of the integral via the Romberg method.
    !>                              <li>    The larger values of `nref` can delay the first estimation of the integral by requiring more function evaluations.
    !>                              <li>    The computational precision of the `real` kind used in this procedure imposes an upper limit on the value of `nref`.<br>
    !>                                      The maximum value for `nref` is roughly equal to `int(log(epsilon(1._RKG)) / log(0.25))`
    !>                                      with `RKG` representing the `real` kind used for the integration.
    !>                              <li>    The maximum `nref` for `real32`, `real64`, `real128` are respectively `12`, `26`, `56`.<br>
    !>                              <li>    If the specified `nref` is larger than the maximum possible value, the integration will fail to converge.<br>
    !>                              <li>    The number `nref = 2` corresponds to the famous Simpson integration rule.
    !>                              <li>    A number between 4-6 is frequently a reasonable choice.
    !>                          </ul>
    !>  \param[out] relerr  :   The output scalar of the same type and kind as the output `quadRomb` containing the final estimated relative error in the result.<br>
    !>                          By definition, this is **always a positive** value **if the integration converges**.<br>
    !>                          **Specify the `relerr` optional output argument to monitor convergence**.<br>
    !>                          If `relerr < 0.`, then the integration has failed to converge.<br>
    !>                          (**optional**. If missing and the integration fails to converge, the program will halt by calling `error stop`.)<br>
    !>  \param[out] neval   :   The output scalar `integer` of default kind \IK, representing the number of function evaluations made within the integrator.<br>
    !>                          (**optional**. It can be present <b>if and only if</b> `relerr` argument is also present.)
    !>
    !>  \return
    !>  `quadRomb`          :   The output scalar `real` of kind \RKALL, containing the integration result.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_quadRomb, only: getQuadRomb
    !>
    !>      quadRomb = getQuadRomb(getFunc, lb, ub, tol, nref)
    !>      quadRomb = getQuadRomb(getFunc, lb, ub, tol, nref, relerr)
    !>      quadRomb = getQuadRomb(getFunc, lb, ub, tol, nref, relerr, neval)
    !>
    !>      quadRomb = getQuadRomb(getFunc, lb, ub, tol, nref, interval)
    !>      quadRomb = getQuadRomb(getFunc, lb, ub, tol, nref, interval, relerr)
    !>      quadRomb = getQuadRomb(getFunc, lb, ub, tol, nref, interval, relerr, neval)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The procedures of this generic interface will behave differently if the integration fails to converge:
    !>  <ul>
    !>      <li>    If the optional `relerr` output argument is missing, the program will halt by calling `error stop` upon integration convergence failure.
    !>      <li>    If the optional `relerr` output argument is present, the program will return `relerr < 0.` to indicate integration convergence failure.
    !>  </ul>
    !>
    !>  \see
    !>  [getQuadErr](@ref pm_quadPack::getQuadErr)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_quadRomb/getQuadRomb/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_quadRomb/getQuadRomb/main.out.F90
    !>
    !>  \test
    !>  [test_pm_quadRomb](@ref test_pm_quadRomb)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    !>  \JoshuaOsborne, May 28, 2020, 8:58 PM, Arlington, TX
    interface getQuadRomb

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive impure module function getQR_Clos_EM_NM_RK5(getFunc, lb, ub, tol, nref) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Clos_EM_NM_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK4_ENABLED
    recursive impure module function getQR_Clos_EM_NM_RK4(getFunc, lb, ub, tol, nref) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Clos_EM_NM_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK3_ENABLED
    recursive impure module function getQR_Clos_EM_NM_RK3(getFunc, lb, ub, tol, nref) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Clos_EM_NM_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK2_ENABLED
    recursive impure module function getQR_Clos_EM_NM_RK2(getFunc, lb, ub, tol, nref) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Clos_EM_NM_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK1_ENABLED
    recursive impure module function getQR_Clos_EM_NM_RK1(getFunc, lb, ub, tol, nref) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Clos_EM_NM_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        real(RKG)                               :: quadRomb
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive impure module function getQR_Clos_EP_NM_RK5(getFunc, lb, ub, tol, nref, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Clos_EP_NM_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK4_ENABLED
    recursive impure module function getQR_Clos_EP_NM_RK4(getFunc, lb, ub, tol, nref, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Clos_EP_NM_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK3_ENABLED
    recursive impure module function getQR_Clos_EP_NM_RK3(getFunc, lb, ub, tol, nref, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Clos_EP_NM_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK2_ENABLED
    recursive impure module function getQR_Clos_EP_NM_RK2(getFunc, lb, ub, tol, nref, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Clos_EP_NM_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK1_ENABLED
    recursive impure module function getQR_Clos_EP_NM_RK1(getFunc, lb, ub, tol, nref, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Clos_EP_NM_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive impure module function getQR_Clos_EP_NP_RK5(getFunc, lb, ub, tol, nref, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Clos_EP_NP_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK4_ENABLED
    recursive impure module function getQR_Clos_EP_NP_RK4(getFunc, lb, ub, tol, nref, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Clos_EP_NP_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK3_ENABLED
    recursive impure module function getQR_Clos_EP_NP_RK3(getFunc, lb, ub, tol, nref, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Clos_EP_NP_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK2_ENABLED
    recursive impure module function getQR_Clos_EP_NP_RK2(getFunc, lb, ub, tol, nref, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Clos_EP_NP_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK1_ENABLED
    recursive impure module function getQR_Clos_EP_NP_RK1(getFunc, lb, ub, tol, nref, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Clos_EP_NP_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive impure module function getQR_Open_EM_NM_RK5(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Open_EM_NM_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(open_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK4_ENABLED
    recursive impure module function getQR_Open_EM_NM_RK4(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Open_EM_NM_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(open_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK3_ENABLED
    recursive impure module function getQR_Open_EM_NM_RK3(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Open_EM_NM_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(open_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK2_ENABLED
    recursive impure module function getQR_Open_EM_NM_RK2(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Open_EM_NM_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(open_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK1_ENABLED
    recursive impure module function getQR_Open_EM_NM_RK1(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Open_EM_NM_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(open_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive impure module function getQR_Open_EP_NM_RK5(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Open_EP_NM_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(open_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK4_ENABLED
    recursive impure module function getQR_Open_EP_NM_RK4(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Open_EP_NM_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(open_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK3_ENABLED
    recursive impure module function getQR_Open_EP_NM_RK3(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Open_EP_NM_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(open_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK2_ENABLED
    recursive impure module function getQR_Open_EP_NM_RK2(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Open_EP_NM_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(open_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK1_ENABLED
    recursive impure module function getQR_Open_EP_NM_RK1(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Open_EP_NM_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(open_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive impure module function getQR_Open_EP_NP_RK5(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Open_EP_NP_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(open_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK4_ENABLED
    recursive impure module function getQR_Open_EP_NP_RK4(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Open_EP_NP_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(open_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK3_ENABLED
    recursive impure module function getQR_Open_EP_NP_RK3(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Open_EP_NP_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(open_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK2_ENABLED
    recursive impure module function getQR_Open_EP_NP_RK2(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Open_EP_NP_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(open_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK1_ENABLED
    recursive impure module function getQR_Open_EP_NP_RK1(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_Open_EP_NP_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(open_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive impure module function getQR_PWRL_EM_NM_RK5(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PWRL_EM_NM_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pwrl_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK4_ENABLED
    recursive impure module function getQR_PWRL_EM_NM_RK4(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PWRL_EM_NM_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pwrl_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK3_ENABLED
    recursive impure module function getQR_PWRL_EM_NM_RK3(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PWRL_EM_NM_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pwrl_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK2_ENABLED
    recursive impure module function getQR_PWRL_EM_NM_RK2(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PWRL_EM_NM_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pwrl_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK1_ENABLED
    recursive impure module function getQR_PWRL_EM_NM_RK1(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PWRL_EM_NM_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pwrl_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive impure module function getQR_PWRL_EP_NM_RK5(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PWRL_EP_NM_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pwrl_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK4_ENABLED
    recursive impure module function getQR_PWRL_EP_NM_RK4(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PWRL_EP_NM_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pwrl_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK3_ENABLED
    recursive impure module function getQR_PWRL_EP_NM_RK3(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PWRL_EP_NM_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pwrl_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK2_ENABLED
    recursive impure module function getQR_PWRL_EP_NM_RK2(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PWRL_EP_NM_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pwrl_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK1_ENABLED
    recursive impure module function getQR_PWRL_EP_NM_RK1(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PWRL_EP_NM_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pwrl_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive impure module function getQR_PWRL_EP_NP_RK5(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PWRL_EP_NP_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pwrl_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK4_ENABLED
    recursive impure module function getQR_PWRL_EP_NP_RK4(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PWRL_EP_NP_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pwrl_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK3_ENABLED
    recursive impure module function getQR_PWRL_EP_NP_RK3(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PWRL_EP_NP_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pwrl_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK2_ENABLED
    recursive impure module function getQR_PWRL_EP_NP_RK2(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PWRL_EP_NP_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pwrl_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK1_ENABLED
    recursive impure module function getQR_PWRL_EP_NP_RK1(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PWRL_EP_NP_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pwrl_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive impure module function getQR_NEXP_EM_NM_RK5(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_NEXP_EM_NM_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(nexp_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK4_ENABLED
    recursive impure module function getQR_NEXP_EM_NM_RK4(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_NEXP_EM_NM_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(nexp_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK3_ENABLED
    recursive impure module function getQR_NEXP_EM_NM_RK3(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_NEXP_EM_NM_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(nexp_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK2_ENABLED
    recursive impure module function getQR_NEXP_EM_NM_RK2(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_NEXP_EM_NM_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(nexp_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK1_ENABLED
    recursive impure module function getQR_NEXP_EM_NM_RK1(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_NEXP_EM_NM_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(nexp_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive impure module function getQR_NEXP_EP_NM_RK5(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_NEXP_EP_NM_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(nexp_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK4_ENABLED
    recursive impure module function getQR_NEXP_EP_NM_RK4(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_NEXP_EP_NM_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(nexp_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK3_ENABLED
    recursive impure module function getQR_NEXP_EP_NM_RK3(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_NEXP_EP_NM_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(nexp_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK2_ENABLED
    recursive impure module function getQR_NEXP_EP_NM_RK2(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_NEXP_EP_NM_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(nexp_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK1_ENABLED
    recursive impure module function getQR_NEXP_EP_NM_RK1(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_NEXP_EP_NM_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(nexp_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive impure module function getQR_NEXP_EP_NP_RK5(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_NEXP_EP_NP_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(nexp_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK4_ENABLED
    recursive impure module function getQR_NEXP_EP_NP_RK4(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_NEXP_EP_NP_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(nexp_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK3_ENABLED
    recursive impure module function getQR_NEXP_EP_NP_RK3(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_NEXP_EP_NP_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(nexp_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK2_ENABLED
    recursive impure module function getQR_NEXP_EP_NP_RK2(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_NEXP_EP_NP_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(nexp_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK1_ENABLED
    recursive impure module function getQR_NEXP_EP_NP_RK1(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_NEXP_EP_NP_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(nexp_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive impure module function getQR_PEXP_EM_NM_RK5(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PEXP_EM_NM_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pexp_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK4_ENABLED
    recursive impure module function getQR_PEXP_EM_NM_RK4(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PEXP_EM_NM_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pexp_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK3_ENABLED
    recursive impure module function getQR_PEXP_EM_NM_RK3(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PEXP_EM_NM_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pexp_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK2_ENABLED
    recursive impure module function getQR_PEXP_EM_NM_RK2(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PEXP_EM_NM_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pexp_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK1_ENABLED
    recursive impure module function getQR_PEXP_EM_NM_RK1(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PEXP_EM_NM_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pexp_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive impure module function getQR_PEXP_EP_NM_RK5(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PEXP_EP_NM_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pexp_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK4_ENABLED
    recursive impure module function getQR_PEXP_EP_NM_RK4(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PEXP_EP_NM_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pexp_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK3_ENABLED
    recursive impure module function getQR_PEXP_EP_NM_RK3(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PEXP_EP_NM_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pexp_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK2_ENABLED
    recursive impure module function getQR_PEXP_EP_NM_RK2(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PEXP_EP_NM_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pexp_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK1_ENABLED
    recursive impure module function getQR_PEXP_EP_NM_RK1(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PEXP_EP_NM_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pexp_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive impure module function getQR_PEXP_EP_NP_RK5(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PEXP_EP_NP_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pexp_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK4_ENABLED
    recursive impure module function getQR_PEXP_EP_NP_RK4(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PEXP_EP_NP_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pexp_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK3_ENABLED
    recursive impure module function getQR_PEXP_EP_NP_RK3(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PEXP_EP_NP_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pexp_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK2_ENABLED
    recursive impure module function getQR_PEXP_EP_NP_RK2(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PEXP_EP_NP_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pexp_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK1_ENABLED
    recursive impure module function getQR_PEXP_EP_NP_RK1(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_PEXP_EP_NP_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(pexp_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive impure module function getQR_LBIS_EM_NM_RK5(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_LBIS_EM_NM_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(lbis_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK4_ENABLED
    recursive impure module function getQR_LBIS_EM_NM_RK4(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_LBIS_EM_NM_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(lbis_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK3_ENABLED
    recursive impure module function getQR_LBIS_EM_NM_RK3(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_LBIS_EM_NM_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(lbis_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK2_ENABLED
    recursive impure module function getQR_LBIS_EM_NM_RK2(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_LBIS_EM_NM_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(lbis_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK1_ENABLED
    recursive impure module function getQR_LBIS_EM_NM_RK1(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_LBIS_EM_NM_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(lbis_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive impure module function getQR_LBIS_EP_NM_RK5(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_LBIS_EP_NM_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(lbis_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK4_ENABLED
    recursive impure module function getQR_LBIS_EP_NM_RK4(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_LBIS_EP_NM_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(lbis_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK3_ENABLED
    recursive impure module function getQR_LBIS_EP_NM_RK3(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_LBIS_EP_NM_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(lbis_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK2_ENABLED
    recursive impure module function getQR_LBIS_EP_NM_RK2(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_LBIS_EP_NM_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(lbis_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK1_ENABLED
    recursive impure module function getQR_LBIS_EP_NM_RK1(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_LBIS_EP_NM_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(lbis_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive impure module function getQR_LBIS_EP_NP_RK5(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_LBIS_EP_NP_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(lbis_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK4_ENABLED
    recursive impure module function getQR_LBIS_EP_NP_RK4(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_LBIS_EP_NP_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(lbis_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK3_ENABLED
    recursive impure module function getQR_LBIS_EP_NP_RK3(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_LBIS_EP_NP_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(lbis_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK2_ENABLED
    recursive impure module function getQR_LBIS_EP_NP_RK2(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_LBIS_EP_NP_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(lbis_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK1_ENABLED
    recursive impure module function getQR_LBIS_EP_NP_RK1(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_LBIS_EP_NP_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(lbis_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive impure module function getQR_UBIS_EM_NM_RK5(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_UBIS_EM_NM_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(ubis_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK4_ENABLED
    recursive impure module function getQR_UBIS_EM_NM_RK4(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_UBIS_EM_NM_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(ubis_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK3_ENABLED
    recursive impure module function getQR_UBIS_EM_NM_RK3(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_UBIS_EM_NM_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(ubis_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK2_ENABLED
    recursive impure module function getQR_UBIS_EM_NM_RK2(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_UBIS_EM_NM_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(ubis_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK1_ENABLED
    recursive impure module function getQR_UBIS_EM_NM_RK1(getFunc, lb, ub, tol, nref, interval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_UBIS_EM_NM_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(ubis_type)                         :: interval
        real(RKG)                               :: quadRomb
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive impure module function getQR_UBIS_EP_NM_RK5(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_UBIS_EP_NM_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(ubis_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK4_ENABLED
    recursive impure module function getQR_UBIS_EP_NM_RK4(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_UBIS_EP_NM_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(ubis_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK3_ENABLED
    recursive impure module function getQR_UBIS_EP_NM_RK3(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_UBIS_EP_NM_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(ubis_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK2_ENABLED
    recursive impure module function getQR_UBIS_EP_NM_RK2(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_UBIS_EP_NM_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(ubis_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK1_ENABLED
    recursive impure module function getQR_UBIS_EP_NM_RK1(getFunc, lb, ub, tol, nref, interval, relerr) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_UBIS_EP_NM_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(ubis_type)                         :: interval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive impure module function getQR_UBIS_EP_NP_RK5(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_UBIS_EP_NP_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(ubis_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK4_ENABLED
    recursive impure module function getQR_UBIS_EP_NP_RK4(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_UBIS_EP_NP_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(ubis_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK3_ENABLED
    recursive impure module function getQR_UBIS_EP_NP_RK3(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_UBIS_EP_NP_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(ubis_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK2_ENABLED
    recursive impure module function getQR_UBIS_EP_NP_RK2(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_UBIS_EP_NP_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(ubis_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

#if RK1_ENABLED
    recursive impure module function getQR_UBIS_EP_NP_RK1(getFunc, lb, ub, tol, nref, interval, relerr, neval) result(quadRomb)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQR_UBIS_EP_NP_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                    :: getFunc
        real(RKG)   , intent(in)                :: lb, ub, tol
        integer(IK) , intent(in)                :: nref
        type(ubis_type)                         :: interval
        integer(IK) , intent(out)               :: neval
        real(RKG)   , intent(out)               :: relerr
        real(RKG)                               :: quadRomb
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    !>  \brief
!    !>  Generate and return the \f$n\f$th stage refinement of the integral of the input function `getFunc()` in the closed range `[lb, ub]`
!    !>  using the **Extended Trapezoidal Rule**.
!    !>
!    !>  \details
!    !>  The coarsest implementation of the trapezoidal rule is to average the function at its endpoints `lb` and `ub`.<br>
!    !>  The first stage of the refinement is to add to this average the value of the function at the halfway point.<br>
!    !>  The second stage of refinement is to add the values at the `1/4` and `3/4` points.<br>
!    !>  This procedure continues until the user is satisfied with the achieved accuracy of the integral.<br>
!    !>  With each new sequential call to `setQuadTrap()`, the accuracy of the output integral will be improved
!    !>  by adding \f42^{\ms{nref} - 2} extra function evaluation points between the integration bounds.<br>
!    !>
!    !>  \param[in]      getFunc     :   The input function to be integrated (i.e., the integrand).
!    !>                                      -#  On entry, it must take an input scalar of the same type and kind as `quadTrap`.<br>
!    !>                                      -#  On exit, it must generate an input scalar of the same type and kind as `quadTrap`, representing the corresponding function value.<br>
!    !>  \param[in]      lb          :   The input scalar of the same type and kind as the output `quadTrap`, containing the lower bound of the integration.
!    !>  \param[in]      ub          :   The input scalar of the same type and kind as the output `quadTrap`, containing the upper bound of the integration.
!    !>  \param[in]      nref        :   The input scalar `integer` of default kind \IK, representing the integration refinement stage.<br>
!    !>                                  The value of `nref` is supposed to be `1` in the first call to `setQuadTrap()` and subsequently increase by one in each subsequent call.<br>
!    !>  \param[inout]   quadTrap    :   The input/output scalar `real` of kind \RKALL, representing the integration result.<br>
!    !>                                  On input, it must contain the integration result at previous refinements.<br>
!    !>                                  On output, it is overwritten with the most recent integral estimate via the most recent function evaluations.<br>
!    !>                                  On the first call to `setQuadTrap()`, it must be set to `0`.<br>
!    !>  \param[out]     neval       :   The output scalar `integer` of default kind \IK, representing the number of function evaluations made within the integrator.<br>
!    !>                                  (**optional**. It can be present <b>if and only if</b> `relerr` is also present.)
!    !>  \param[out]     relerr      :   The output scalar of the same type and kind as the output `quadTrap` containing the final estimated relative error in the result.<br>
!    !>                                  By definition, this is always smaller than the specified input `tol`.<br>
!    !>                                  (**optional**. It can be present <b>if and only if</b> `neval` is also present.)
!    !>
!    !>  \interface
!    !>  \code{.F90}
!    !>
!    !>      use pm_quadRomb, only: setQuadTrap
!    !>
!    !>      call setQuadTrap(getFunc, lb, ub, nref, quadTrap)
!    !>      call setQuadTrap(getFunc, lb, ub, nref, quadTrap, relerr, neval)
!    !>
!    !>  \endcode
!    !>
!    !>  \see
!    !>  [getQuadRombOpen](@ref pm_quadRomb::getQuadRombOpen)<br>
!    !>
!    !>  \example
!    !>  \include{lineno} example/pm_quadRomb/setQuadTrap/main.F90
!    !>  \compilef
!    !>  \output
!    !>  \include{lineno} example/pm_quadRomb/setQuadTrap/main.out.F90
!    !>  \postproc
!    !>  \include{lineno} example/pm_quadRomb/setQuadTrap/main.py
!    !>  \vis
!    !>  \image html pm_quadRomb/setQuadTrap/setQuadTrap.png width=700
!    !>
!    !>  \test
!    !>  [test_pm_quadRomb](@ref test_pm_quadRomb)
!    !>
!    !>  \final
!    !>
!    !>  \author
!    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
!    interface setQuadTrap
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    recursive impure module subroutine setQuadTrap_RK5(getFunc, lb, ub, nref, quadTrap)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setQuadTrap_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        procedure(real(RKG))                    :: getFunc
!        real(RKG)   , intent(in)                :: lb, ub
!        integer(IK) , intent(in)                :: nref
!        real(RKG)                               :: quadTrap
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    recursive impure module subroutine setQuadTrap_RK4(getFunc, lb, ub, nref, quadTrap)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setQuadTrap_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        procedure(real(RKG))                    :: getFunc
!        real(RKG)   , intent(in)                :: lb, ub
!        integer(IK) , intent(in)                :: nref
!        real(RKG)                               :: quadTrap
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    recursive impure module subroutine setQuadTrap_RK3(getFunc, lb, ub, nref, quadTrap)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setQuadTrap_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        procedure(real(RKG))                    :: getFunc
!        real(RKG)   , intent(in)                :: lb, ub
!        integer(IK) , intent(in)                :: nref
!        real(RKG)                               :: quadTrap
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    recursive impure module subroutine setQuadTrap_RK2(getFunc, lb, ub, nref, quadTrap)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setQuadTrap_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        procedure(real(RKG))                    :: getFunc
!        real(RKG)   , intent(in)                :: lb, ub
!        integer(IK) , intent(in)                :: nref
!        real(RKG)                               :: quadTrap
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    recursive impure module subroutine setQuadTrap_RK1(getFunc, lb, ub, nref, quadTrap)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setQuadTrap_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        procedure(real(RKG))                    :: getFunc
!        real(RKG)   , intent(in)                :: lb, ub
!        integer(IK) , intent(in)                :: nref
!        real(RKG)                               :: quadTrap
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    recursive impure module subroutine setQuadTrap_EP_NM_RK5(getFunc, lb, ub, nref, quadTrap, relerr, neval)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setQuadTrap_EP_NM_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        procedure(real(RKG))                    :: getFunc
!        real(RKG)   , intent(in)                :: lb, ub
!        integer(IK) , intent(in)                :: nref
!        integer(IK) , intent(out)               :: neval
!        real(RKG)   , intent(out)               :: relerr
!        real(RKG)                               :: quadTrap
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    recursive impure module subroutine setQuadTrap_EP_NM_RK4(getFunc, lb, ub, nref, quadTrap, relerr, neval)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setQuadTrap_EP_NM_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        procedure(real(RKG))                    :: getFunc
!        real(RKG)   , intent(in)                :: lb, ub
!        integer(IK) , intent(in)                :: nref
!        integer(IK) , intent(out)               :: neval
!        real(RKG)   , intent(out)               :: relerr
!        real(RKG)                               :: quadTrap
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    recursive impure module subroutine setQuadTrap_EP_NM_RK3(getFunc, lb, ub, nref, quadTrap, relerr, neval)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setQuadTrap_EP_NM_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        procedure(real(RKG))                    :: getFunc
!        real(RKG)   , intent(in)                :: lb, ub
!        integer(IK) , intent(in)                :: nref
!        integer(IK) , intent(out)               :: neval
!        real(RKG)   , intent(out)               :: relerr
!        real(RKG)                               :: quadTrap
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    recursive impure module subroutine setQuadTrap_EP_NM_RK2(getFunc, lb, ub, nref, quadTrap, relerr, neval)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setQuadTrap_EP_NM_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        procedure(real(RKG))                    :: getFunc
!        real(RKG)   , intent(in)                :: lb, ub
!        integer(IK) , intent(in)                :: nref
!        integer(IK) , intent(out)               :: neval
!        real(RKG)   , intent(out)               :: relerr
!        real(RKG)   , intent(inout)             :: quadTrap
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    recursive impure module subroutine setQuadTrap_EP_NM_RK1(getFunc, lb, ub, nref, quadTrap, relerr, neval)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setQuadTrap_EP_NM_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        procedure(real(RKG))                    :: getFunc
!        real(RKG)   , intent(in)                :: lb, ub
!        integer(IK) , intent(in)                :: nref
!        integer(IK) , intent(out)               :: neval
!        real(RKG)   , intent(out)               :: relerr
!        real(RKG)                               :: quadTrap
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    recursive impure setQuadTrap(getFunc, lb, ub, quadTrap, refinementStage, neval)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setQuadTrap
!#endif
!        implicit none
!        integer(IK), intent(in)     :: refinementStage
!        real(RK), intent(in)        :: lb,ub
!        real(RK), intent(inout)     :: quadTrap
!        integer(IK), intent(out)    :: neval
!        integer(IK)                 :: iFuncEval
!        real(RK)                    :: del,sum,tnm,x
!        procedure(integrand_proc)   :: getFunc
!        if (refinementStage == 1_IK) then
!            neval = 2_IK
!            quadTrap = 0.5_RK * (ub - lb) * (getFunc(lb) + getFunc(ub))
!        else
!            neval = 2**(refinementStage - 2_IK)
!            tnm = real(neval, RK)
!            del = (ub - lb) / tnm
!            x = lb + 0.5_RK * del
!            sum = 0._RK
!            do iFuncEval = 1_IK, neval
!                sum = sum + getFunc(x)
!                x = x + del
!            end do
!            quadTrap = 0.5_RK * (quadTrap + (ub-lb) * sum / tnm)
!        endif
!    end subroutine setQuadTrap

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!>  \cond excluded
#if 0
contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Return the refinement of the integration of an exponentially-decaying function on a semi-infinite.
    !>  This function is an integrator driver to be passed to [doQuadRombOpen](@ref doQuadRombOpen).
    !>
    !> \param[in]       getFunc         :   The input function to be integrated. It must have the interface specified by
    !!                                      [integrand_proc](@ref integrand_proc).
    !> \param[in]       lb              :   The lower limit of integration.
    !> \param[in]       ub              :   The upper limit of integration (typically set to `huge(1._RK)` to represent \f$+\infty\f$).
    !> \param[inout]    integral        :   The result of integration.
    !> \param[in]       refinementStage :   The number of refinements since the first call to the integrator.
    !> \param[out]      neval           :   The number of function evaluations made.
    !>
    !> \remark
    !> It is expected that `ub > lb > 0.0` must hold for this integrator to function properly.
    !>
    !> \remark
    !> Tested by Joshua Osborne on 5/28/2020 at 8:58 pm.
    recursive subroutine midexp(getFunc, lb, ub, integral, refinementStage, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: midexp
#endif
        use pm_except, only: LOGHUGE_RK
        implicit none
        integer(IK) , intent(in)    :: refinementStage
        real(RK)    , intent(in)    :: lb,ub
        integer(IK) , intent(out)   :: neval
        real(RK)    , intent(inout) :: integral
        procedure(integrand_proc)   :: getFunc
        real(RK)                    :: ddel,del,summ,x,lbTrans,ubTrans
        real(RK)                    :: inverseThreeNumFuncEval
        integer(IK)                 :: iFuncEval
        ubTrans = exp(-lb)
        lbTrans = 0._RK; if (ub<LOGHUGE_RK) lbTrans = exp(-ub)
        if (refinementStage==1_IK) then
            neval = 1_IK
            integral = (ubTrans-lbTrans)*getTransFunc(0.5_RK*(lbTrans+ubTrans))
        else
            neval = 3**(refinementStage-2)
            inverseThreeNumFuncEval = ONE_THIRD / neval
            del = (ubTrans-lbTrans) * inverseThreeNumFuncEval
            ddel = del + del
            x = lbTrans + 0.5_RK*del
            summ = 0._RK
            do iFuncEval = 1,neval
                summ = summ + getTransFunc(x)
                x = x + ddel
                summ = summ + getTransFunc(x)
                x = x + del
            end do
            integral = ONE_THIRD * integral + (ubTrans-lbTrans) * summ * inverseThreeNumFuncEval
            neval = 2_IK * neval
        end if
    contains
        !!>  \cond excluded
        function getTransFunc(x)
            real(RK), intent(in)    :: x
            real(RK)                :: getTransFunc
            getTransFunc = getFunc(-log(x)) / x
        end function getTransFunc
        !!>  \endcond excluded
    end subroutine midexp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This routine is an exact replacement for [midpnt](@ref midpnt), i.e., returns as `integral` the nth stage of refinement
    !> of the integral of funk from `lb` to `ub`, except that the function is evaluated at evenly spaced
    !> points in `1/x` rather than in `x`. This allows the upper limit `ub` to be as large and positive
    !> as the computer allows, or the lower limit `lb` to be as large and negative, but not both.
    !>
    !> \param[in]       getFunc         :   The input function to be integrated. It must have the interface specified by
    !!                                      [integrand_proc](@ref integrand_proc).
    !> \param[in]       lb        :   The lower limit of integration.
    !> \param[in]       ub        :   The upper limit of integration.
    !> \param[inout]    integral        :   The result of integration.
    !> \param[in]       refinementStage :   The number of refinements since the first call to the integrator.
    !> \param[out]      neval     :   The number of function evaluations made.
    !>
    !> \warning
    !> `lb * ub > 0.0` must hold for this integrator to function properly.
    !>
    !> \remark
    !> Tested by Joshua Osborne on 5/28/2020 at 8:58 pm.
    subroutine midinf(getFunc,lb,ub,integral,refinementStage,neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: midinf
#endif
        implicit none
        real(RK)    , intent(in)    :: lb,ub
        integer(IK) , intent(in)    :: refinementStage
        integer(IK) , intent(out)   :: neval
        real(RK)    , intent(inout) :: integral
        procedure(integrand_proc)   :: getFunc
        real(RK)                    :: lbTrans, ubTrans, del, ddel, summ, x
        real(RK)                    :: inverseThreeNumFuncEval
        integer(IK)                 :: iFuncEval
        ubTrans = 1.0_RK / lb
        lbTrans = 1.0_RK / ub
        if (refinementStage == 1_IK) then
            neval = 1_IK
            integral = (ubTrans-lbTrans) * getTransFunc(0.5_RK * (lbTrans+ubTrans))
        else
            neval = 3**(refinementStage-2)
            inverseThreeNumFuncEval = ONE_THIRD / neval
            del = (ubTrans-lbTrans) * inverseThreeNumFuncEval
            ddel = del + del
            x = lbTrans + 0.5_RK * del
            summ = 0._RK
            do iFuncEval = 1, neval
                summ = summ + getTransFunc(x)
                x = x + ddel
                summ = summ + getTransFunc(x)
                x = x + del
            end do
            integral = ONE_THIRD * integral + (ubTrans-lbTrans) * summ * inverseThreeNumFuncEval
            neval = 2_IK * neval
        end if
    contains
        !!>  \cond excluded
        function getTransFunc(x) result(transFunc)
            implicit none
            real(RK), intent(in)    :: x
            real(RK)                :: transFunc
            transFunc = getFunc(1._RK/x) / x**2
        end function getTransFunc
        !!>  \endcond excluded
    end subroutine midinf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This routine computes the nth stage of refinement of an extended midpoint rule.
    !> When called with `n = 1`, the routine returns as `integral` the crudest estimate of \f$\int_a^b f(x) ~ dx\f$.
    !> Subsequent calls with `n = 2, 3, ...` (in that sequential order) will improve the accuracy of `integral` by adding
    !> `(2/3) * 3n-1` additional interior points.
    !>
    !> \param[in]       getFunc         :   The input function to be integrated. It must have the interface specified by
    !!                                      [integrand_proc](@ref integrand_proc).
    !> \param[in]       lb        :   The lower limit of integration.
    !> \param[in]       ub        :   The upper limit of integration.
    !> \param[inout]    integral        :   The result of integration.
    !> \param[in]       refinementStage :   The number of refinements since the first call to the integrator.
    !> \param[out]      neval     :   The number of function evaluations made.
    !>
    !> \warning
    !> The argument `integral` should not be modified between sequential calls.
    !>
    !> \remark
    !> Tested by Joshua Osborne on 5/28/2020 at 8:55 pm.
    subroutine midpnt(getFunc,lb,ub,integral,refinementStage,neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: midpnt
#endif
        implicit none
        integer(IK) , intent(in)    :: refinementStage
        real(RK)    , intent(in)    :: lb, ub
        real(RK)    , intent(inout) :: integral
        integer(IK) , intent(out)   :: neval
        procedure(integrand_proc)   :: getFunc
        integer(IK)                 :: iFuncEval
        real(RK)                    :: ddel,del,summ,x
        real(RK)                    :: inverseThreeNumFuncEval
        if (refinementStage==1) then
            neval = 1_IK
            integral = (ub-lb) * getFunc( 0.5_RK * (lb+ub) )
        else
            neval = 3_IK**(refinementStage-2)
            inverseThreeNumFuncEval = ONE_THIRD / neval
            del = (ub-lb) * inverseThreeNumFuncEval
            ddel = del+del
            x = lb + 0.5_RK * del
            summ = 0._RK
            do iFuncEval = 1, neval
                summ = summ + getFunc(x)
                x = x + ddel
                summ = summ + getFunc(x)
                x = x + del
            end do
            integral = ONE_THIRD * integral + (ub-lb) * summ * inverseThreeNumFuncEval
            neval = 2_IK * neval
        end if
    end subroutine midpnt

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#endif
!>  \endcond excluded

end module pm_quadRomb