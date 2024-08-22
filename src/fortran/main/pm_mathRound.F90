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
!>  This module contains procedures and generic interfaces for rounding real values to whole **integers**.<br>
!>
!>  \test
!>  [test_pm_mathRound](@ref test_pm_mathRound)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Saturday July 20, 2024, 7:20 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_mathRound

    use pm_kind, only: SK, IK
    use pm_distUnif, only: xoshiro256ssw_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathRound"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **probabilistically-rounded** to the nearest integer value of the input real number.
    !>
    !>  \details
    !>  Unlike the intrinsic Fortran `nint()` function, the procedures of this generic
    !>  interface round the input real number to the nearest integer value **probabilistically**.<br>
    !>  The closer the input `val` is to its `nint(val)` result, the more likely it will be rounded to it.<br>
    !>  See below for example usage.<br>
    !>
    !>  \param[in]  val     :   The input scalar or array of the same rank as other array-like arguments, of<br>
    !>                          <ol>
    !>                              <li>    type `real` of kind \RKALL
    !>                          </ol>
    !>                          containing the number to be rounded probabilistically.<br>
    !>
    !>  \return
    !>  `whole`             :   The output scalar or array of the same rank and shape as the input arguments,
    !>                          of default `integer` kind \IK, containing the **probabilistically-rounded** input value to the nearest integer.<br>
    !>
    !>  \interface{pnint}
    !>  \code{.F90}
    !>
    !>      use pm_mathRound, only: pnint
    !>
    !>      whole = pnint(val)
    !>      whole = pnint(rng, val)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLog1p](@ref pm_mathLog1p::getLog1p)<br>
    !>  [get1mexp](@ref pm_math1mexp::get1mexp)<br>
    !>
    !>  \example{pnint}
    !>  \include{lineno} example/pm_mathRound/pnint/main.F90
    !>  \compilef{pnint}
    !>  \output{pnint}
    !>  \include{lineno} example/pm_mathRound/pnint/main.out.F90
    !>
    !>  \test
    !>  [test_pm_mathIntSqrt](@ref test_pm_mathIntSqrt)
    !>
    !>  \todo
    !>  \plow
    !>  This generic interface can be extended to output `integer` values of non-default kind.<br>
    !>
    !>  \todo
    !>  \plow
    !>  This generic interface can be extended to accept custom random number generators.<br>
    !>
    !>  \final{pnint}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface pnint

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function pnint_RNGD_IKD_RK5(val) result(whole)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: pnint_RNGD_IKD_RK5
#endif
        use pm_kind, only: IKG => IK, RKG => RK5
        real(RKG)               , intent(in)    :: val
        integer(IKG)                            :: whole
    end function
#endif

#if RK4_ENABLED
    impure elemental module function pnint_RNGD_IKD_RK4(val) result(whole)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: pnint_RNGD_IKD_RK4
#endif
        use pm_kind, only: IKG => IK, RKG => RK4
        real(RKG)               , intent(in)    :: val
        integer(IKG)                            :: whole
    end function
#endif

#if RK3_ENABLED
    impure elemental module function pnint_RNGD_IKD_RK3(val) result(whole)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: pnint_RNGD_IKD_RK3
#endif
        use pm_kind, only: IKG => IK, RKG => RK3
        real(RKG)               , intent(in)    :: val
        integer(IKG)                            :: whole
    end function
#endif

#if RK2_ENABLED
    impure elemental module function pnint_RNGD_IKD_RK2(val) result(whole)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: pnint_RNGD_IKD_RK2
#endif
        use pm_kind, only: IKG => IK, RKG => RK2
        real(RKG)               , intent(in)    :: val
        integer(IKG)                            :: whole
    end function
#endif

#if RK1_ENABLED
    impure elemental module function pnint_RNGD_IKD_RK1(val) result(whole)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: pnint_RNGD_IKD_RK1
#endif
        use pm_kind, only: IKG => IK, RKG => RK1
        real(RKG)               , intent(in)    :: val
        integer(IKG)                            :: whole
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function pnint_RNGX_IKD_RK5(rng, val) result(whole)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: pnint_RNGX_IKD_RK5
#endif
        use pm_kind, only: IKG => IK, RKG => RK5
        type(xoshiro256ssw_type), intent(inout) :: rng
        real(RKG)               , intent(in)    :: val
        integer(IKG)                            :: whole
    end function
#endif

#if RK4_ENABLED
    impure elemental module function pnint_RNGX_IKD_RK4(rng, val) result(whole)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: pnint_RNGX_IKD_RK4
#endif
        use pm_kind, only: IKG => IK, RKG => RK4
        type(xoshiro256ssw_type), intent(inout) :: rng
        real(RKG)               , intent(in)    :: val
        integer(IKG)                            :: whole
    end function
#endif

#if RK3_ENABLED
    impure elemental module function pnint_RNGX_IKD_RK3(rng, val) result(whole)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: pnint_RNGX_IKD_RK3
#endif
        use pm_kind, only: IKG => IK, RKG => RK3
        type(xoshiro256ssw_type), intent(inout) :: rng
        real(RKG)               , intent(in)    :: val
        integer(IKG)                            :: whole
    end function
#endif

#if RK2_ENABLED
    impure elemental module function pnint_RNGX_IKD_RK2(rng, val) result(whole)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: pnint_RNGX_IKD_RK2
#endif
        use pm_kind, only: IKG => IK, RKG => RK2
        type(xoshiro256ssw_type), intent(inout) :: rng
        real(RKG)               , intent(in)    :: val
        integer(IKG)                            :: whole
    end function
#endif

#if RK1_ENABLED
    impure elemental module function pnint_RNGX_IKD_RK1(rng, val) result(whole)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: pnint_RNGX_IKD_RK1
#endif
        use pm_kind, only: IKG => IK, RKG => RK1
        type(xoshiro256ssw_type), intent(inout) :: rng
        real(RKG)               , intent(in)    :: val
        integer(IKG)                            :: whole
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathRound