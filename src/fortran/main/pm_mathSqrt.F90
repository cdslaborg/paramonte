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
!>  This module contains procedures and generic interfaces for computing the **square root** of **integers**.
!>
!>  \details
!>  In number theory, the integer square root (**getSqrt**) of a non-negative integer \f$n\f$ is the non-negative integer \f$m\f$
!>  which is the greatest integer less than or equal to the square root of \f$n\f$,<br>
!>  \f{equation}{
!>      \mbox{getSqrt}(n) = \lfloor{\sqrt{n}}\rfloor ~.
!>  \f}
!>  For example,
!>  \f{equation}{
!>      \mbox{getSqrt}(27) = \lfloor{\sqrt{27}}\rfloor = \lfloor 5.19615242270663\ldots \rfloor = 5 ~.
!>  \f}
!>
!>  Two methods of finding the integer square root include linear and binary search,
!>  both of which are implemented by the procedures of this module.<br>
!>
!>  \note
!>  While the procedures of this module return `getSqrt(n)`, the greatest integer less than or equal to the square root of a given integer,
!>  one can readily compute the smallest integer `a` greater than or equal to the square root of `n` as `a = 1 + getSqrt(n - 1)`.<br>
!>  This number is equivalent to the ceiling of the exact square root of `n`.<br>
!>
!>  \benchmarks
!>
!>  \benchmark{getSqrt, The runtime performance of [getSqrt](@ref pm_mathSqrt::getSqrt) for `integer` vs. `real` input `ndim`.}
!>  \include{lineno} benchmark/pm_mathSqrt/getSqrt/main.F90
!>  \compilefb{getSqrt}
!>  \postprocb{getSqrt}
!>  \include{lineno} benchmark/pm_mathSqrt/getSqrt/main.py
!>  \visb{getSqrt}
!>  \image html benchmark/pm_mathSqrt/getSqrt/benchmark.getSqrt.runtime.png width=1000
!>  \image html benchmark/pm_mathSqrt/getSqrt/benchmark.getSqrt.runtime.ratio.png width=1000
!>  \moralb{getSqrt}
!>      -#  The benchmark procedures named `getSqrtBinary` and `getSqrtLinear` call the generic interface [getSqrt](@ref pm_mathSqrt::getSqrt)
!>          with a `method` argument of type [binary_type](@ref pm_search::binary_type) and [linear_type](@ref pm_search::linear_type) respectively.<br>
!>      -#  From the benchmark results it appears that the binary search algorithm for computing
!>          the integer square root is significantly faster than the linear search method.<br>
!>          The binary search algorithm is also faster than the naive approach 
!>          based on the formula `floor(sqrt(real(posint, RK)))`.<br>
!>
!>  \test
!>  [test_pm_mathSqrt](@ref test_pm_mathSqrt)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_mathSqrt

    use pm_kind, only: SK, IK
    use pm_search, only: linear, linear_type
    use pm_search, only: binary, binary_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathSqrt"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the integer square root of an input non-negative integer.
    !>
    !>  \details
    !>  See the documentation of [pm_mathSqrt](@ref pm_mathSqrt) for more information on integer square root.
    !>
    !>  \param[in]  posint  :   The input scalar of type `integer` of kind \IKALL
    !>                          containing the non-negative integer whose factoring is to be computed on return.
    !>
    !>  \return
    !>  `intSqrt`           :   The output `allocatable` vector of the same type and kind as the input `posint`,
    !>                          containing the **prime factors** of the input integer `posint`.<br>
    !>
    !>  \interface{getSqrt}
    !>  \code{.F90}
    !>
    !>      use pm_mathSqrt, only: getSqrt
    !>
    !>      sqrt = getSqrt(posint)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= posint` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warning
    !>  The input argument `posint` has the `value` attribute.<br>
    !>
    !>  \warning
    !>  The computation of the integer square root using the [linear search](@ref pm_search::linear_type) method
    !>  can become extremely costly for large input `posint` values, typically larger than `10**5`.<br>
    !>  For more information, see the [relevant benchmarks](@ref pm_mathSqrt).<br>
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  The procedures of this generic interface compute the integer square root method without undue overflow.<br>
    !>
    !>  \see
    !>  [getLog1p](@ref pm_mathLog1p::getLog1p)<br>
    !>  [get1mexp](@ref pm_math1mexp::get1mexp)<br>
    !>
    !>  \example{getSqrt}
    !>  \include{lineno} example/pm_mathSqrt/getSqrt/main.F90
    !>  \compilef{getSqrt}
    !>  \output{getSqrt}
    !>  \include{lineno} example/pm_mathSqrt/getSqrt/main.out.F90
    !>
    !>  \test
    !>  [test_pm_mathIntSqrt](@ref test_pm_mathIntSqrt)
    !>
    !>  \todo
    !>  \plow
    !>  A binary-representation calculation of the integer square root should be added in the future.<br>
    !>
    !>  \final{getSqrt}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface getSqrt

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE elemental module function getSqrtDef_IK5(posint) result(intSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSqrtDef_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)        , value         :: posint
        integer(IKG)                        :: intSqrt
    end function
#endif

#if IK4_ENABLED
    PURE elemental module function getSqrtDef_IK4(posint) result(intSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSqrtDef_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)        , value         :: posint
        integer(IKG)                        :: intSqrt
    end function
#endif

#if IK3_ENABLED
    PURE elemental module function getSqrtDef_IK3(posint) result(intSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSqrtDef_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)        , value         :: posint
        integer(IKG)                        :: intSqrt
    end function
#endif

#if IK2_ENABLED
    PURE elemental module function getSqrtDef_IK2(posint) result(intSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSqrtDef_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)        , value         :: posint
        integer(IKG)                        :: intSqrt
    end function
#endif

#if IK1_ENABLED
    PURE elemental module function getSqrtDef_IK1(posint) result(intSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSqrtDef_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)        , value         :: posint
        integer(IKG)                        :: intSqrt
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE elemental module function getSqrtBin_IK5(posint, method) result(intSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSqrtBin_IK5
#endif
        use pm_kind, only: IKG => IK5
        type(binary_type)   , intent(in)    :: method
        integer(IKG)        , value         :: posint
        integer(IKG)                        :: intSqrt
    end function
#endif

#if IK4_ENABLED
    PURE elemental module function getSqrtBin_IK4(posint, method) result(intSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSqrtBin_IK4
#endif
        use pm_kind, only: IKG => IK4
        type(binary_type)   , intent(in)    :: method
        integer(IKG)        , value         :: posint
        integer(IKG)                        :: intSqrt
    end function
#endif

#if IK3_ENABLED
    PURE elemental module function getSqrtBin_IK3(posint, method) result(intSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSqrtBin_IK3
#endif
        use pm_kind, only: IKG => IK3
        type(binary_type)   , intent(in)    :: method
        integer(IKG)        , value         :: posint
        integer(IKG)                        :: intSqrt
    end function
#endif

#if IK2_ENABLED
    PURE elemental module function getSqrtBin_IK2(posint, method) result(intSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSqrtBin_IK2
#endif
        use pm_kind, only: IKG => IK2
        type(binary_type)   , intent(in)    :: method
        integer(IKG)        , value         :: posint
        integer(IKG)                        :: intSqrt
    end function
#endif

#if IK1_ENABLED
    PURE elemental module function getSqrtBin_IK1(posint, method) result(intSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSqrtBin_IK1
#endif
        use pm_kind, only: IKG => IK1
        type(binary_type)   , intent(in)    :: method
        integer(IKG)        , value         :: posint
        integer(IKG)                        :: intSqrt
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE elemental module function getSqrtLin_IK5(posint, method) result(intSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSqrtLin_IK5
#endif
        use pm_kind, only: IKG => IK5
        type(linear_type)   , intent(in)    :: method
        integer(IKG)        , value         :: posint
        integer(IKG)                        :: intSqrt
    end function
#endif

#if IK4_ENABLED
    PURE elemental module function getSqrtLin_IK4(posint, method) result(intSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSqrtLin_IK4
#endif
        use pm_kind, only: IKG => IK4
        type(linear_type)   , intent(in)    :: method
        integer(IKG)        , value         :: posint
        integer(IKG)                        :: intSqrt
    end function
#endif

#if IK3_ENABLED
    PURE elemental module function getSqrtLin_IK3(posint, method) result(intSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSqrtLin_IK3
#endif
        use pm_kind, only: IKG => IK3
        type(linear_type)   , intent(in)    :: method
        integer(IKG)        , value         :: posint
        integer(IKG)                        :: intSqrt
    end function
#endif

#if IK2_ENABLED
    PURE elemental module function getSqrtLin_IK2(posint, method) result(intSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSqrtLin_IK2
#endif
        use pm_kind, only: IKG => IK2
        type(linear_type)   , intent(in)    :: method
        integer(IKG)        , value         :: posint
        integer(IKG)                        :: intSqrt
    end function
#endif

#if IK1_ENABLED
    PURE elemental module function getSqrtLin_IK1(posint, method) result(intSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSqrtLin_IK1
#endif
        use pm_kind, only: IKG => IK1
        type(linear_type)   , intent(in)    :: method
        integer(IKG)        , value         :: posint
        integer(IKG)                        :: intSqrt
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathSqrt
