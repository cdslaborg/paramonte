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
!>  This module contains procedures, generic interfaces, and types for generating default values for optional arguments.
!>
!>
!>
!>  \test
!>  [test_pm_arrayInsert](@ref test_pm_arrayInsert)
!>
!>  \final
!>  This module is inspired by a similar functionality in the Fortran stdlib.<br>
!>  If you use or redistribute this module you should also acknowledge and cite the [Fortran stdlib](https://github.com/fortran-lang/stdlib).<br>
!>
!>  \author
!>  \AmirShahmoradi, Saturday 10:17 PM, December 5, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_option

    use pm_kind, only: SK, IK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_option"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the value of the `optional` input argument if it is present, otherwise, return the input `default` value.
    !>
    !>  \param[in]  default     :   The input scalar, are array of arbitrary rank, of either<br>
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL, or <br>
    !>                                  <li>    type `logical` of kind \LKALL, or <br>
    !>                                  <li>    type `integer` of kind \IKALL, or <br>
    !>                                  <li>    type `complex` of kind \CKALL, or <br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ol>
    !>                              or,
    !>                              <ol>
    !>                                  <li>    a scalar assumed-length `character` of kind \SKALL, <br>
    !>                              </ol>
    !>                              whose value will be returned as the `output` if `optional` is missing among the input arguments.<br>
    !>  \param[in]  optional    :   The input value of the same type, kind, and rank as the input `default` argument,
    !>                              representing the `optional` input argument whose value will be returned if present.<br>
    !>                              (**optional**, default = `default`)
    !>
    !>  \return
    !>  `value`                 :   The output value of the same type, kind, and rank as the input `default`.<br>
    !>                              <ol>
    !>                                  <li>    If `optional` is present as an input argument, `value` takes the value of `optional`.<br>
    !>                                  <li>    If `optional` is missing as an input argument, `value` takes the value of `default`.<br>
    !>                              </ol>
    !>                              If `value` is of type `character`, then its length type parameter is `len(default)`.<br>
    !>
    !>  \interface{getOption}
    !>  \code{.F90}
    !>
    !>      use pm_option, only: getOption
    !>
    !>      defaultValue = getOption(default, optional)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `len(optional) <= len(default)` must hold when the input arguments are of type `character`.<br>
    !>  This is to ensure proper full assignment of the `optional` value to the output whose length is that of `default`.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  The functions under this generic interface provide are solely intended to a provide convenient shortcut
    !>  to evaluating whether an optional argument is present or not, and if missing, to return a default value.<br>
    !>  Due to significant performance penalty associated with the use of these procedures, their usage should be
    !>  limited to non-performance-critical sections of code.<br>
    !>
    !>  \example{getOption}
    !>  \include{lineno} example/pm_option/getOption/main.F90
    !>  \compilef{getOption}
    !>  \output{getOption}
    !>  \include{lineno} example/pm_option/getOption/main.out.F90
    !>
    !>  \benchmarks
    !>
    !>  \benchmark{getOption_vs_direct, The runtime performance of [getOption](@ref pm_option::getOption) vs. direct optional choice}
    !>  \include{lineno} benchmark/pm_option/getOption_vs_direct/main.F90
    !>  \compilefb{getOption_vs_direct}
    !>  \postprocb{getOption_vs_direct}
    !>  \include{lineno} benchmark/pm_option/getOption_vs_direct/main.py
    !>  \visb{getOption_vs_direct}
    !>  \image html benchmark/pm_option/getOption_vs_direct/benchmark.getOption_vs_direct.runtime.png width=1000
    !>  \image html benchmark/pm_option/getOption_vs_direct/benchmark.getOption_vs_direct.runtime.ratio.png width=1000
    !>  \moralb{getOption_vs_direct}
    !>      -#  The benchmark results indicate that explicit if-blocks for optional arguments tends to be 100-1000 times
    !>          faster than calling the convenience functions under the generic interface [getOption](@ref pm_option::getOption).<br>
    !>          This performance difference tends to be about a factor of 10 times for scalar `optional` arguments and grows
    !>          substantially to larger factors with switching to increasing-size array-like `optional` dummy arguments.<br>
    !>      -#  Despite significant performance degradation with the use of [getOption](@ref pm_option::getOption), note that
    !>          the overall average overhead of calling [getOption](@ref pm_option::getOption) remains extremely small as of 2022 AD,
    !>          on the order of nano to micro seconds depending on the size of the optional argument from scalar to large arrays.<br>
    !>      -#  Additionally, the significance of performance degradation, if any, depends entirely on the ability of the compiler to inline the procedure.<br>
    !>          The observed performance degradations report here occur without any compiler inlining attempts.<br>
    !>      -#  Consequently, **there is practically no harm in using [getOption](@ref pm_option::getOption) in non-performance-critical
    !>          sections of a codebase**, that is, **parts of a codebase that are to be called less than billions of times at runtime**.<br>
    !>
    !>  \test
    !>  [test_pm_option](@ref test_pm_option)
    !>
    !>  \final{getOption}
    !>  The functions under the generic interface [getOption](@ref pm_option::getOption) are inspired by and further extend the
    !>  functionalities implemented in [optval](https://github.com/fortran-lang/stdlib/blob/stdlib-fpm/src/stdlib_optval.F90)
    !>  of the [Fortran stdlib](https://github.com/fortran-lang/stdlib).
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2020, 01:12 AM, Dallas, Texas
    interface getOption

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE elemental module function getOption_D0_SK5(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)                , intent(in)            :: default
        character(*,SKG)                , intent(in), optional  :: optional
        character(len(default, IK),SKG)                         :: value
    end function
#endif

#if SK4_ENABLED
    PURE elemental module function getOption_D0_SK4(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)                , intent(in)            :: default
        character(*,SKG)                , intent(in), optional  :: optional
        character(len(default, IK),SKG)                         :: value
    end function
#endif

#if SK3_ENABLED
    PURE elemental module function getOption_D0_SK3(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)                , intent(in)            :: default
        character(*,SKG)                , intent(in), optional  :: optional
        character(len(default, IK),SKG)                         :: value
    end function
#endif

#if SK2_ENABLED
    PURE elemental module function getOption_D0_SK2(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)                , intent(in)            :: default
        character(*,SKG)                , intent(in), optional  :: optional
        character(len(default, IK),SKG)                         :: value
    end function
#endif

#if SK1_ENABLED
    PURE elemental module function getOption_D0_SK1(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)                , intent(in)            :: default
        character(*,SKG)                , intent(in), optional  :: optional
        character(len(default, IK),SKG)                         :: value
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure elemental module function getOption_D0_IK5(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                    , intent(in)            :: default
        integer(IKG)                    , intent(in), optional  :: optional
        integer(IKG)                                            :: value
    end function
#endif

#if IK4_ENABLED
    pure elemental module function getOption_D0_IK4(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                    , intent(in)            :: default
        integer(IKG)                    , intent(in), optional  :: optional
        integer(IKG)                                            :: value
    end function
#endif

#if IK3_ENABLED
    pure elemental module function getOption_D0_IK3(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                    , intent(in)            :: default
        integer(IKG)                    , intent(in), optional  :: optional
        integer(IKG)                                            :: value
    end function
#endif

#if IK2_ENABLED
    pure elemental module function getOption_D0_IK2(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                    , intent(in)            :: default
        integer(IKG)                    , intent(in), optional  :: optional
        integer(IKG)                                            :: value
    end function
#endif

#if IK1_ENABLED
    pure elemental module function getOption_D0_IK1(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                    , intent(in)            :: default
        integer(IKG)                    , intent(in), optional  :: optional
        integer(IKG)                                            :: value
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure elemental module function getOption_D0_LK5(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                    , intent(in)            :: default
        logical(LKG)                    , intent(in), optional  :: optional
        logical(LKG)                                            :: value
    end function
#endif

#if LK4_ENABLED
    pure elemental module function getOption_D0_LK4(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                    , intent(in)            :: default
        logical(LKG)                    , intent(in), optional  :: optional
        logical(LKG)                                            :: value
    end function
#endif

#if LK3_ENABLED
    pure elemental module function getOption_D0_LK3(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                    , intent(in)            :: default
        logical(LKG)                    , intent(in), optional  :: optional
        logical(LKG)                                            :: value
    end function
#endif

#if LK2_ENABLED
    pure elemental module function getOption_D0_LK2(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                    , intent(in)            :: default
        logical(LKG)                    , intent(in), optional  :: optional
        logical(LKG)                                            :: value
    end function
#endif

#if LK1_ENABLED
    pure elemental module function getOption_D0_LK1(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                    , intent(in)            :: default
        logical(LKG)                    , intent(in), optional  :: optional
        logical(LKG)                                            :: value
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function getOption_D0_CK5(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                    , intent(in)            :: default
        complex(CKG)                    , intent(in), optional  :: optional
        complex(CKG)                                            :: value
    end function
#endif

#if CK4_ENABLED
    pure elemental module function getOption_D0_CK4(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                    , intent(in)            :: default
        complex(CKG)                    , intent(in), optional  :: optional
        complex(CKG)                                            :: value
    end function
#endif

#if CK3_ENABLED
    pure elemental module function getOption_D0_CK3(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                    , intent(in)            :: default
        complex(CKG)                    , intent(in), optional  :: optional
        complex(CKG)                                            :: value
    end function
#endif

#if CK2_ENABLED
    pure elemental module function getOption_D0_CK2(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                    , intent(in)            :: default
        complex(CKG)                    , intent(in), optional  :: optional
        complex(CKG)                                            :: value
    end function
#endif

#if CK1_ENABLED
    pure elemental module function getOption_D0_CK1(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                    , intent(in)            :: default
        complex(CKG)                    , intent(in), optional  :: optional
        complex(CKG)                                            :: value
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function getOption_D0_RK5(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                       , intent(in)            :: default
        real(RKG)                       , intent(in), optional  :: optional
        real(RKG)                                               :: value
    end function
#endif

#if RK4_ENABLED
    pure elemental module function getOption_D0_RK4(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                       , intent(in)            :: default
        real(RKG)                       , intent(in), optional  :: optional
        real(RKG)                                               :: value
    end function
#endif

#if RK3_ENABLED
    pure elemental module function getOption_D0_RK3(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                       , intent(in)            :: default
        real(RKG)                       , intent(in), optional  :: optional
        real(RKG)                                               :: value
    end function
#endif

#if RK2_ENABLED
    pure elemental module function getOption_D0_RK2(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                       , intent(in)            :: default
        real(RKG)                       , intent(in), optional  :: optional
        real(RKG)                                               :: value
    end function
#endif

#if RK1_ENABLED
    pure elemental module function getOption_D0_RK1(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                       , intent(in)            :: default
        real(RKG)                       , intent(in), optional  :: optional
        real(RKG)                                               :: value
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getOption_D1_SK5(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)                , intent(in)            :: default(:)
        character(*,SKG)                , intent(in), optional  :: optional(:)
        character(len(default),SKG)                             :: value(size(default, 1, IK))
    end function
#endif

#if SK4_ENABLED
    PURE module function getOption_D1_SK4(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)                , intent(in)            :: default(:)
        character(*,SKG)                , intent(in), optional  :: optional(:)
        character(len(default),SKG)                             :: value(size(default, 1, IK))
    end function
#endif

#if SK3_ENABLED
    PURE module function getOption_D1_SK3(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)                , intent(in)            :: default(:)
        character(*,SKG)                , intent(in), optional  :: optional(:)
        character(len(default),SKG)                             :: value(size(default, 1, IK))
    end function
#endif

#if SK2_ENABLED
    PURE module function getOption_D1_SK2(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)                , intent(in)            :: default(:)
        character(*,SKG)                , intent(in), optional  :: optional(:)
        character(len(default),SKG)                             :: value(size(default, 1, IK))
    end function
#endif

#if SK1_ENABLED
    PURE module function getOption_D1_SK1(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)                , intent(in)            :: default(:)
        character(*,SKG)                , intent(in), optional  :: optional(:)
        character(len(default),SKG)                             :: value(size(default, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module function getOption_D1_IK5(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                    , intent(in)            :: default(:)
        integer(IKG)                    , intent(in), optional  :: optional(:)
        integer(IKG)                                            :: value(size(default, 1, IK))
    end function
#endif

#if IK4_ENABLED
    pure module function getOption_D1_IK4(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                    , intent(in)            :: default(:)
        integer(IKG)                    , intent(in), optional  :: optional(:)
        integer(IKG)                                            :: value(size(default, 1, IK))
    end function
#endif

#if IK3_ENABLED
    pure module function getOption_D1_IK3(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                    , intent(in)            :: default(:)
        integer(IKG)                    , intent(in), optional  :: optional(:)
        integer(IKG)                                            :: value(size(default, 1, IK))
    end function
#endif

#if IK2_ENABLED
    pure module function getOption_D1_IK2(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                    , intent(in)            :: default(:)
        integer(IKG)                    , intent(in), optional  :: optional(:)
        integer(IKG)                                            :: value(size(default, 1, IK))
    end function
#endif

#if IK1_ENABLED
    pure module function getOption_D1_IK1(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                    , intent(in)            :: default(:)
        integer(IKG)                    , intent(in), optional  :: optional(:)
        integer(IKG)                                            :: value(size(default, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module function getOption_D1_LK5(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                    , intent(in)            :: default(:)
        logical(LKG)                    , intent(in), optional  :: optional(:)
        logical(LKG)                                            :: value(size(default, 1, IK))
    end function
#endif

#if LK4_ENABLED
    pure module function getOption_D1_LK4(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                    , intent(in)            :: default(:)
        logical(LKG)                    , intent(in), optional  :: optional(:)
        logical(LKG)                                            :: value(size(default, 1, IK))
    end function
#endif

#if LK3_ENABLED
    pure module function getOption_D1_LK3(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                    , intent(in)            :: default(:)
        logical(LKG)                    , intent(in), optional  :: optional(:)
        logical(LKG)                                            :: value(size(default, 1, IK))
    end function
#endif

#if LK2_ENABLED
    pure module function getOption_D1_LK2(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                    , intent(in)            :: default(:)
        logical(LKG)                    , intent(in), optional  :: optional(:)
        logical(LKG)                                            :: value(size(default, 1, IK))
    end function
#endif

#if LK1_ENABLED
    pure module function getOption_D1_LK1(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                    , intent(in)            :: default(:)
        logical(LKG)                    , intent(in), optional  :: optional(:)
        logical(LKG)                                            :: value(size(default, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function getOption_D1_CK5(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                    , intent(in)            :: default(:)
        complex(CKG)                    , intent(in), optional  :: optional(:)
        complex(CKG)                                            :: value(size(default, 1, IK))
    end function
#endif

#if CK4_ENABLED
    pure module function getOption_D1_CK4(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                    , intent(in)            :: default(:)
        complex(CKG)                    , intent(in), optional  :: optional(:)
        complex(CKG)                                            :: value(size(default, 1, IK))
    end function
#endif

#if CK3_ENABLED
    pure module function getOption_D1_CK3(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                    , intent(in)            :: default(:)
        complex(CKG)                    , intent(in), optional  :: optional(:)
        complex(CKG)                                            :: value(size(default, 1, IK))
    end function
#endif

#if CK2_ENABLED
    pure module function getOption_D1_CK2(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                    , intent(in)            :: default(:)
        complex(CKG)                    , intent(in), optional  :: optional(:)
        complex(CKG)                                            :: value(size(default, 1, IK))
    end function
#endif

#if CK1_ENABLED
    pure module function getOption_D1_CK1(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                    , intent(in)            :: default(:)
        complex(CKG)                    , intent(in), optional  :: optional(:)
        complex(CKG)                                            :: value(size(default, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module function getOption_D1_RK5(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                       , intent(in)            :: default(:)
        real(RKG)                       , intent(in), optional  :: optional(:)
        real(RKG)                                               :: value(size(default, 1, IK))
    end function
#endif

#if RK4_ENABLED
    pure module function getOption_D1_RK4(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                       , intent(in)            :: default(:)
        real(RKG)                       , intent(in), optional  :: optional(:)
        real(RKG)                                               :: value(size(default, 1, IK))
    end function
#endif

#if RK3_ENABLED
    pure module function getOption_D1_RK3(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                       , intent(in)            :: default(:)
        real(RKG)                       , intent(in), optional  :: optional(:)
        real(RKG)                                               :: value(size(default, 1, IK))
    end function
#endif

#if RK2_ENABLED
    pure module function getOption_D1_RK2(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                       , intent(in)            :: default(:)
        real(RKG)                       , intent(in), optional  :: optional(:)
        real(RKG)                                               :: value(size(default, 1, IK))
    end function
#endif

#if RK1_ENABLED
    pure module function getOption_D1_RK1(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                       , intent(in)            :: default(:)
        real(RKG)                       , intent(in), optional  :: optional(:)
        real(RKG)                                               :: value(size(default, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getOption_D2_SK5(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)                , intent(in)            :: default(:,:)
        character(*,SKG)                , intent(in), optional  :: optional(:,:)
        character(len(default),SKG)                             :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

#if SK4_ENABLED
    PURE module function getOption_D2_SK4(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)                , intent(in)            :: default(:,:)
        character(*,SKG)                , intent(in), optional  :: optional(:,:)
        character(len(default),SKG)                             :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

#if SK3_ENABLED
    PURE module function getOption_D2_SK3(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)                , intent(in)            :: default(:,:)
        character(*,SKG)                , intent(in), optional  :: optional(:,:)
        character(len(default),SKG)                             :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

#if SK2_ENABLED
    PURE module function getOption_D2_SK2(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)                , intent(in)            :: default(:,:)
        character(*,SKG)                , intent(in), optional  :: optional(:,:)
        character(len(default),SKG)                             :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

#if SK1_ENABLED
    PURE module function getOption_D2_SK1(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)                , intent(in)            :: default(:,:)
        character(*,SKG)                , intent(in), optional  :: optional(:,:)
        character(len(default),SKG)                             :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module function getOption_D2_IK5(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                    , intent(in)            :: default(:,:)
        integer(IKG)                    , intent(in), optional  :: optional(:,:)
        integer(IKG)                                            :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

#if IK4_ENABLED
    pure module function getOption_D2_IK4(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                    , intent(in)            :: default(:,:)
        integer(IKG)                    , intent(in), optional  :: optional(:,:)
        integer(IKG)                                            :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

#if IK3_ENABLED
    pure module function getOption_D2_IK3(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                    , intent(in)            :: default(:,:)
        integer(IKG)                    , intent(in), optional  :: optional(:,:)
        integer(IKG)                                            :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

#if IK2_ENABLED
    pure module function getOption_D2_IK2(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                    , intent(in)            :: default(:,:)
        integer(IKG)                    , intent(in), optional  :: optional(:,:)
        integer(IKG)                                            :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

#if IK1_ENABLED
    pure module function getOption_D2_IK1(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                    , intent(in)            :: default(:,:)
        integer(IKG)                    , intent(in), optional  :: optional(:,:)
        integer(IKG)                                            :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module function getOption_D2_LK5(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                    , intent(in)            :: default(:,:)
        logical(LKG)                    , intent(in), optional  :: optional(:,:)
        logical(LKG)                                            :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

#if LK4_ENABLED
    pure module function getOption_D2_LK4(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                    , intent(in)            :: default(:,:)
        logical(LKG)                    , intent(in), optional  :: optional(:,:)
        logical(LKG)                                            :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

#if LK3_ENABLED
    pure module function getOption_D2_LK3(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                    , intent(in)            :: default(:,:)
        logical(LKG)                    , intent(in), optional  :: optional(:,:)
        logical(LKG)                                            :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

#if LK2_ENABLED
    pure module function getOption_D2_LK2(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                    , intent(in)            :: default(:,:)
        logical(LKG)                    , intent(in), optional  :: optional(:,:)
        logical(LKG)                                            :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

#if LK1_ENABLED
    pure module function getOption_D2_LK1(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                    , intent(in)            :: default(:,:)
        logical(LKG)                    , intent(in), optional  :: optional(:,:)
        logical(LKG)                                            :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function getOption_D2_CK5(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                    , intent(in)            :: default(:,:)
        complex(CKG)                    , intent(in), optional  :: optional(:,:)
        complex(CKG)                                            :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

#if CK4_ENABLED
    pure module function getOption_D2_CK4(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                    , intent(in)            :: default(:,:)
        complex(CKG)                    , intent(in), optional  :: optional(:,:)
        complex(CKG)                                            :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

#if CK3_ENABLED
    pure module function getOption_D2_CK3(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                    , intent(in)            :: default(:,:)
        complex(CKG)                    , intent(in), optional  :: optional(:,:)
        complex(CKG)                                            :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

#if CK2_ENABLED
    pure module function getOption_D2_CK2(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                    , intent(in)            :: default(:,:)
        complex(CKG)                    , intent(in), optional  :: optional(:,:)
        complex(CKG)                                            :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

#if CK1_ENABLED
    pure module function getOption_D2_CK1(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                    , intent(in)            :: default(:,:)
        complex(CKG)                    , intent(in), optional  :: optional(:,:)
        complex(CKG)                                            :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module function getOption_D2_RK5(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                       , intent(in)            :: default(:,:)
        real(RKG)                       , intent(in), optional  :: optional(:,:)
        real(RKG)                                               :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

#if RK4_ENABLED
    pure module function getOption_D2_RK4(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                       , intent(in)            :: default(:,:)
        real(RKG)                       , intent(in), optional  :: optional(:,:)
        real(RKG)                                               :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

#if RK3_ENABLED
    pure module function getOption_D2_RK3(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                       , intent(in)            :: default(:,:)
        real(RKG)                       , intent(in), optional  :: optional(:,:)
        real(RKG)                                               :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

#if RK2_ENABLED
    pure module function getOption_D2_RK2(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                       , intent(in)            :: default(:,:)
        real(RKG)                       , intent(in), optional  :: optional(:,:)
        real(RKG)                                               :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

#if RK1_ENABLED
    pure module function getOption_D2_RK1(default, optional) result(value)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getOption_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                       , intent(in)            :: default(:,:)
        real(RKG)                       , intent(in), optional  :: optional(:,:)
        real(RKG)                                               :: value(size(default, 1, IK), size(default,2,IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_option