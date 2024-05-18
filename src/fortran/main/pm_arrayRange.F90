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
!>  This module contains procedures and generic interfaces for generating ranges
!>  of discrete `character`, `integer`, or `real` -valued sequences with minimum-possible or user-specified fixed linear spacings.<br>
!>
!>  \see
!>  [pm_arrayChange](@ref pm_arrayChange)<br>
!>  [pm_arraySpace](@ref pm_arraySpace)<br>
!>  [pm_arraySpace](@ref pm_arraySpace)<br>
!>
!>  \test
!>  [test_pm_arraySpace](@ref test_pm_arraySpace)<br>
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Tuesday 08:49 PM, August 10, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayRange

    use pm_kind, only: SK, IK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_arrayRange"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate minimally-spaced `character`, `integer`, `real` sequences or sequences at fixed intervals of size `step` from `start` to `stop`.<br>
    !>
    !>  \param[in]  start   :   The input scalar of the same type and kind as the output `range`,
    !>                          representing the start point of the output `range`.<br>
    !>                          If `start` is of type `character`, then it must be of length type parameter `1`.<br>
    !>  \param[in]  stop    :   The input scalar of the same type and kind as the output `range` representing the end point of the output `range`.<br>
    !>                          If `stop` is of type `character`, then it must be of length type parameter `1`.<br>
    !>  \param[in]  step    :   The **non-zero** input scalar representing the size of the interval between adjacent values in the output `range`.<br>
    !>                          <ol>
    !>                              <li>    If `range` is of type `character`, then `step` must have the default `integer` kind \IK.<br>
    !>                              <li>    If `range` is of type `integer`, or `real` then `step` must have the same type and kind parameter as `range`.<br>
    !>                          </ol>
    !>                          (**optional**, default = dynamically set to the smallest possible spacing between any two values between `start` and `stop`.<br>
    !>                          This means `step = sign(1, stop - start)` for output `range` of type `character` or `integer` and a variable `step`
    !>                          set by the intrinsic `nearest(range(i), stop - start)` for the element `i` of range.)
    !>
    !>  \return
    !>  `range`             :   The output `allocatable` of
    !>                          <ol>
    !>                              <li>    type `character` of kind \SKALL of length type parameter `1`,<br>
    !>                              <li>    type `integer` of kind \IKALL,<br>
    !>                              <li>    type `real` of kind \RKALL,<br>
    !>                          </ol>
    !>                          containing the minimally or evenly -spaced sequence within the interval starting with `start` and `stop`.<br>
    !>                          Note that the specified value of `stop` may not necessarily be the ending value in the sequence.<br>
    !>                          <ol>
    !>                              <li>    If `range` is of type `character`, then `range` has the shape `(1 : 1 + (ichar(stop) - ichar(start)) / step)`.<br>
    !>                              <li>    If `range` is of type `integer`, then `range` has the shape `(1 : 1 + (stop - start) / step)`.<br>
    !>                              <li>    If `range` is of type `real` and `step` is missing, then `range` is an `allocatable` of shape `(:)`
    !>                                      whose size is determined at runtime based on the `real` number representation model of the computer.<br>
    !>                              <li>    If `range` is of type `real`, then `range` has the shape `(1 : 1 + (stop - start) / step)`.<br>
    !>                          </ol>
    !>
    !>  \interface{getRange}
    !>  \code{.F90}
    !>
    !>      use pm_arrayRange, only: getRange
    !>
    !>      ! real range.
    !>
    !>      range = getRange(start, stop) ! allocatable (unknown size)
    !>      range(1 : 1 + (stop - start) / step) = getRange(start, stop, step)
    !>
    !>      ! integer range.
    !>
    !>      range(1 : 1 +  stop - start) = getRange(start, stop)
    !>      range(1 : 1 + (stop - start) / step) = getRange(start, stop, step)
    !>
    !>      ! character range.
    !>
    !>      range(1 : 1 +  ichar(stop) - ichar(start)) = getRange(start, stop)
    !>      range(1 : 1 + (ichar(stop) - ichar(start)) / step) = getRange(start, stop, step)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `step /= 0` must hold for the corresponding input arguments.<br>
    !>  The condition `len(stop) == 1` must hold for the corresponding input `character` arguments.<br>
    !>  The condition `len(start) == 1` must hold for the corresponding input `character` arguments.<br>
    !>  \vericons
    !>
    !>  \warning
    !>  Beware that the output `real` sequences with default minimum spacings with arbitrary `start` and `stop` can readily overflow the computer memory.<br>
    !>  This is because the number of representable `real` values between any two given `start` and `stop` is generally extremely large.<br>
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getRange](@ref pm_arrayRange::getRange)<br>
    !>  [setRange](@ref pm_arrayRange::setRange)<br>
    !>  [getLinSpace](@ref pm_arraySpace::getLinSpace)<br>
    !>  [setLinSpace](@ref pm_arraySpace::setLinSpace)<br>
    !>  [getLogSpace](@ref pm_arraySpace::getLogSpace)<br>
    !>  [setLogSpace](@ref pm_arraySpace::setLogSpace)<br>
    !>
    !>  \example{getRange}
    !>  \include{lineno} example/pm_arrayRange/getRange/main.F90
    !>  \compilef{getRange}
    !>  \output{getRange}
    !>  \include{lineno} example/pm_arrayRange/getRange/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayRange](@ref test_pm_arrayRange)
    !>
    !>  \final{getRange}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 08:49 PM, August 10, 2021, Dallas, TX
    interface getRange

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRangeUnit_D0_SK5(start, stop) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeUnit_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(1,SKC)        , intent(in)                                        :: start, stop
        character(abs(ichar(stop) - ichar(start)) + 1,SKC)                          :: range
    end function
#endif

#if SK4_ENABLED
    PURE module function getRangeUnit_D0_SK4(start, stop) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeUnit_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(1,SKC)        , intent(in)                                        :: start, stop
        character(abs(ichar(stop) - ichar(start)) + 1,SKC)                          :: range
    end function
#endif

#if SK3_ENABLED
    PURE module function getRangeUnit_D0_SK3(start, stop) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeUnit_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(1,SKC)        , intent(in)                                        :: start, stop
        character(abs(ichar(stop) - ichar(start)) + 1,SKC)                          :: range
    end function
#endif

#if SK2_ENABLED
    PURE module function getRangeUnit_D0_SK2(start, stop) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeUnit_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(1,SKC)        , intent(in)                                        :: start, stop
        character(abs(ichar(stop) - ichar(start)) + 1,SKC)                          :: range
    end function
#endif

#if SK1_ENABLED
    PURE module function getRangeUnit_D0_SK1(start, stop) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeUnit_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(1,SKC)        , intent(in)                                        :: start, stop
        character(abs(ichar(stop) - ichar(start)) + 1,SKC)                          :: range
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getRangeUnit_D1_IK5(start, stop) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeUnit_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)                                        :: start, stop
        integer(IKC)                                                                :: range(abs(stop - start) + 1_IKC)
    end function
#endif

#if IK4_ENABLED
    PURE module function getRangeUnit_D1_IK4(start, stop) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeUnit_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)                                        :: start, stop
        integer(IKC)                                                                :: range(abs(stop - start) + 1_IKC)
    end function
#endif

#if IK3_ENABLED
    PURE module function getRangeUnit_D1_IK3(start, stop) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeUnit_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)                                        :: start, stop
        integer(IKC)                                                                :: range(abs(stop - start) + 1_IKC)
    end function
#endif

#if IK2_ENABLED
    PURE module function getRangeUnit_D1_IK2(start, stop) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeUnit_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)                                        :: start, stop
        integer(IKC)                                                                :: range(abs(stop - start) + 1_IKC)
    end function
#endif

#if IK1_ENABLED
    PURE module function getRangeUnit_D1_IK1(start, stop) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeUnit_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)                                        :: start, stop
        integer(IKC)                                                                :: range(abs(stop - start) + 1_IKC)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getRangeUnit_D1_RK5(start, stop) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeUnit_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , value                                             :: start, stop
        real(RKC)               , allocatable                                       :: range(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getRangeUnit_D1_RK4(start, stop) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeUnit_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , value                                             :: start, stop
        real(RKC)               , allocatable                                       :: range(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getRangeUnit_D1_RK3(start, stop) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeUnit_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , value                                             :: start, stop
        real(RKC)               , allocatable                                       :: range(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getRangeUnit_D1_RK2(start, stop) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeUnit_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , value                                             :: start, stop
        real(RKC)               , allocatable                                       :: range(:)
    end function
#endif

#if RK1_ENABLED
    PURE module function getRangeUnit_D1_RK1(start, stop) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeUnit_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , value                                             :: start, stop
        real(RKC)               , allocatable                                       :: range(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRangeStep_D0_SK5(start, stop, step) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeStep_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        integer(IK)             , intent(in)                                        :: step
        character(1,SKC)        , intent(in)                                        :: start, stop
        character(max(0, 1 + floor(real(ichar(stop) - ichar(start)) / step)),SKC)   :: range
    end function
#endif

#if SK4_ENABLED
    PURE module function getRangeStep_D0_SK4(start, stop, step) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeStep_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        integer(IK)             , intent(in)                                        :: step
        character(1,SKC)        , intent(in)                                        :: start, stop
        character(max(0, 1 + floor(real(ichar(stop) - ichar(start)) / step)),SKC)   :: range
    end function
#endif

#if SK3_ENABLED
    PURE module function getRangeStep_D0_SK3(start, stop, step) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeStep_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        integer(IK)             , intent(in)                                        :: step
        character(1,SKC)        , intent(in)                                        :: start, stop
        character(max(0, 1 + floor(real(ichar(stop) - ichar(start)) / step)),SKC)   :: range
    end function
#endif

#if SK2_ENABLED
    PURE module function getRangeStep_D0_SK2(start, stop, step) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeStep_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        integer(IK)             , intent(in)                                        :: step
        character(1,SKC)        , intent(in)                                        :: start, stop
        character(max(0, 1 + floor(real(ichar(stop) - ichar(start)) / step)),SKC)   :: range
    end function
#endif

#if SK1_ENABLED
    PURE module function getRangeStep_D0_SK1(start, stop, step) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeStep_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        integer(IK)             , intent(in)                                        :: step
        character(1,SKC)        , intent(in)                                        :: start, stop
        character(max(0, 1 + floor(real(ichar(stop) - ichar(start)) / step)),SKC)   :: range
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getRangeStep_D1_IK5(start, stop, step) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeStep_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)                                        :: start, stop, step
        integer(IKC)                                                                :: range(max(0_IKC, 1_IKC + floor(real(stop - start) / real(step), kind = IKC)))
        ! (merge(1_IKC + (stop - start) / step, 0_IKC, (stop - start >= 0_IKC .and. step > 0_IKC) .or. (stop - start <= 0_IKC .and. step < 0_IKC)))
    end function
#endif

#if IK4_ENABLED
    PURE module function getRangeStep_D1_IK4(start, stop, step) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeStep_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)                                        :: start, stop, step
        integer(IKC)                                                                :: range(max(0_IKC, 1_IKC + floor(real(stop - start) / real(step), kind = IKC)))
    end function
#endif

#if IK3_ENABLED
    PURE module function getRangeStep_D1_IK3(start, stop, step) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeStep_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)                                        :: start, stop, step
        integer(IKC)                                                                :: range(max(0_IKC, 1_IKC + floor(real(stop - start) / real(step), kind = IKC)))
    end function
#endif

#if IK2_ENABLED
    PURE module function getRangeStep_D1_IK2(start, stop, step) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeStep_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)                                        :: start, stop, step
        integer(IKC)                                                                :: range(max(0_IKC, 1_IKC + floor(real(stop - start) / real(step), kind = IKC)))
    end function
#endif

#if IK1_ENABLED
    PURE module function getRangeStep_D1_IK1(start, stop, step) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeStep_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)                                        :: start, stop, step
        integer(IKC)                                                                :: range(max(0_IKC, 1_IKC + floor(real(stop - start) / real(step), kind = IKC)))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getRangeStep_D1_RK5(start, stop, step) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeStep_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)                                        :: start, stop, step
        real(RKC)                                                                   :: range(max(0_IK, 1_IK + floor((stop - start) / step, kind = IK)))
    end function
#endif

#if RK4_ENABLED
    PURE module function getRangeStep_D1_RK4(start, stop, step) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeStep_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)                                        :: start, stop, step
        real(RKC)                                                                   :: range(max(0_IK, 1_IK + floor((stop - start) / step, kind = IK)))
    end function
#endif

#if RK3_ENABLED
    PURE module function getRangeStep_D1_RK3(start, stop, step) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeStep_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)                                        :: start, stop, step
        real(RKC)                                                                   :: range(max(0_IK, 1_IK + floor((stop - start) / step, kind = IK)))
    end function
#endif

#if RK2_ENABLED
    PURE module function getRangeStep_D1_RK2(start, stop, step) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeStep_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)                                        :: start, stop, step
        real(RKC)                                                                   :: range(max(0_IK, 1_IK + floor((stop - start) / step, kind = IK)))
    end function
#endif

#if RK1_ENABLED
    PURE module function getRangeStep_D1_RK1(start, stop, step) result(range)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeStep_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)                                        :: start, stop, step
        real(RKC)                                                                   :: range(max(0_IK, 1_IK + floor((stop - start) / step, kind = IK)))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module function getRangeUnit_D2_IK5(start, stop) result(range)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeUnit_D2_IK5
!#endif
!        use pm_kind, only: IKC => IK5
!        integer(IKC)            , intent(in)                :: start(2), stop(2)
!        integer(IKC)                                        :: range(2, max(0_IK, stop(1) - start(1) + 1_IK), max(0_IK, stop(2) - start(2) + 1_IK))
!    end function
!#endif
!
!#if IK4_ENABLED
!    PURE module function getRangeUnit_D2_IK4(start, stop) result(range)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeUnit_D2_IK4
!#endif
!        use pm_kind, only: IKC => IK4
!        integer(IKC)            , intent(in)                :: start(2), stop(2)
!        integer(IKC)                                        :: range(2, max(0_IK, stop(1) - start(1) + 1_IK), max(0_IK, stop(2) - start(2) + 1_IK))
!    end function
!#endif
!
!#if IK3_ENABLED
!    PURE module function getRangeUnit_D2_IK3(start, stop) result(range)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeUnit_D2_IK3
!#endif
!        use pm_kind, only: IKC => IK3
!        integer(IKC)            , intent(in)                :: start(2), stop(2)
!        integer(IKC)                                        :: range(2, max(0_IK, stop(1) - start(1) + 1_IK), max(0_IK, stop(2) - start(2) + 1_IK))
!    end function
!#endif
!
!#if IK2_ENABLED
!    PURE module function getRangeUnit_D2_IK2(start, stop) result(range)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeUnit_D2_IK2
!#endif
!        use pm_kind, only: IKC => IK2
!        integer(IKC)            , intent(in)                :: start(2), stop(2)
!        integer(IKC)                                        :: range(2, max(0_IK, stop(1) - start(1) + 1_IK), max(0_IK, stop(2) - start(2) + 1_IK))
!    end function
!#endif
!
!#if IK1_ENABLED
!    PURE module function getRangeUnit_D2_IK1(start, stop) result(range)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeUnit_D2_IK1
!#endif
!        use pm_kind, only: IKC => IK1
!        integer(IKC)            , intent(in)                :: start(2), stop(2)
!        integer(IKC)                                        :: range(2, max(0_IK, stop(1) - start(1) + 1_IK), max(0_IK, stop(2) - start(2) + 1_IK))
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module function getRangeStep_D2_IK5(start, stop, step) result(range)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeStep_D2_IK5
!#endif
!        use pm_kind, only: IKC => IK5
!        integer(IKC)            , intent(in)                :: start(2), stop(2), step(2)
!        integer(IKC)                                        :: range(2, max(0_IKC, 1_IKC + floor(real(stop(1) - start(1)) / real(step(1)), kind = IKC)), max(0_IKC, 1_IKC + floor(real(stop(2) - start(2)) / real(step(2)), kind = IKC)))
!        ! (merge(1_IKC + (stop(1) - start(1)) / step, 0_IKC, (stop(1) - start(1) >= 0_IKC .and. step > 0_IKC) .or. (stop(1) - start(1) <= 0_IKC .and. step < 0_IKC)))
!    end function
!#endif
!
!#if IK4_ENABLED
!    PURE module function getRangeStep_D2_IK4(start, stop, step) result(range)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeStep_D2_IK4
!#endif
!        use pm_kind, only: IKC => IK4
!        integer(IKC)            , intent(in)                :: start(2), stop(2), step(2)
!        integer(IKC)                                        :: range(2, max(0_IKC, 1_IKC + floor(real(stop(1) - start(1)) / real(step(1)), kind = IKC)), max(0_IKC, 1_IKC + floor(real(stop(2) - start(2)) / real(step(2)), kind = IKC)))
!    end function
!#endif
!
!#if IK3_ENABLED
!    PURE module function getRangeStep_D2_IK3(start, stop, step) result(range)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeStep_D2_IK3
!#endif
!        use pm_kind, only: IKC => IK3
!        integer(IKC)            , intent(in)                :: start(2), stop(2), step(2)
!        integer(IKC)                                        :: range(2, max(0_IKC, 1_IKC + floor(real(stop(1) - start(1)) / real(step(1)), kind = IKC)), max(0_IKC, 1_IKC + floor(real(stop(2) - start(2)) / real(step(2)), kind = IKC)))
!    end function
!#endif
!
!#if IK2_ENABLED
!    PURE module function getRangeStep_D2_IK2(start, stop, step) result(range)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeStep_D2_IK2
!#endif
!        use pm_kind, only: IKC => IK2
!        integer(IKC)            , intent(in)                :: start(2), stop(2), step(2)
!        integer(IKC)                                        :: range(2, max(0_IKC, 1_IKC + floor(real(stop(1) - start(1)) / real(step(1)), kind = IKC)), max(0_IKC, 1_IKC + floor(real(stop(2) - start(2)) / real(step(2)), kind = IKC)))
!    end function
!#endif
!
!#if IK1_ENABLED
!    PURE module function getRangeStep_D2_IK1(start, stop, step) result(range)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getRangeStep_D2_IK1
!#endif
!        use pm_kind, only: IKC => IK1
!        integer(IKC)            , intent(in)                :: start(2), stop(2), step(2)
!        integer(IKC)                                        :: range(2, max(0_IKC, 1_IKC + floor(real(stop(1) - start(1)) / real(step(1)), kind = IKC)), max(0_IKC, 1_IKC + floor(real(stop(2) - start(2)) / real(step(2)), kind = IKC)))
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return evenly spaced `integer` or `character` sequence starting at `start` with jump sizes of `step`.<br>
    !>
    !>  \details
    !>  The value of the last element of the output array is determined by the
    !>  size of the array: `range(size(range)) = start + step * (size(range) - 1)`.<br>
    !>  If `range` is of type `character`, then replace `size` with `len` in the above formula.<br>
    !>
    !>  \param[out] range   :   The output sequence that can be,
    !>                          <ol>
    !>                              <li>    a scalar `character` of kind \SKALL of length type parameter `1`,<br>
    !>                              <li>    a `contiguous` array of rank `1` of type `integer` of kind \IKALL,<br>
    !>                          </ol>
    !>                          containing the evenly-spaced sequence within the interval starting with `start` and with jumps of size `step`.<br>
    !>  \param[in]  stop    :   The input scalar of the same type and kind as `start` representing the end point of the output `range`.<br>
    !>                          If `stop` is of type `character`, then it must be of length type parameter `1`.<br>
    !>  \param[in]  step    :   The **non-zero** input scalar of type `integer` representing the size of the
    !>                          interval between adjacent values in the output `range`.<br>
    !>                          <ol>
    !>                              <li>    If `range` is of type `character`, then `step` must have the default `integer` kind \IK.<br>
    !>                              <li>    If `range` is of type `integer`, then `step` must have the same type kind parameter as `range`.<br>
    !>                          </ol>
    !>                          (**optional**, default = `1`)
    !>
    !>  \interface{setRange}
    !>  \code{.F90}
    !>
    !>      use pm_arrayRange, only: setRange
    !>
    !>      ! Integer sequence.
    !>
    !>      call setRange(range(:), start)
    !>      call setRange(range(:), start, step)
    !>
    !>      ! Character sequence.
    !>
    !>      call setRange(range, start)
    !>      call setRange(range, start, step)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `step /= 0` must hold for the corresponding input arguments.<br>
    !>  The condition `len(start) == 1` must hold for the corresponding input `character` arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getRange](@ref pm_arrayRange::getRange)<br>
    !>  [setRange](@ref pm_arrayRange::setRange)<br>
    !>  [getLinSpace](@ref pm_arraySpace::getLinSpace)<br>
    !>  [setLinSpace](@ref pm_arraySpace::setLinSpace)<br>
    !>  [getLogSpace](@ref pm_arraySpace::getLogSpace)<br>
    !>  [setLogSpace](@ref pm_arraySpace::setLogSpace)<br>
    !>
    !>  \example{setRange}
    !>  \include{lineno} example/pm_arrayRange/setRange/main.F90
    !>  \compilef{setRange}
    !>  \output{setRange}
    !>  \include{lineno} example/pm_arrayRange/setRange/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayRange](@ref test_pm_arrayRange)<br>
    !>
    !>  \final{setRange}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 08:49 PM, August 10, 2021, Dallas, TX
    interface setRange

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRangeUnit_D0_SK5(range, start)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeUnit_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(out)                   :: range
        character(1,SKC)        , intent(in)                    :: start
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRangeUnit_D0_SK4(range, start)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeUnit_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(out)                   :: range
        character(1,SKC)        , intent(in)                    :: start
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRangeUnit_D0_SK3(range, start)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeUnit_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(out)                   :: range
        character(1,SKC)        , intent(in)                    :: start
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRangeUnit_D0_SK2(range, start)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeUnit_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(out)                   :: range
        character(1,SKC)        , intent(in)                    :: start
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRangeUnit_D0_SK1(range, start)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeUnit_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(out)                   :: range
        character(1,SKC)        , intent(in)                    :: start
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRangeUnit_D1_IK5(range, start)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeUnit_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(out)   , contiguous    :: range(:)
        integer(IKC)            , intent(in)                    :: start
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRangeUnit_D1_IK4(range, start)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeUnit_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(out)   , contiguous    :: range(:)
        integer(IKC)            , intent(in)                    :: start
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRangeUnit_D1_IK3(range, start)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeUnit_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(out)   , contiguous    :: range(:)
        integer(IKC)            , intent(in)                    :: start
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRangeUnit_D1_IK2(range, start)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeUnit_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(out)   , contiguous    :: range(:)
        integer(IKC)            , intent(in)                    :: start
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRangeUnit_D1_IK1(range, start)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeUnit_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(out)   , contiguous    :: range(:)
        integer(IKC)            , intent(in)                    :: start
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRangeUnit_D1_RK5(range, start)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeUnit_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(out)   , contiguous    :: range(:)
        real(RKC)               , intent(in)                    :: start
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRangeUnit_D1_RK4(range, start)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeUnit_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(out)   , contiguous    :: range(:)
        real(RKC)               , intent(in)                    :: start
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRangeUnit_D1_RK3(range, start)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeUnit_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(out)   , contiguous    :: range(:)
        real(RKC)               , intent(in)                    :: start
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRangeUnit_D1_RK2(range, start)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeUnit_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(out)   , contiguous    :: range(:)
        real(RKC)               , intent(in)                    :: start
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRangeUnit_D1_RK1(range, start)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeUnit_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(out)   , contiguous    :: range(:)
        real(RKC)               , intent(in)                    :: start
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRangeStep_D0_SK5(range, start, step)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeStep_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(out)                   :: range
        character(1,SKC)        , intent(in)                    :: start
        integer(IK)             , intent(in)                    :: step
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRangeStep_D0_SK4(range, start, step)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeStep_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(out)                   :: range
        character(1,SKC)        , intent(in)                    :: start
        integer(IK)             , intent(in)                    :: step
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRangeStep_D0_SK3(range, start, step)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeStep_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(out)                   :: range
        character(1,SKC)        , intent(in)                    :: start
        integer(IK)             , intent(in)                    :: step
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRangeStep_D0_SK2(range, start, step)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeStep_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(out)                   :: range
        character(1,SKC)        , intent(in)                    :: start
        integer(IK)             , intent(in)                    :: step
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRangeStep_D0_SK1(range, start, step)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeStep_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(out)                   :: range
        character(1,SKC)        , intent(in)                    :: start
        integer(IK)             , intent(in)                    :: step
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRangeStep_D1_IK5(range, start, step)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeStep_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(out)   , contiguous    :: range(:)
        integer(IKC)            , intent(in)                    :: start, step
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRangeStep_D1_IK4(range, start, step)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeStep_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(out)   , contiguous    :: range(:)
        integer(IKC)            , intent(in)                    :: start, step
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRangeStep_D1_IK3(range, start, step)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeStep_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(out)   , contiguous    :: range(:)
        integer(IKC)            , intent(in)                    :: start, step
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRangeStep_D1_IK2(range, start, step)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeStep_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(out)   , contiguous    :: range(:)
        integer(IKC)            , intent(in)                    :: start, step
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRangeStep_D1_IK1(range, start, step)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeStep_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(out)   , contiguous    :: range(:)
        integer(IKC)            , intent(in)                    :: start, step
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRangeStep_D1_RK5(range, start, step)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeStep_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(out)   , contiguous    :: range(:)
        real(RKC)               , intent(in)                    :: start, step
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRangeStep_D1_RK4(range, start, step)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeStep_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(out)   , contiguous    :: range(:)
        real(RKC)               , intent(in)                    :: start, step
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRangeStep_D1_RK3(range, start, step)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeStep_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(out)   , contiguous    :: range(:)
        real(RKC)               , intent(in)                    :: start, step
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRangeStep_D1_RK2(range, start, step)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeStep_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(out)   , contiguous    :: range(:)
        real(RKC)               , intent(in)                    :: start, step
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRangeStep_D1_RK1(range, start, step)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRangeStep_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(out)   , contiguous    :: range(:)
        real(RKC)               , intent(in)                    :: start, step
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayRange ! LCOV_EXCL_LINE