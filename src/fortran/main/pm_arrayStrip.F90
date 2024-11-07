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
!>  This module contains procedures and generic interfaces for stripping a given pattern
!>  from the left and right ends of an array of arbitrary intrinsic type and kind.<br>
!>
!>  \details
!>  The functionalities of the [getStripped](@ref pm_arrayStrip::getStripped)
!>  generic interface of this module are similar to the `strip()` method of Python strings.<br>
!>  However, keep in mind that returning an `allocatable` stripped copy of the input can be an expensive operation.<br>
!>  If performance matters, use instead the generic interfaces [getSIL](@ref pm_arrayStrip::getSIL) for left-stripping
!>  [getSIR](@ref pm_arrayStrip::getSIR) for right-stripping.<br>
!>  These **indexing** avoid the unnecessary and costly `allocation`s,
!>  potentially leading to much faster runtime performance.<br>
!>
!>  \devnote
!>  The `subroutine` equivalents of the `function` interface
!>  [getStripped](@ref pm_arrayStrip::getStripped) of this module were deemed unnecessary.<br>
!>  The functional interfaces exist solely for improved flexibility.<br>
!>  If performance matters, [getSIL](@ref pm_arrayStrip::getSIL) and [getSIR](@ref pm_arrayStrip::getSIR)
!>  can be readily used to eliminate the avoidable allocations within the procedures.<br>
!>
!>  \see
!>  [pm_arraySplit](@ref pm_arraySplit)<br>
!>  [pm_arrayRemove](@ref pm_arrayRemove)<br>
!>
!>  \test
!>  [test_pm_arrayStrip](@ref test_pm_arrayStrip)
!>
!>  \todo
!>  \plow
!>  A benchmark comparing the performance of indexing vs. allocating a modified array should be provided here.
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, 3:48 AM, October 13, 2017, Institute for Computational Engineering and Sciences, Austin, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayStrip

    use pm_kind, only: SK, IK, LK
    use pm_array, only: leftRight_type, leftRight
    use pm_array, only: right_type, right
    use pm_array, only: left_type, left

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_arrayStrip"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the input `array` modified such that all instances of
    !>  the input `pattern` are stripped from its both sides or only the side specified.
    !>
    !>  \param[in]  array       :   The input `contiguous` vector of either <br>
    !>                              <ul>
    !>                                  <li>    type `character` of kind \SKALL, or <br>
    !>                                  <li>    type `logical` of kind \LKALL, or <br>
    !>                                  <li>    type `integer` of kind \IKALL, or <br>
    !>                                  <li>    type `complex` of kind \CKALL, or <br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ul>
    !>                              or,
    !>                              <ul>
    !>                                  <li>    a scalar assumed-length `character` of kind \SKALL, <br>
    !>                              </ul>
    !>  \param[in]  pattern     :   The input object of the same type and kind as the input array, of rank similar
    !>                              to or lower than that of `array`, whose value is to be stripped from the specified `side` of the `array`.<br>
    !>  \param      iseq        :   The `external` user-specified function that takes either two input assumed-length `character` arguments
    !>                              (if the input `array` is also an assumed-length `character`) or two array-valued **explicit-shape**
    !>                              arguments of the same type and kind as the input `array`.<br>
    !>                              It must return a scalar of type `logical` of default kind \LK that is `.true.` if all elements of
    !>                              the two input arguments are equivalent (e.g., equal) according to the user-defined criterion, otherwise, it is `.false.`.<br>
    !>                              The the input `array` is array-valued, then the last argument to `iseq` is the length of the input `pattern`.<br>
    !>                              The following illustrates the generic interface of `iseq` where `pattern` is array-valued,
    !>                              \code{.F90}
    !>                                  function iseq(Segment, pattern, lenPattern) result(equivalent)
    !>                                      use pm_kind, only: IK, LK
    !>                                      integer(IK) , intent(in)    :: lenPattern
    !>                                      TYPE(KIND)  , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                      logical(LK)                 :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
    !>                              \code{.F90}
    !>                                  character(*, SK), intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  integer(IK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  logical(LK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  complex(CK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  real(RK)        , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                              \endcode
    !>                              where the kinds `SKG`, `IKG`, `LKG`, `CKG`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              The following illustrates the generic interface of `iseq` where `pattern` is scalar-valued (**including Fortran scalar strings**),
    !>                              \code{.F90}
    !>                                  function iseq(segment, pattern, side) result(equivalent)
    !>                                      use pm_kind, only: LK
    !>                                      TYPE(KIND)  , intent(in)    :: segment, pattern
    !>                                      logical(LK)                 :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
    !>                              \code{.F90}
    !>                                  use pm_kind, only: SK, IK, LK, CK, RK
    !>                                  character(*, SK), intent(in)    :: segment, pattern
    !>                                  integer(IK)     , intent(in)    :: segment, pattern
    !>                                  logical(LK)     , intent(in)    :: segment, pattern
    !>                                  complex(CK)     , intent(in)    :: segment, pattern
    !>                                  real(RK)        , intent(in)    :: segment, pattern
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              This user-defined equivalence check is extremely useful where an equivalence test other than exact identity is needed,
    !>                              for example, when the array segments should match the input `pattern` only within a given threshold or,
    !>                              when the case-sensitivity in character comparisons do not matter.<br>
    !>                              In such cases, user can define a custom equivalence criterion within the user-defined external function `iseq` to achieve the goal.<br>
    !>                              (**optional**, the default equivalence operator is `.eqv.` if the input `array` is `logical`, otherwise `==`.)
    !>  \param      side        :   The input scalar constant that be any of the following:<br>
    !>                              <ol>
    !>                                  <li>    The scalar constant [left](@ref pm_array::left) or an object of type [left_type](@ref pm_array::left_type),
    !>                                          implying that only the left side of the input array should be stripped.<br>
    !>                                  <li>    The scalar constant [right](@ref pm_array::right) or an object of type [right_type](@ref pm_array::right_type),
    !>                                          implying that only the right side of the input array should be stripped.<br>
    !>                                  <li>    The scalar constant [leftRight](@ref pm_array::leftRight) or an object of type [leftRight_type](@ref pm_array::leftRight_type),
    !>                                          implying that both left and right sides of the input array should be stripped.<br>
    !>                              </ol>
    !>                              (**optional**, default = [leftRight](@ref pm_array::leftRight))
    !>
    !>  \return
    !>  `arrayStripped`         :   The output `allocatable` array of the same type, kind, rank, shape, and size as the input `array`,
    !>                              containing the input `array` where all instances of `pattern` are stripped from the specified `side`.<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_arrayStrip, only: getStripped, left, right, leftRight
    !>
    !>      arrayStripped = getStripped(array, pattern, side = side) ! scalar strings
    !>      arrayStripped = getStripped(array(:), pattern, side = side) ! all intrinsic types
    !>      arrayStripped = getStripped(array(:), pattern(:), side = side) ! all intrinsic types
    !>      arrayStripped = getStripped(array, pattern, iseq, side = side) ! scalar strings
    !>      arrayStripped = getStripped(array(:), pattern, iseq, side = side) ! all intrinsic types
    !>      arrayStripped = getStripped(array(:), pattern(:), iseq, side = side) ! all intrinsic types
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The procedures under this generic interface are `impure` when the user-specified `external` procedure `iseq` is specified as input argument.
    !>
    !>  \warning
    !>  Note that in Fortran, trailing blanks are ignored in character equivalence checks, that is, `"Fortran" == "Fortran "` yields `.true.`.
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getSIL](@ref pm_arrayStrip::getSIL)<br>
    !>  [getSIR](@ref pm_arrayStrip::getSIR)<br>
    !>  [getStripped](@ref pm_arrayStrip::getStripped)<br>
    !>  [setReplaced](@ref pm_arrayReplace::setReplaced)<br>
    !>  [getReplaced](@ref pm_arrayReplace::getReplaced)<br>
    !>  [setInserted](@ref pm_arrayInsert::setInserted)<br>
    !>  [setSplit](@ref pm_arraySplit::setSplit)<br>
    !>  [setLoc](@ref pm_arrayFind::setLoc)<br>
    !>  [getLoc](@ref pm_arrayFind::getLoc)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_arrayStrip/getStripped/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_arrayStrip/getStripped/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayStrip](@ref test_pm_arrayStrip)
    !>
    !>  \final{getStripped}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! side: both

    interface getStripped

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getStrippedDefComB_D0_D0_SK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK4_ENABLED
    PURE module function getStrippedDefComB_D0_D0_SK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK3_ENABLED
    PURE module function getStrippedDefComB_D0_D0_SK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK2_ENABLED
    PURE module function getStrippedDefComB_D0_D0_SK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK1_ENABLED
    PURE module function getStrippedDefComB_D0_D0_SK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getStrippedDefComB_D1_D0_SK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK4_ENABLED
    PURE module function getStrippedDefComB_D1_D0_SK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK3_ENABLED
    PURE module function getStrippedDefComB_D1_D0_SK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK2_ENABLED
    PURE module function getStrippedDefComB_D1_D0_SK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK1_ENABLED
    PURE module function getStrippedDefComB_D1_D0_SK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getStrippedDefComB_D1_D0_IK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if IK4_ENABLED
    PURE module function getStrippedDefComB_D1_D0_IK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if IK3_ENABLED
    PURE module function getStrippedDefComB_D1_D0_IK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if IK2_ENABLED
    PURE module function getStrippedDefComB_D1_D0_IK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if IK1_ENABLED
    PURE module function getStrippedDefComB_D1_D0_IK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getStrippedDefComB_D1_D0_LK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if LK4_ENABLED
    PURE module function getStrippedDefComB_D1_D0_LK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if LK3_ENABLED
    PURE module function getStrippedDefComB_D1_D0_LK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if LK2_ENABLED
    PURE module function getStrippedDefComB_D1_D0_LK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if LK1_ENABLED
    PURE module function getStrippedDefComB_D1_D0_LK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getStrippedDefComB_D1_D0_CK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if CK4_ENABLED
    PURE module function getStrippedDefComB_D1_D0_CK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if CK3_ENABLED
    PURE module function getStrippedDefComB_D1_D0_CK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if CK2_ENABLED
    PURE module function getStrippedDefComB_D1_D0_CK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if CK1_ENABLED
    PURE module function getStrippedDefComB_D1_D0_CK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getStrippedDefComB_D1_D0_RK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if RK4_ENABLED
    PURE module function getStrippedDefComB_D1_D0_RK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if RK3_ENABLED
    PURE module function getStrippedDefComB_D1_D0_RK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if RK2_ENABLED
    PURE module function getStrippedDefComB_D1_D0_RK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if RK1_ENABLED
    PURE module function getStrippedDefComB_D1_D0_RK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getStrippedCusComLR_D0_D0_SK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK4_ENABLED
    module function getStrippedCusComLR_D0_D0_SK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK3_ENABLED
    module function getStrippedCusComLR_D0_D0_SK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK2_ENABLED
    module function getStrippedCusComLR_D0_D0_SK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK1_ENABLED
    module function getStrippedCusComLR_D0_D0_SK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getStrippedCusComLR_D1_D0_SK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK4_ENABLED
    module function getStrippedCusComLR_D1_D0_SK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK3_ENABLED
    module function getStrippedCusComLR_D1_D0_SK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK2_ENABLED
    module function getStrippedCusComLR_D1_D0_SK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK1_ENABLED
    module function getStrippedCusComLR_D1_D0_SK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getStrippedCusComLR_D1_D0_IK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if IK4_ENABLED
    module function getStrippedCusComLR_D1_D0_IK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if IK3_ENABLED
    module function getStrippedCusComLR_D1_D0_IK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if IK2_ENABLED
    module function getStrippedCusComLR_D1_D0_IK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if IK1_ENABLED
    module function getStrippedCusComLR_D1_D0_IK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getStrippedCusComLR_D1_D0_LK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if LK4_ENABLED
    module function getStrippedCusComLR_D1_D0_LK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if LK3_ENABLED
    module function getStrippedCusComLR_D1_D0_LK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if LK2_ENABLED
    module function getStrippedCusComLR_D1_D0_LK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if LK1_ENABLED
    module function getStrippedCusComLR_D1_D0_LK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getStrippedCusComLR_D1_D0_CK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if CK4_ENABLED
    module function getStrippedCusComLR_D1_D0_CK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if CK3_ENABLED
    module function getStrippedCusComLR_D1_D0_CK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if CK2_ENABLED
    module function getStrippedCusComLR_D1_D0_CK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if CK1_ENABLED
    module function getStrippedCusComLR_D1_D0_CK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getStrippedCusComLR_D1_D0_RK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if RK4_ENABLED
    module function getStrippedCusComLR_D1_D0_RK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if RK3_ENABLED
    module function getStrippedCusComLR_D1_D0_RK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if RK2_ENABLED
    module function getStrippedCusComLR_D1_D0_RK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if RK1_ENABLED
    module function getStrippedCusComLR_D1_D0_RK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getStrippedDefComB_D1_D1_SK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK4_ENABLED
    PURE module function getStrippedDefComB_D1_D1_SK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK3_ENABLED
    PURE module function getStrippedDefComB_D1_D1_SK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK2_ENABLED
    PURE module function getStrippedDefComB_D1_D1_SK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK1_ENABLED
    PURE module function getStrippedDefComB_D1_D1_SK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getStrippedDefComB_D1_D1_IK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if IK4_ENABLED
    PURE module function getStrippedDefComB_D1_D1_IK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if IK3_ENABLED
    PURE module function getStrippedDefComB_D1_D1_IK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if IK2_ENABLED
    PURE module function getStrippedDefComB_D1_D1_IK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if IK1_ENABLED
    PURE module function getStrippedDefComB_D1_D1_IK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getStrippedDefComB_D1_D1_LK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if LK4_ENABLED
    PURE module function getStrippedDefComB_D1_D1_LK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if LK3_ENABLED
    PURE module function getStrippedDefComB_D1_D1_LK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if LK2_ENABLED
    PURE module function getStrippedDefComB_D1_D1_LK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if LK1_ENABLED
    PURE module function getStrippedDefComB_D1_D1_LK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getStrippedDefComB_D1_D1_CK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if CK4_ENABLED
    PURE module function getStrippedDefComB_D1_D1_CK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if CK3_ENABLED
    PURE module function getStrippedDefComB_D1_D1_CK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if CK2_ENABLED
    PURE module function getStrippedDefComB_D1_D1_CK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if CK1_ENABLED
    PURE module function getStrippedDefComB_D1_D1_CK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getStrippedDefComB_D1_D1_RK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if RK4_ENABLED
    PURE module function getStrippedDefComB_D1_D1_RK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if RK3_ENABLED
    PURE module function getStrippedDefComB_D1_D1_RK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if RK2_ENABLED
    PURE module function getStrippedDefComB_D1_D1_RK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if RK1_ENABLED
    PURE module function getStrippedDefComB_D1_D1_RK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComB_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getStrippedCusComLR_D1_D1_SK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK4_ENABLED
    module function getStrippedCusComLR_D1_D1_SK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK3_ENABLED
    module function getStrippedCusComLR_D1_D1_SK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK2_ENABLED
    module function getStrippedCusComLR_D1_D1_SK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if SK1_ENABLED
    module function getStrippedCusComLR_D1_D1_SK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getStrippedCusComLR_D1_D1_IK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if IK4_ENABLED
    module function getStrippedCusComLR_D1_D1_IK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if IK3_ENABLED
    module function getStrippedCusComLR_D1_D1_IK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if IK2_ENABLED
    module function getStrippedCusComLR_D1_D1_IK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if IK1_ENABLED
    module function getStrippedCusComLR_D1_D1_IK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getStrippedCusComLR_D1_D1_LK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if LK4_ENABLED
    module function getStrippedCusComLR_D1_D1_LK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if LK3_ENABLED
    module function getStrippedCusComLR_D1_D1_LK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if LK2_ENABLED
    module function getStrippedCusComLR_D1_D1_LK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if LK1_ENABLED
    module function getStrippedCusComLR_D1_D1_LK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getStrippedCusComLR_D1_D1_CK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if CK4_ENABLED
    module function getStrippedCusComLR_D1_D1_CK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if CK3_ENABLED
    module function getStrippedCusComLR_D1_D1_CK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if CK2_ENABLED
    module function getStrippedCusComLR_D1_D1_CK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if CK1_ENABLED
    module function getStrippedCusComLR_D1_D1_CK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getStrippedCusComLR_D1_D1_RK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(logical(LK))                              :: iseq
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if RK4_ENABLED
    module function getStrippedCusComLR_D1_D1_RK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(logical(LK))                              :: iseq
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if RK3_ENABLED
    module function getStrippedCusComLR_D1_D1_RK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(logical(LK))                              :: iseq
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if RK2_ENABLED
    module function getStrippedCusComLR_D1_D1_RK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(logical(LK))                              :: iseq
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

#if RK1_ENABLED
    module function getStrippedCusComLR_D1_D1_RK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComLR_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(logical(LK))                              :: iseq
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(leftRight_type)    , intent(in), optional      :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! side: left

    interface getStripped

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getStrippedDefComSL_D0_D0_SK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK4_ENABLED
    PURE module function getStrippedDefComSL_D0_D0_SK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK3_ENABLED
    PURE module function getStrippedDefComSL_D0_D0_SK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK2_ENABLED
    PURE module function getStrippedDefComSL_D0_D0_SK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK1_ENABLED
    PURE module function getStrippedDefComSL_D0_D0_SK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(left_type)         , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_SK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK4_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_SK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK3_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_SK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK2_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_SK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK1_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_SK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_IK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if IK4_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_IK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if IK3_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_IK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if IK2_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_IK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if IK1_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_IK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_LK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if LK4_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_LK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if LK3_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_LK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if LK2_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_LK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if LK1_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_LK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_CK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if CK4_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_CK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if CK3_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_CK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if CK2_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_CK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if CK1_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_CK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_RK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if RK4_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_RK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if RK3_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_RK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if RK2_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_RK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if RK1_ENABLED
    PURE module function getStrippedDefComSL_D1_D0_RK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getStrippedCusComSL_D0_D0_SK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK4_ENABLED
    module function getStrippedCusComSL_D0_D0_SK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK3_ENABLED
    module function getStrippedCusComSL_D0_D0_SK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK2_ENABLED
    module function getStrippedCusComSL_D0_D0_SK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK1_ENABLED
    module function getStrippedCusComSL_D0_D0_SK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(left_type)         , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getStrippedCusComSL_D1_D0_SK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK4_ENABLED
    module function getStrippedCusComSL_D1_D0_SK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK3_ENABLED
    module function getStrippedCusComSL_D1_D0_SK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK2_ENABLED
    module function getStrippedCusComSL_D1_D0_SK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK1_ENABLED
    module function getStrippedCusComSL_D1_D0_SK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getStrippedCusComSL_D1_D0_IK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if IK4_ENABLED
    module function getStrippedCusComSL_D1_D0_IK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if IK3_ENABLED
    module function getStrippedCusComSL_D1_D0_IK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if IK2_ENABLED
    module function getStrippedCusComSL_D1_D0_IK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if IK1_ENABLED
    module function getStrippedCusComSL_D1_D0_IK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getStrippedCusComSL_D1_D0_LK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if LK4_ENABLED
    module function getStrippedCusComSL_D1_D0_LK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if LK3_ENABLED
    module function getStrippedCusComSL_D1_D0_LK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if LK2_ENABLED
    module function getStrippedCusComSL_D1_D0_LK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if LK1_ENABLED
    module function getStrippedCusComSL_D1_D0_LK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getStrippedCusComSL_D1_D0_CK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if CK4_ENABLED
    module function getStrippedCusComSL_D1_D0_CK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if CK3_ENABLED
    module function getStrippedCusComSL_D1_D0_CK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if CK2_ENABLED
    module function getStrippedCusComSL_D1_D0_CK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if CK1_ENABLED
    module function getStrippedCusComSL_D1_D0_CK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getStrippedCusComSL_D1_D0_RK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if RK4_ENABLED
    module function getStrippedCusComSL_D1_D0_RK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if RK3_ENABLED
    module function getStrippedCusComSL_D1_D0_RK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if RK2_ENABLED
    module function getStrippedCusComSL_D1_D0_RK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if RK1_ENABLED
    module function getStrippedCusComSL_D1_D0_RK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_SK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK4_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_SK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK3_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_SK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK2_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_SK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK1_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_SK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_IK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if IK4_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_IK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if IK3_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_IK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if IK2_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_IK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if IK1_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_IK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_LK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if LK4_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_LK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if LK3_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_LK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if LK2_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_LK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if LK1_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_LK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_CK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if CK4_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_CK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if CK3_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_CK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if CK2_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_CK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if CK1_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_CK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_RK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if RK4_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_RK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if RK3_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_RK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if RK2_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_RK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if RK1_ENABLED
    PURE module function getStrippedDefComSL_D1_D1_RK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSL_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getStrippedCusComSL_D1_D1_SK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK4_ENABLED
    module function getStrippedCusComSL_D1_D1_SK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK3_ENABLED
    module function getStrippedCusComSL_D1_D1_SK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK2_ENABLED
    module function getStrippedCusComSL_D1_D1_SK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if SK1_ENABLED
    module function getStrippedCusComSL_D1_D1_SK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getStrippedCusComSL_D1_D1_IK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if IK4_ENABLED
    module function getStrippedCusComSL_D1_D1_IK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if IK3_ENABLED
    module function getStrippedCusComSL_D1_D1_IK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if IK2_ENABLED
    module function getStrippedCusComSL_D1_D1_IK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if IK1_ENABLED
    module function getStrippedCusComSL_D1_D1_IK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getStrippedCusComSL_D1_D1_LK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if LK4_ENABLED
    module function getStrippedCusComSL_D1_D1_LK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if LK3_ENABLED
    module function getStrippedCusComSL_D1_D1_LK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if LK2_ENABLED
    module function getStrippedCusComSL_D1_D1_LK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if LK1_ENABLED
    module function getStrippedCusComSL_D1_D1_LK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getStrippedCusComSL_D1_D1_CK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if CK4_ENABLED
    module function getStrippedCusComSL_D1_D1_CK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if CK3_ENABLED
    module function getStrippedCusComSL_D1_D1_CK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if CK2_ENABLED
    module function getStrippedCusComSL_D1_D1_CK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if CK1_ENABLED
    module function getStrippedCusComSL_D1_D1_CK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getStrippedCusComSL_D1_D1_RK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(logical(LK))                              :: iseq
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if RK4_ENABLED
    module function getStrippedCusComSL_D1_D1_RK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(logical(LK))                              :: iseq
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if RK3_ENABLED
    module function getStrippedCusComSL_D1_D1_RK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(logical(LK))                              :: iseq
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if RK2_ENABLED
    module function getStrippedCusComSL_D1_D1_RK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(logical(LK))                              :: iseq
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

#if RK1_ENABLED
    module function getStrippedCusComSL_D1_D1_RK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSL_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(logical(LK))                              :: iseq
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(left_type)         , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! side: right

    interface getStripped

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getStrippedDefComSR_D0_D0_SK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK4_ENABLED
    PURE module function getStrippedDefComSR_D0_D0_SK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK3_ENABLED
    PURE module function getStrippedDefComSR_D0_D0_SK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK2_ENABLED
    PURE module function getStrippedDefComSR_D0_D0_SK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK1_ENABLED
    PURE module function getStrippedDefComSR_D0_D0_SK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(right_type)        , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_SK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK4_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_SK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK3_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_SK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK2_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_SK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK1_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_SK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_IK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if IK4_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_IK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if IK3_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_IK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if IK2_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_IK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if IK1_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_IK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_LK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if LK4_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_LK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if LK3_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_LK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if LK2_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_LK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if LK1_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_LK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_CK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if CK4_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_CK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if CK3_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_CK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if CK2_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_CK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if CK1_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_CK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_RK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if RK4_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_RK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if RK3_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_RK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if RK2_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_RK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if RK1_ENABLED
    PURE module function getStrippedDefComSR_D1_D0_RK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getStrippedCusComSR_D0_D0_SK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK4_ENABLED
    module function getStrippedCusComSR_D0_D0_SK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK3_ENABLED
    module function getStrippedCusComSR_D0_D0_SK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK2_ENABLED
    module function getStrippedCusComSR_D0_D0_SK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK1_ENABLED
    module function getStrippedCusComSR_D0_D0_SK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        character(:,SKG)        , allocatable               :: arrayStripped
        type(right_type)        , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getStrippedCusComSR_D1_D0_SK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK4_ENABLED
    module function getStrippedCusComSR_D1_D0_SK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK3_ENABLED
    module function getStrippedCusComSR_D1_D0_SK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK2_ENABLED
    module function getStrippedCusComSR_D1_D0_SK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK1_ENABLED
    module function getStrippedCusComSR_D1_D0_SK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in)                :: pattern
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getStrippedCusComSR_D1_D0_IK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if IK4_ENABLED
    module function getStrippedCusComSR_D1_D0_IK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if IK3_ENABLED
    module function getStrippedCusComSR_D1_D0_IK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if IK2_ENABLED
    module function getStrippedCusComSR_D1_D0_IK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if IK1_ENABLED
    module function getStrippedCusComSR_D1_D0_IK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getStrippedCusComSR_D1_D0_LK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if LK4_ENABLED
    module function getStrippedCusComSR_D1_D0_LK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if LK3_ENABLED
    module function getStrippedCusComSR_D1_D0_LK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if LK2_ENABLED
    module function getStrippedCusComSR_D1_D0_LK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if LK1_ENABLED
    module function getStrippedCusComSR_D1_D0_LK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getStrippedCusComSR_D1_D0_CK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if CK4_ENABLED
    module function getStrippedCusComSR_D1_D0_CK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if CK3_ENABLED
    module function getStrippedCusComSR_D1_D0_CK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if CK2_ENABLED
    module function getStrippedCusComSR_D1_D0_CK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if CK1_ENABLED
    module function getStrippedCusComSR_D1_D0_CK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getStrippedCusComSR_D1_D0_RK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if RK4_ENABLED
    module function getStrippedCusComSR_D1_D0_RK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if RK3_ENABLED
    module function getStrippedCusComSR_D1_D0_RK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if RK2_ENABLED
    module function getStrippedCusComSR_D1_D0_RK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if RK1_ENABLED
    module function getStrippedCusComSR_D1_D0_RK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_SK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK4_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_SK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK3_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_SK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK2_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_SK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK1_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_SK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_IK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if IK4_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_IK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if IK3_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_IK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if IK2_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_IK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if IK1_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_IK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_LK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if LK4_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_LK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if LK3_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_LK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if LK2_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_LK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if LK1_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_LK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_CK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if CK4_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_CK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if CK3_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_CK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if CK2_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_CK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if CK1_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_CK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_RK5(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if RK4_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_RK4(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if RK3_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_RK3(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if RK2_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_RK2(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if RK1_ENABLED
    PURE module function getStrippedDefComSR_D1_D1_RK1(array, pattern, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedDefComSR_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getStrippedCusComSR_D1_D1_SK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK4_ENABLED
    module function getStrippedCusComSR_D1_D1_SK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK3_ENABLED
    module function getStrippedCusComSR_D1_D1_SK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK2_ENABLED
    module function getStrippedCusComSR_D1_D1_SK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if SK1_ENABLED
    module function getStrippedCusComSR_D1_D1_SK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        procedure(logical(LK))                              :: iseq
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        character(len(array,IK),SKG)        , allocatable   :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getStrippedCusComSR_D1_D1_IK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if IK4_ENABLED
    module function getStrippedCusComSR_D1_D1_IK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if IK3_ENABLED
    module function getStrippedCusComSR_D1_D1_IK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if IK2_ENABLED
    module function getStrippedCusComSR_D1_D1_IK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if IK1_ENABLED
    module function getStrippedCusComSR_D1_D1_IK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        procedure(logical(LK))                              :: iseq
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getStrippedCusComSR_D1_D1_LK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if LK4_ENABLED
    module function getStrippedCusComSR_D1_D1_LK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if LK3_ENABLED
    module function getStrippedCusComSR_D1_D1_LK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if LK2_ENABLED
    module function getStrippedCusComSR_D1_D1_LK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if LK1_ENABLED
    module function getStrippedCusComSR_D1_D1_LK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        procedure(logical(LK))                              :: iseq
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        logical(LKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getStrippedCusComSR_D1_D1_CK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if CK4_ENABLED
    module function getStrippedCusComSR_D1_D1_CK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if CK3_ENABLED
    module function getStrippedCusComSR_D1_D1_CK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if CK2_ENABLED
    module function getStrippedCusComSR_D1_D1_CK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if CK1_ENABLED
    module function getStrippedCusComSR_D1_D1_CK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        procedure(logical(LK))                              :: iseq
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        complex(CKG)            , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getStrippedCusComSR_D1_D1_RK5(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(logical(LK))                              :: iseq
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if RK4_ENABLED
    module function getStrippedCusComSR_D1_D1_RK4(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(logical(LK))                              :: iseq
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if RK3_ENABLED
    module function getStrippedCusComSR_D1_D1_RK3(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(logical(LK))                              :: iseq
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if RK2_ENABLED
    module function getStrippedCusComSR_D1_D1_RK2(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(logical(LK))                              :: iseq
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

#if RK1_ENABLED
    module function getStrippedCusComSR_D1_D1_RK1(array, pattern, iseq, side) result(arrayStripped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrippedCusComSR_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(logical(LK))                              :: iseq
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        real(RKG)               , allocatable               :: arrayStripped(:)
        type(right_type)        , intent(in)                :: side
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `index`, the index of the first element of
    !>  the input `array` such that `array(1:index-1)` contains only full repetitions of the user-specified `pattern`.
    !>
    !>  \param[in]  array       :   The input `contiguous` vector of either <br>
    !>                              <ul>
    !>                                  <li>    type `character` of kind \SKALL, or <br>
    !>                                  <li>    type `logical` of kind \LKALL, or <br>
    !>                                  <li>    type `integer` of kind \IKALL, or <br>
    !>                                  <li>    type `complex` of kind \CKALL, or <br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ul>
    !>                              or,
    !>                              <ul>
    !>                                  <li>    a scalar assumed-length `character` of kind \SKALL, <br>
    !>                              </ul>
    !>  \param[in]  pattern     :   The input object of the same type and kind as the input array, of rank similar to or lower than that of `array`, whose value is to be stripped from the beginning of `array`.
    !>  \param      iseq        :   The `external` user-specified function that takes either two input assumed-length `character` arguments
    !>                              (if the input `array` is also an assumed-length `character`) or two array-valued **explicit-shape**
    !>                              arguments of the same type and kind as the input `array`.<br>
    !>                              It must return a scalar of type `logical` of default kind \LK that is `.true.` if all elements of
    !>                              the two input arguments are equivalent (e.g., equal) according to the user-defined criterion, otherwise, it is `.false.`.<br>
    !>                              The the input `array` is array-valued, then the last argument to `iseq` is the length of the input `pattern`.<br>
    !>                              The following illustrates the generic interface of `iseq` where `pattern` is array-valued,
    !>                              \code{.F90}
    !>                                  function iseq(Segment, pattern, lenPattern) result(equivalent)
    !>                                      use pm_kind, only: IK, LK
    !>                                      integer(IK) , intent(in)    :: lenPattern
    !>                                      TYPE(KIND)  , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                      logical(LK)                 :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
    !>                              \code{.F90}
    !>                                  character(*, SK), intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  integer(IK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  logical(LK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  complex(CK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  real(RK)        , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                              \endcode
    !>                              where the kinds `SKG`, `IKG`, `LKG`, `CKG`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              The following illustrates the generic interface of `iseq` where `pattern` is scalar-valued (**including Fortran scalar strings**),
    !>                              \code{.F90}
    !>                                  function iseq(segment, pattern) result(equivalent)
    !>                                      use pm_kind, only: LK
    !>                                      TYPE(KIND)  , intent(in)    :: segment, pattern
    !>                                      logical(LK)                 :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
    !>                              \code{.F90}
    !>                                  use pm_kind, only: SK, IK, LK, CK, RK
    !>                                  character(*, SK), intent(in)    :: segment, pattern
    !>                                  integer(IK)     , intent(in)    :: segment, pattern
    !>                                  logical(LK)     , intent(in)    :: segment, pattern
    !>                                  complex(CK)     , intent(in)    :: segment, pattern
    !>                                  real(RK)        , intent(in)    :: segment, pattern
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              This user-defined equivalence check is extremely useful where an equivalence test other than exact identity is needed,
    !>                              for example, when the array segments should match the input `pattern` only within a given threshold or,
    !>                              when the case-sensitivity in character comparisons do not matter.<br>
    !>                              In such cases, user can define a custom equivalence criterion within the user-defined external function `iseq` to achieve the goal.<br>
    !>                              (**optional**, the default equivalence operator is `.eqv.` if the input `array` is `logical`, otherwise `==`.)
    !>
    !>  \return
    !>  `index`                 :   The output scalar `integer` of default kind \IK representing the index of the first element of
    !>                              the input `array` such that `array(1:index-1)` contains only full repetitions of the user-specified `pattern`.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_arrayStrip, only: getSIL
    !>      use pm_kind, only: IK
    !>      integer(IK) :: index
    !>
    !>      index = getSIL(array, pattern) ! scalar strings
    !>      index = getSIL(array(:), pattern) ! all intrinsic types
    !>      index = getSIL(array(:), pattern(:)) ! all intrinsic types
    !>      index = getSIL(array, pattern, iseq) ! scalar strings
    !>      index = getSIL(array(:), pattern, iseq) ! all intrinsic types
    !>      index = getSIL(array(:), pattern(:), iseq) ! all intrinsic types
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The procedures under this generic interface are `impure` when the user-specified `external` procedure `iseq` is specified as input argument.
    !>
    !>  \warning
    !>  Note that in Fortran, trailing blanks are ignored in character equivalence checks, that is, `"Fortran" == "Fortran "` yields `.true.`.
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getSIL](@ref pm_arrayStrip::getSIL)<br>
    !>  [getSIR](@ref pm_arrayStrip::getSIR)<br>
    !>  [getStripped](@ref pm_arrayStrip::getStripped)<br>
    !>  [setReplaced](@ref pm_arrayReplace::setReplaced)<br>
    !>  [getReplaced](@ref pm_arrayReplace::getReplaced)<br>
    !>  [setInserted](@ref pm_arrayInsert::setInserted)<br>
    !>  [setSplit](@ref pm_arraySplit::setSplit)<br>
    !>  [setLoc](@ref pm_arrayFind::setLoc)<br>
    !>  [getLoc](@ref pm_arrayFind::getLoc)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_arrayStrip/getSIL/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_arrayStrip/getSIL/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayStrip](@ref test_pm_arrayStrip)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getSIL

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getSILDefCom_D0_D0_SK5(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if SK4_ENABLED
    PURE module function getSILDefCom_D0_D0_SK4(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if SK3_ENABLED
    PURE module function getSILDefCom_D0_D0_SK3(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if SK2_ENABLED
    PURE module function getSILDefCom_D0_D0_SK2(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if SK1_ENABLED
    PURE module function getSILDefCom_D0_D0_SK1(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getSILDefCom_D1_D0_SK5(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if SK4_ENABLED
    PURE module function getSILDefCom_D1_D0_SK4(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if SK3_ENABLED
    PURE module function getSILDefCom_D1_D0_SK3(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if SK2_ENABLED
    PURE module function getSILDefCom_D1_D0_SK2(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if SK1_ENABLED
    PURE module function getSILDefCom_D1_D0_SK1(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getSILDefCom_D1_D0_IK5(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if IK4_ENABLED
    PURE module function getSILDefCom_D1_D0_IK4(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if IK3_ENABLED
    PURE module function getSILDefCom_D1_D0_IK3(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if IK2_ENABLED
    PURE module function getSILDefCom_D1_D0_IK2(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if IK1_ENABLED
    PURE module function getSILDefCom_D1_D0_IK1(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getSILDefCom_D1_D0_LK5(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if LK4_ENABLED
    PURE module function getSILDefCom_D1_D0_LK4(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if LK3_ENABLED
    PURE module function getSILDefCom_D1_D0_LK3(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if LK2_ENABLED
    PURE module function getSILDefCom_D1_D0_LK2(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if LK1_ENABLED
    PURE module function getSILDefCom_D1_D0_LK1(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getSILDefCom_D1_D0_CK5(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if CK4_ENABLED
    PURE module function getSILDefCom_D1_D0_CK4(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if CK3_ENABLED
    PURE module function getSILDefCom_D1_D0_CK3(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if CK2_ENABLED
    PURE module function getSILDefCom_D1_D0_CK2(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if CK1_ENABLED
    PURE module function getSILDefCom_D1_D0_CK1(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getSILDefCom_D1_D0_RK5(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if RK4_ENABLED
    PURE module function getSILDefCom_D1_D0_RK4(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if RK3_ENABLED
    PURE module function getSILDefCom_D1_D0_RK3(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if RK2_ENABLED
    PURE module function getSILDefCom_D1_D0_RK2(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if RK1_ENABLED
    PURE module function getSILDefCom_D1_D0_RK1(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getSILCusCom_D0_D0_SK5(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK4_ENABLED
    module function getSILCusCom_D0_D0_SK4(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK3_ENABLED
    module function getSILCusCom_D0_D0_SK3(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK2_ENABLED
    module function getSILCusCom_D0_D0_SK2(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK1_ENABLED
    module function getSILCusCom_D0_D0_SK1(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getSILCusCom_D1_D0_SK5(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK4_ENABLED
    module function getSILCusCom_D1_D0_SK4(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK3_ENABLED
    module function getSILCusCom_D1_D0_SK3(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK2_ENABLED
    module function getSILCusCom_D1_D0_SK2(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK1_ENABLED
    module function getSILCusCom_D1_D0_SK1(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getSILCusCom_D1_D0_IK5(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if IK4_ENABLED
    module function getSILCusCom_D1_D0_IK4(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if IK3_ENABLED
    module function getSILCusCom_D1_D0_IK3(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if IK2_ENABLED
    module function getSILCusCom_D1_D0_IK2(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if IK1_ENABLED
    module function getSILCusCom_D1_D0_IK1(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getSILCusCom_D1_D0_LK5(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if LK4_ENABLED
    module function getSILCusCom_D1_D0_LK4(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if LK3_ENABLED
    module function getSILCusCom_D1_D0_LK3(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if LK2_ENABLED
    module function getSILCusCom_D1_D0_LK2(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if LK1_ENABLED
    module function getSILCusCom_D1_D0_LK1(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getSILCusCom_D1_D0_CK5(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if CK4_ENABLED
    module function getSILCusCom_D1_D0_CK4(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if CK3_ENABLED
    module function getSILCusCom_D1_D0_CK3(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if CK2_ENABLED
    module function getSILCusCom_D1_D0_CK2(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if CK1_ENABLED
    module function getSILCusCom_D1_D0_CK1(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getSILCusCom_D1_D0_RK5(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if RK4_ENABLED
    module function getSILCusCom_D1_D0_RK4(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if RK3_ENABLED
    module function getSILCusCom_D1_D0_RK3(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if RK2_ENABLED
    module function getSILCusCom_D1_D0_RK2(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if RK1_ENABLED
    module function getSILCusCom_D1_D0_RK1(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getSILDefCom_D1_D1_SK5(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if SK4_ENABLED
    PURE module function getSILDefCom_D1_D1_SK4(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if SK3_ENABLED
    PURE module function getSILDefCom_D1_D1_SK3(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if SK2_ENABLED
    PURE module function getSILDefCom_D1_D1_SK2(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if SK1_ENABLED
    PURE module function getSILDefCom_D1_D1_SK1(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getSILDefCom_D1_D1_IK5(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if IK4_ENABLED
    PURE module function getSILDefCom_D1_D1_IK4(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if IK3_ENABLED
    PURE module function getSILDefCom_D1_D1_IK3(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if IK2_ENABLED
    PURE module function getSILDefCom_D1_D1_IK2(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if IK1_ENABLED
    PURE module function getSILDefCom_D1_D1_IK1(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getSILDefCom_D1_D1_LK5(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if LK4_ENABLED
    PURE module function getSILDefCom_D1_D1_LK4(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if LK3_ENABLED
    PURE module function getSILDefCom_D1_D1_LK3(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if LK2_ENABLED
    PURE module function getSILDefCom_D1_D1_LK2(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if LK1_ENABLED
    PURE module function getSILDefCom_D1_D1_LK1(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getSILDefCom_D1_D1_CK5(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if CK4_ENABLED
    PURE module function getSILDefCom_D1_D1_CK4(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if CK3_ENABLED
    PURE module function getSILDefCom_D1_D1_CK3(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if CK2_ENABLED
    PURE module function getSILDefCom_D1_D1_CK2(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if CK1_ENABLED
    PURE module function getSILDefCom_D1_D1_CK1(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getSILDefCom_D1_D1_RK5(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if RK4_ENABLED
    PURE module function getSILDefCom_D1_D1_RK4(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if RK3_ENABLED
    PURE module function getSILDefCom_D1_D1_RK3(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if RK2_ENABLED
    PURE module function getSILDefCom_D1_D1_RK2(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if RK1_ENABLED
    PURE module function getSILDefCom_D1_D1_RK1(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILDefCom_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getSILCusCom_D1_D1_SK5(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK4_ENABLED
    module function getSILCusCom_D1_D1_SK4(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK3_ENABLED
    module function getSILCusCom_D1_D1_SK3(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK2_ENABLED
    module function getSILCusCom_D1_D1_SK2(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK1_ENABLED
    module function getSILCusCom_D1_D1_SK1(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getSILCusCom_D1_D1_IK5(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if IK4_ENABLED
    module function getSILCusCom_D1_D1_IK4(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if IK3_ENABLED
    module function getSILCusCom_D1_D1_IK3(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if IK2_ENABLED
    module function getSILCusCom_D1_D1_IK2(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if IK1_ENABLED
    module function getSILCusCom_D1_D1_IK1(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getSILCusCom_D1_D1_LK5(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if LK4_ENABLED
    module function getSILCusCom_D1_D1_LK4(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if LK3_ENABLED
    module function getSILCusCom_D1_D1_LK3(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if LK2_ENABLED
    module function getSILCusCom_D1_D1_LK2(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if LK1_ENABLED
    module function getSILCusCom_D1_D1_LK1(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getSILCusCom_D1_D1_CK5(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if CK4_ENABLED
    module function getSILCusCom_D1_D1_CK4(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if CK3_ENABLED
    module function getSILCusCom_D1_D1_CK3(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if CK2_ENABLED
    module function getSILCusCom_D1_D1_CK2(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if CK1_ENABLED
    module function getSILCusCom_D1_D1_CK1(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getSILCusCom_D1_D1_RK5(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if RK4_ENABLED
    module function getSILCusCom_D1_D1_RK4(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if RK3_ENABLED
    module function getSILCusCom_D1_D1_RK3(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if RK2_ENABLED
    module function getSILCusCom_D1_D1_RK2(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if RK1_ENABLED
    module function getSILCusCom_D1_D1_RK1(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSILCusCom_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `index`, the index of the last element of
    !>  the input `array` such that `array(index+1:)` contains only full repetitions of the user-specified `pattern`.
    !>
    !>  \param[in]  array       :   The input `contiguous` vector of either <br>
    !>                              <ul>
    !>                                  <li>    type `character` of kind \SKALL, or <br>
    !>                                  <li>    type `logical` of kind \LKALL, or <br>
    !>                                  <li>    type `integer` of kind \IKALL, or <br>
    !>                                  <li>    type `complex` of kind \CKALL, or <br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ul>
    !>                              or,
    !>                              <ul>
    !>                                  <li>    a scalar assumed-length `character` of kind \SKALL, <br>
    !>                              </ul>
    !>  \param[in]  pattern     :   The input object of the same type and kind as the input array, of rank similar to or lower than that of `array`, whose value is to be stripped from the beginning of `array`.
    !>  \param      iseq        :   The `external` user-specified function that takes either two input assumed-length `character` arguments
    !>                              (if the input `array` is also an assumed-length `character`) or two array-valued **explicit-shape**
    !>                              arguments of the same type and kind as the input `array`.<br>
    !>                              It must return a scalar of type `logical` of default kind \LK that is `.true.` if all elements of
    !>                              the two input arguments are equivalent (e.g., equal) according to the user-defined criterion, otherwise, it is `.false.`.<br>
    !>                              The the input `array` is array-valued, then the last argument to `iseq` is the length of the input `pattern`.<br>
    !>                              The following illustrates the generic interface of `iseq` where `pattern` is array-valued,<br>
    !>                              \code{.F90}
    !>                                  function iseq(Segment, pattern, lenPattern) result(equivalent)
    !>                                      use pm_kind, only: IK, LK
    !>                                      integer(IK) , intent(in)    :: lenPattern
    !>                                      TYPE(KIND)  , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                      logical(LK)                 :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,<br>
    !>                              \code{.F90}
    !>                                  character(*, SK), intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  integer(IK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  logical(LK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  complex(CK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  real(RK)        , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                              \endcode
    !>                              where the kinds `SKG`, `IKG`, `LKG`, `CKG`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              The following illustrates the generic interface of `iseq` where `pattern` is scalar-valued (**including Fortran scalar strings**),<br>
    !>                              \code{.F90}
    !>                                  function iseq(segment, pattern) result(equivalent)
    !>                                      use pm_kind, only: LK
    !>                                      TYPE(KIND)  , intent(in)    :: segment, pattern
    !>                                      logical(LK)                 :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,<br>
    !>                              \code{.F90}
    !>                                  use pm_kind, only: SK, IK, LK, CK, RK
    !>                                  character(*, SK), intent(in)    :: segment, pattern
    !>                                  integer(IK)     , intent(in)    :: segment, pattern
    !>                                  logical(LK)     , intent(in)    :: segment, pattern
    !>                                  complex(CK)     , intent(in)    :: segment, pattern
    !>                                  real(RK)        , intent(in)    :: segment, pattern
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              This user-defined equivalence check is extremely useful where an equivalence test other than exact identity is needed,
    !>                              for example, when the array segments should match the input `pattern` only within a given threshold or,
    !>                              when the case-sensitivity in character comparisons do not matter.<br>
    !>                              In such cases, user can define a custom equivalence criterion within the user-defined external function `iseq` to achieve the goal.<br>
    !>                              (**optional**, the default equivalence operator is `.eqv.` if the input `array` is `logical`, otherwise `==`)
    !>
    !>  \return
    !>  `index`                 :   The output scalar `integer` of default kind \IK representing the index of the last element of
    !>                              the input `array` such that `array(index+1:)` contains only full repetitions of the user-specified `pattern` .
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_arrayStrip, only: getSIR
    !>      use pm_kind, only: IK
    !>      integer(IK) :: index
    !>
    !>      index = getSIR(array, pattern) ! scalar strings
    !>      index = getSIR(array(:), pattern) ! all intrinsic types
    !>      index = getSIR(array(:), pattern(:)) ! all intrinsic types
    !>      index = getSIR(array, pattern, iseq) ! scalar strings
    !>      index = getSIR(array(:), pattern, iseq) ! all intrinsic types
    !>      index = getSIR(array(:), pattern(:), iseq) ! all intrinsic types
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The procedures under this generic interface are `impure` when the user-specified `external` procedure `iseq` is specified as input argument.
    !>
    !>  \warning
    !>  Note that in Fortran, trailing blanks are ignored in character equivalence checks, that is, `"Fortran" == "Fortran "` yields `.true.`.
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getSIL](@ref pm_arrayStrip::getSIL)<br>
    !>  [getSIR](@ref pm_arrayStrip::getSIR)<br>
    !>  [getStripped](@ref pm_arrayStrip::getStripped)<br>
    !>  [setReplaced](@ref pm_arrayReplace::setReplaced)<br>
    !>  [getReplaced](@ref pm_arrayReplace::getReplaced)<br>
    !>  [setInserted](@ref pm_arrayInsert::setInserted)<br>
    !>  [setSplit](@ref pm_arraySplit::setSplit)<br>
    !>  [setLoc](@ref pm_arrayFind::setLoc)<br>
    !>  [getLoc](@ref pm_arrayFind::getLoc)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_arrayStrip/getSIR/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_arrayStrip/getSIR/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayStrip](@ref test_pm_arrayStrip)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getSIR

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getSIRDefCom_D0_D0_SK5(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if SK4_ENABLED
    PURE module function getSIRDefCom_D0_D0_SK4(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if SK3_ENABLED
    PURE module function getSIRDefCom_D0_D0_SK3(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if SK2_ENABLED
    PURE module function getSIRDefCom_D0_D0_SK2(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if SK1_ENABLED
    PURE module function getSIRDefCom_D0_D0_SK1(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getSIRDefCom_D1_D0_SK5(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if SK4_ENABLED
    PURE module function getSIRDefCom_D1_D0_SK4(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if SK3_ENABLED
    PURE module function getSIRDefCom_D1_D0_SK3(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if SK2_ENABLED
    PURE module function getSIRDefCom_D1_D0_SK2(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if SK1_ENABLED
    PURE module function getSIRDefCom_D1_D0_SK1(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getSIRDefCom_D1_D0_IK5(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if IK4_ENABLED
    PURE module function getSIRDefCom_D1_D0_IK4(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if IK3_ENABLED
    PURE module function getSIRDefCom_D1_D0_IK3(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if IK2_ENABLED
    PURE module function getSIRDefCom_D1_D0_IK2(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if IK1_ENABLED
    PURE module function getSIRDefCom_D1_D0_IK1(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getSIRDefCom_D1_D0_LK5(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if LK4_ENABLED
    PURE module function getSIRDefCom_D1_D0_LK4(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if LK3_ENABLED
    PURE module function getSIRDefCom_D1_D0_LK3(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if LK2_ENABLED
    PURE module function getSIRDefCom_D1_D0_LK2(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if LK1_ENABLED
    PURE module function getSIRDefCom_D1_D0_LK1(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getSIRDefCom_D1_D0_CK5(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if CK4_ENABLED
    PURE module function getSIRDefCom_D1_D0_CK4(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if CK3_ENABLED
    PURE module function getSIRDefCom_D1_D0_CK3(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if CK2_ENABLED
    PURE module function getSIRDefCom_D1_D0_CK2(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if CK1_ENABLED
    PURE module function getSIRDefCom_D1_D0_CK1(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getSIRDefCom_D1_D0_RK5(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if RK4_ENABLED
    PURE module function getSIRDefCom_D1_D0_RK4(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if RK3_ENABLED
    PURE module function getSIRDefCom_D1_D0_RK3(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if RK2_ENABLED
    PURE module function getSIRDefCom_D1_D0_RK2(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

#if RK1_ENABLED
    PURE module function getSIRDefCom_D1_D0_RK1(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getSIRCusCom_D0_D0_SK5(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK4_ENABLED
    module function getSIRCusCom_D0_D0_SK4(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK3_ENABLED
    module function getSIRCusCom_D0_D0_SK3(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK2_ENABLED
    module function getSIRCusCom_D0_D0_SK2(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK1_ENABLED
    module function getSIRCusCom_D0_D0_SK1(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getSIRCusCom_D1_D0_SK5(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK4_ENABLED
    module function getSIRCusCom_D1_D0_SK4(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK3_ENABLED
    module function getSIRCusCom_D1_D0_SK3(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK2_ENABLED
    module function getSIRCusCom_D1_D0_SK2(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK1_ENABLED
    module function getSIRCusCom_D1_D0_SK1(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getSIRCusCom_D1_D0_IK5(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if IK4_ENABLED
    module function getSIRCusCom_D1_D0_IK4(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if IK3_ENABLED
    module function getSIRCusCom_D1_D0_IK3(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if IK2_ENABLED
    module function getSIRCusCom_D1_D0_IK2(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if IK1_ENABLED
    module function getSIRCusCom_D1_D0_IK1(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getSIRCusCom_D1_D0_LK5(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if LK4_ENABLED
    module function getSIRCusCom_D1_D0_LK4(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if LK3_ENABLED
    module function getSIRCusCom_D1_D0_LK3(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if LK2_ENABLED
    module function getSIRCusCom_D1_D0_LK2(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if LK1_ENABLED
    module function getSIRCusCom_D1_D0_LK1(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getSIRCusCom_D1_D0_CK5(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if CK4_ENABLED
    module function getSIRCusCom_D1_D0_CK4(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if CK3_ENABLED
    module function getSIRCusCom_D1_D0_CK3(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if CK2_ENABLED
    module function getSIRCusCom_D1_D0_CK2(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if CK1_ENABLED
    module function getSIRCusCom_D1_D0_CK1(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getSIRCusCom_D1_D0_RK5(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if RK4_ENABLED
    module function getSIRCusCom_D1_D0_RK4(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if RK3_ENABLED
    module function getSIRCusCom_D1_D0_RK3(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if RK2_ENABLED
    module function getSIRCusCom_D1_D0_RK2(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if RK1_ENABLED
    module function getSIRCusCom_D1_D0_RK1(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getSIRDefCom_D1_D1_SK5(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if SK4_ENABLED
    PURE module function getSIRDefCom_D1_D1_SK4(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if SK3_ENABLED
    PURE module function getSIRDefCom_D1_D1_SK3(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if SK2_ENABLED
    PURE module function getSIRDefCom_D1_D1_SK2(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if SK1_ENABLED
    PURE module function getSIRDefCom_D1_D1_SK1(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getSIRDefCom_D1_D1_IK5(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if IK4_ENABLED
    PURE module function getSIRDefCom_D1_D1_IK4(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if IK3_ENABLED
    PURE module function getSIRDefCom_D1_D1_IK3(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if IK2_ENABLED
    PURE module function getSIRDefCom_D1_D1_IK2(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if IK1_ENABLED
    PURE module function getSIRDefCom_D1_D1_IK1(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getSIRDefCom_D1_D1_LK5(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if LK4_ENABLED
    PURE module function getSIRDefCom_D1_D1_LK4(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if LK3_ENABLED
    PURE module function getSIRDefCom_D1_D1_LK3(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if LK2_ENABLED
    PURE module function getSIRDefCom_D1_D1_LK2(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if LK1_ENABLED
    PURE module function getSIRDefCom_D1_D1_LK1(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getSIRDefCom_D1_D1_CK5(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if CK4_ENABLED
    PURE module function getSIRDefCom_D1_D1_CK4(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if CK3_ENABLED
    PURE module function getSIRDefCom_D1_D1_CK3(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if CK2_ENABLED
    PURE module function getSIRDefCom_D1_D1_CK2(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if CK1_ENABLED
    PURE module function getSIRDefCom_D1_D1_CK1(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getSIRDefCom_D1_D1_RK5(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if RK4_ENABLED
    PURE module function getSIRDefCom_D1_D1_RK4(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if RK3_ENABLED
    PURE module function getSIRDefCom_D1_D1_RK3(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if RK2_ENABLED
    PURE module function getSIRDefCom_D1_D1_RK2(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

#if RK1_ENABLED
    PURE module function getSIRDefCom_D1_D1_RK1(array, pattern) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRDefCom_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getSIRCusCom_D1_D1_SK5(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK4_ENABLED
    module function getSIRCusCom_D1_D1_SK4(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK3_ENABLED
    module function getSIRCusCom_D1_D1_SK3(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK2_ENABLED
    module function getSIRCusCom_D1_D1_SK2(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if SK1_ENABLED
    module function getSIRCusCom_D1_D1_SK1(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getSIRCusCom_D1_D1_IK5(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if IK4_ENABLED
    module function getSIRCusCom_D1_D1_IK4(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if IK3_ENABLED
    module function getSIRCusCom_D1_D1_IK3(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if IK2_ENABLED
    module function getSIRCusCom_D1_D1_IK2(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if IK1_ENABLED
    module function getSIRCusCom_D1_D1_IK1(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getSIRCusCom_D1_D1_LK5(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if LK4_ENABLED
    module function getSIRCusCom_D1_D1_LK4(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if LK3_ENABLED
    module function getSIRCusCom_D1_D1_LK3(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if LK2_ENABLED
    module function getSIRCusCom_D1_D1_LK2(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if LK1_ENABLED
    module function getSIRCusCom_D1_D1_LK1(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getSIRCusCom_D1_D1_CK5(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if CK4_ENABLED
    module function getSIRCusCom_D1_D1_CK4(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if CK3_ENABLED
    module function getSIRCusCom_D1_D1_CK3(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if CK2_ENABLED
    module function getSIRCusCom_D1_D1_CK2(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if CK1_ENABLED
    module function getSIRCusCom_D1_D1_CK1(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getSIRCusCom_D1_D1_RK5(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if RK4_ENABLED
    module function getSIRCusCom_D1_D1_RK4(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if RK3_ENABLED
    module function getSIRCusCom_D1_D1_RK3(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if RK2_ENABLED
    module function getSIRCusCom_D1_D1_RK2(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

#if RK1_ENABLED
    module function getSIRCusCom_D1_D1_RK1(array, pattern, iseq) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSIRCusCom_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)                                         :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayStrip ! LCOV_EXCL_LINE