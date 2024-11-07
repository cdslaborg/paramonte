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
!>  This module contains procedures and generic interfaces for remapping arrays of various types.
!>
!>  \benchmarks
!>
!>  \benchmark{setRemapped_getRemapped_direct, The runtime performance of [setRemapped](@ref pm_arrayRemap::setRemapped) vs. direct remapping}
!>  \include{lineno} benchmark/pm_arrayRemap/setRemapped_getRemapped_direct/main.F90
!>  \compilefb{setRemapped_getRemapped_direct}
!>  \postprocb{setRemapped_getRemapped_direct}
!>  \include{lineno} benchmark/pm_arrayRemap/setRemapped_getRemapped_direct/main.py
!>  \visb{setRemapped_getRemapped_direct}
!>  \image html benchmark/pm_arrayRemap/setRemapped_getRemapped_direct/benchmark.setRemapped_getRemapped_direct.runtime.png width=1000
!>  \image html benchmark/pm_arrayRemap/setRemapped_getRemapped_direct/benchmark.setRemapped_getRemapped_direct.runtime.ratio.png width=1000
!>  \moralb{setRemapped_getRemapped_direct}
!>      -#  The procedures under the generic interface [setRemapped](@ref pm_arrayRemap::setRemapped) tend to be significantly faster
!>          than directly remapping arrays. This likely only true for remapping of **allocatable** arrays.<br>
!>          The primary reason for the better performance of [setRemapped](@ref pm_arrayRemap::setRemapped) is that
!>          [setRemapped](@ref pm_arrayRemap::setRemapped) avoids a final data copy from a dummy array to the original
!>          array by copying the allocation descriptor instead of the whole remapped array to the original array.<br>
!>      -#  Note that the observed performance benefit slightly diminishes if the remapping is not in `action =` [reverse](@ref pm_array::reverse) mode.<br>
!>      -#  The other benefit of [setRemapped](@ref pm_arrayRemap::setRemapped) is that it provides a unified seamless generic interface
!>          for reversing all intrinsic types and kinds of arrays as well as assumed-length and assumed-shape characters.<br>
!>
!>  \test
!>  [test_pm_arrayRemap](@ref test_pm_arrayRemap)
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX
module pm_arrayRemap

    use pm_kind, only: SK, IK, LK
    use pm_array, only: reverse_type, reverse

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_arrayRemap"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate a copy of the input `array` whose elements are reordered according to the input `index` array such that,<br>
    !>  &nbsp; `arrayNew = array(index)` or, <br>
    !>  &nbsp; `arrayNew = array(index(size(index):1:-1))` if `action =` [reverse](@ref pm_array::reverse) or, <br>
    !>  holds for the output.
    !>
    !>  \param[in]  array       :   The **input contiguous** array of shape `(:)` of either <br>
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
    !>                              whose elements will be reordered according to the input `index`.<br>
    !>                              If the input argument `arrayNew` is specified, then `array` has `intent(in)` and will not change.<br>
    !>                              If the input argument `arrayNew` is missing, then `array` will be reallocated and reordered on output.<br>
    !>                              The the `allocatable` attribute and the reallocation of `array` is essential for efficient fast reordering of `array`.<br>
    !>  \param[in]  index       :   The input `contiguous` array of shape `(:)` of type `integer` of default kind \IK of the same size as `array`, containing the reordering indices.<br>
    !>  \param[in]  action      :   The input scalar object that can be
    !>                              <ol>
    !>                                  <li>    the constant [reverse](@ref pm_array::reverse) or equivalently, an object of type [reverse_type](@ref pm_array::reverse_type).<br>
    !>                                          Specifying this value implies that the input `index` must be reversed (`index(size(index):1:-1)`) prior to being used for remapping.<br>
    !>                                          On output, the input `array` will be reordered to `array(index(size(index):1:-1))`.<br>
    !>                              </ol>
    !>                              (**optional**, If missing, the input `index` will be used as is, such that `array` will be reordered to `array(index)` on output.)
    !>
    !>  \return
    !>  `arrayNew`              :   The output array of the same type, kind, shape, and size as the input `array` that will contain the reordered array
    !>
    !>  \interface{getRemapped}
    !>  \code{.F90}
    !>
    !>      use pm_arrayRemap, only: getRemapped
    !>
    !>      arrayNew = getRemapped(array, index)            ! `array` must be allocatable.
    !>      arrayNew = getRemapped(array, index, action)    ! `array` must be allocatable.
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \warning
    !>  The sizes of `array`, `arrayNew`, and `index` arguments must be equal.<br>
    !>  All elements of `index` must be within the lower and upper bounds of `array`.<br>
    !>  \vericons
    !>
    !>  \remark
    !>  The primary purpose of the procedures under this generic interface is to provide a convenient bug-free method of generating a remapped copy
    !>  of an array of any intrinsic type and kind, similar to the what is implicitly done by current compilers with the regular Fortran syntax.<br>
    !>  With further compiler and language template enhancements in the future, the need for the procedures under this generic interface might be resolved in the future.<br>
    !>  See [pm_arrayRemap](@ref pm_arrayRemap) for the relevant benchmarks.<br>
    !>
    !>  \see
    !>  [setRemapped](@ref pm_arrayRemap::setRemapped)<br>
    !>  [setShuffled](@ref pm_arrayShuffle::setShuffled)<br>
    !>  [getReversed](@ref pm_arrayReverse::getReversed)<br>
    !>  [setReversed](@ref pm_arrayReverse::setReversed)<br>
    !>
    !>  \example{getRemapped}
    !>  \include{lineno} example/pm_arrayRemap/getRemapped/main.F90
    !>  \compilef{getRemapped}
    !>  \output{getRemapped}
    !>  \include{lineno} example/pm_arrayRemap/getRemapped/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayRemap](@ref test_pm_arrayRemap)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to 1D input objects.<br>
    !>
    !>  \todo
    !>  \pvhigh The gfortran bugs in the implementations of this generic interface must be resolved in the future.<br>
    !>
    !>  \final{getRemapped}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getRemapped

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRemappedFor_D0_SK5(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        character(len(array,IK),SKG)                            :: arrayNew
    end function
#endif

#if SK4_ENABLED
    PURE module function getRemappedFor_D0_SK4(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        character(len(array,IK),SKG)                            :: arrayNew
    end function
#endif

#if SK3_ENABLED
    PURE module function getRemappedFor_D0_SK3(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        character(len(array,IK),SKG)                            :: arrayNew
    end function
#endif

#if SK2_ENABLED
    PURE module function getRemappedFor_D0_SK2(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        character(len(array,IK),SKG)                            :: arrayNew
    end function
#endif

#if SK1_ENABLED
    PURE module function getRemappedFor_D0_SK1(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        character(len(array,IK),SKG)                            :: arrayNew
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRemappedFor_D1_SK5(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        character(len(array,IK),SKG)                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if SK4_ENABLED
    PURE module function getRemappedFor_D1_SK4(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        character(len(array,IK),SKG)                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if SK3_ENABLED
    PURE module function getRemappedFor_D1_SK3(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        character(len(array,IK),SKG)                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if SK2_ENABLED
    PURE module function getRemappedFor_D1_SK2(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        character(len(array,IK),SKG)                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if SK1_ENABLED
    PURE module function getRemappedFor_D1_SK1(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        character(len(array,IK),SKG)                            :: arrayNew(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getRemappedFor_D1_IK5(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        integer(IKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if IK4_ENABLED
    PURE module function getRemappedFor_D1_IK4(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        integer(IKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if IK3_ENABLED
    PURE module function getRemappedFor_D1_IK3(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        integer(IKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if IK2_ENABLED
    PURE module function getRemappedFor_D1_IK2(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        integer(IKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if IK1_ENABLED
    PURE module function getRemappedFor_D1_IK1(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        integer(IKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getRemappedFor_D1_LK5(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        logical(LKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if LK4_ENABLED
    PURE module function getRemappedFor_D1_LK4(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        logical(LKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if LK3_ENABLED
    PURE module function getRemappedFor_D1_LK3(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        logical(LKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if LK2_ENABLED
    PURE module function getRemappedFor_D1_LK2(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        logical(LKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if LK1_ENABLED
    PURE module function getRemappedFor_D1_LK1(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        logical(LKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getRemappedFor_D1_CK5(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        complex(CKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getRemappedFor_D1_CK4(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        complex(CKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getRemappedFor_D1_CK3(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        complex(CKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getRemappedFor_D1_CK2(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        complex(CKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getRemappedFor_D1_CK1(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        complex(CKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getRemappedFor_D1_RK5(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        real(RKG)                                               :: arrayNew(size(array, kind = IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getRemappedFor_D1_RK4(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        real(RKG)                                               :: arrayNew(size(array, kind = IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getRemappedFor_D1_RK3(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        real(RKG)                                               :: arrayNew(size(array, kind = IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getRemappedFor_D1_RK2(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        real(RKG)                                               :: arrayNew(size(array, kind = IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getRemappedFor_D1_RK1(array, index) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedFor_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        real(RKG)                                               :: arrayNew(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRemappedRev_D0_SK5(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        character(len(array,IK),SKG)                            :: arrayNew
    end function
#endif

#if SK4_ENABLED
    PURE module function getRemappedRev_D0_SK4(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        character(len(array,IK),SKG)                            :: arrayNew
    end function
#endif

#if SK3_ENABLED
    PURE module function getRemappedRev_D0_SK3(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        character(len(array,IK),SKG)                            :: arrayNew
    end function
#endif

#if SK2_ENABLED
    PURE module function getRemappedRev_D0_SK2(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        character(len(array,IK),SKG)                            :: arrayNew
    end function
#endif

#if SK1_ENABLED
    PURE module function getRemappedRev_D0_SK1(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        character(len(array,IK),SKG)                            :: arrayNew
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRemappedRev_D1_SK5(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        character(len(array,IK),SKG)                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if SK4_ENABLED
    PURE module function getRemappedRev_D1_SK4(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        character(len(array,IK),SKG)                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if SK3_ENABLED
    PURE module function getRemappedRev_D1_SK3(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        character(len(array,IK),SKG)                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if SK2_ENABLED
    PURE module function getRemappedRev_D1_SK2(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        character(len(array,IK),SKG)                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if SK1_ENABLED
    PURE module function getRemappedRev_D1_SK1(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        character(len(array,IK),SKG)                            :: arrayNew(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getRemappedRev_D1_IK5(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        integer(IKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if IK4_ENABLED
    PURE module function getRemappedRev_D1_IK4(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        integer(IKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if IK3_ENABLED
    PURE module function getRemappedRev_D1_IK3(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        integer(IKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if IK2_ENABLED
    PURE module function getRemappedRev_D1_IK2(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        integer(IKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if IK1_ENABLED
    PURE module function getRemappedRev_D1_IK1(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        integer(IKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getRemappedRev_D1_LK5(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        logical(LKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if LK4_ENABLED
    PURE module function getRemappedRev_D1_LK4(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        logical(LKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if LK3_ENABLED
    PURE module function getRemappedRev_D1_LK3(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        logical(LKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if LK2_ENABLED
    PURE module function getRemappedRev_D1_LK2(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        logical(LKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if LK1_ENABLED
    PURE module function getRemappedRev_D1_LK1(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        logical(LKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getRemappedRev_D1_CK5(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        complex(CKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getRemappedRev_D1_CK4(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        complex(CKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getRemappedRev_D1_CK3(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        complex(CKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getRemappedRev_D1_CK2(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        complex(CKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getRemappedRev_D1_CK1(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        complex(CKG)                                            :: arrayNew(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getRemappedRev_D1_RK5(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        real(RKG)                                               :: arrayNew(size(array, kind = IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getRemappedRev_D1_RK4(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        real(RKG)                                               :: arrayNew(size(array, kind = IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getRemappedRev_D1_RK3(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        real(RKG)                                               :: arrayNew(size(array, kind = IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getRemappedRev_D1_RK2(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        real(RKG)                                               :: arrayNew(size(array, kind = IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getRemappedRev_D1_RK1(array, index, action) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemappedRev_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        real(RKG)                                               :: arrayNew(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Reorder the elements of the input `array` according to the input `index` array, such that <br>
    !>  &nbsp; `array = array(index)` or, <br>
    !>  &nbsp; `array = array(index(size(index):1:-1))` if `action =` [reverse](@ref pm_array::reverse) or, <br>
    !>  &nbsp; `arrayNew = array(index)` if `arrayNew` is specified as input argument or, <br>
    !>  &nbsp; `arrayNew = array(index(size(index):1:-1))` if `arrayNew` and `action =` [reverse](@ref pm_array::reverse) are specified as input arguments, <br>
    !>
    !>  \param[inout]   array       :   The **input contiguous** or **input/output allocatable** array of shape `(:)` of either <br>
    !>                                  <ul>
    !>                                      <li>    type `character` of kind \SKALL, or <br>
    !>                                      <li>    type `integer` of kind \IKALL, or <br>
    !>                                      <li>    type `logical` of kind \LKALL, or <br>
    !>                                      <li>    type `complex` of kind \CKALL, or <br>
    !>                                      <li>    type `real` of kind \RKALL, <br>
    !>                                  </ul>
    !>                                  or,
    !>                                  <ul>
    !>                                      <li>    a scalar assumed-length `character` of kind \SKALL, <br>
    !>                                  </ul>
    !>                                  whose elements will be reordered according to the input `index`.<br>
    !>                                  If the input argument `arrayNew` is specified, then `array` has `intent(in)` and will not change.<br>
    !>                                  If the input argument `arrayNew` is missing, then `array` will be reallocated and reordered on output.<br>
    !>                                  The the `allocatable` attribute and the reallocation of `array` is essential for efficient fast reordering of `array`.<br>
    !>  \param[in]      index       :   The input `contiguous` array of shape `(:)` of type `integer` of default kind \IK of the same size as `array`, containing the reordering indices.<br>
    !>  \param[in]      action      :   The input scalar object that can be
    !>                                  <ol>
    !>                                      <li>    the constant [reverse](@ref pm_array::reverse) or equivalently, an object of type [reverse_type](@ref pm_array::reverse_type).<br>
    !>                                              Specifying this value implies that the input `index` must be reversed (`index(size(index):1:-1)`) prior to being used for remapping.<br>
    !>                                              On output, the input `array` will be reordered to `array(index(size(index):1:-1))`.<br>
    !>                                  </ol>
    !>                                  (**optional**, If missing, the input `index` will be used as is, such that `array` will be reordered to `array(index)` on output.)
    !>  \param[out]     arrayNew    :   The output `contiguous` array of the same type, kind, shape, and size as the input `array` that will contain the reordered array.<br>
    !>                                  (**optional**, if missing, the ordered array will stored and output in the input allocatable `array` **in-place**)
    !>
    !>  \interface{setRemapped}
    !>  \code{.F90}
    !>
    !>      use pm_arrayRemap, only: setRemapped
    !>
    !>      call setRemapped(array, index)                      ! `array` must be allocatable.
    !>      call setRemapped(array, index, action)              ! `array` must be allocatable.
    !>      call setRemapped(array, index, arrayNew)            ! `array` can be any supported contiguous entity.
    !>      call setRemapped(array, index, action, arrayNew)    ! `array` can be any supported contiguous entity.
    !>      !
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \warning
    !>  The sizes of `array`, `arrayNew`, and `index` arguments must be equal.<br>
    !>  All elements of `index` must be within the lower and upper bounds of `array`.<br>
    !>  \vericons
    !>
    !>  \remark
    !>  The primary purpose of the procedures under this generic interface is to provide a more efficient faster method of remapping an
    !>  array of any intrinsic type and kind, in place, without an extra copy that is implicitly done by current compilers with the regular Fortran syntax.<br>
    !>  This generic interface also provides a comfortable generic way to reverse arrays of any intrinsic type and kind. This is particularly useful in case of scalar character.<br>
    !>  With further compiler and language template enhancements in the future, the need for the procedures under this generic interface might be resolved in the future.<br>
    !>  See [pm_arrayRemap](@ref pm_arrayRemap) for the relevant benchmarks.<br>
    !>
    !>  \note
    !>  Upon return, if `array` is an `allocatable` with `intent(inout)`, then it is guaranteed to have the same lower bound as before.
    !>  This happens when the output argument `arrayNew` is missing.
    !>
    !>  \see
    !>  [getRemapped](@ref pm_arrayRemap::getRemapped)<br>
    !>  [setShuffled](@ref pm_arrayShuffle::setShuffled)<br>
    !>  [getReversed](@ref pm_arrayReverse::getReversed)<br>
    !>  [setReversed](@ref pm_arrayReverse::setReversed)<br>
    !>
    !>  \example{setRemapped}
    !>  \include{lineno} example/pm_arrayRemap/setRemapped/main.F90
    !>  \compilef{setRemapped}
    !>  \output{setRemapped}
    !>  \include{lineno} example/pm_arrayRemap/setRemapped/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayRemap](@ref test_pm_arrayRemap)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{10.3}
    !>  \desc
    !>  There is still a bug in gfortran as of version 10.3, where the compiler fails to properly
    !>  deallocate the array of origin in a call to intrinsic `move_alloc`:
    !>  \code{.F90}
    !>      call move_alloc(from = arrayNew, to = array)
    !>  \endcode
    !>  \remedy
    !>  Currently, a preprocessor fence is added specifically for gfortran to explicitly deallocate `arrayNew` after the call to `move_alloc()`.<br>
    !>  This is going to yield an extra tiny performance penalty. Upon the necessary corrections to gfortran, this extra preprocessing fence must be removed.<br>
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{10.3-12}
    !>  \desc
    !>  \gfortran currently does not accept deferred-length allocatable characters as
    !>  actual argument to procedures with assumed-length allocatable dummy arguments.<br>
    !>  \remedy
    !>  Avoid passing deferred-length `allocatable` scalar `character` arguments in the interfaces of the library until this \gfortran bug is resolved.<br>
    !>
    !>  \bug
    !>  \status \unresolved, possibly \resolved in \ifort{2023}.
    !>  \source \ifort{2021.4}
    !>  \desc
    !>  The \ifort cannot compile interfaces with contiguous arguments whose lower bounds
    !>  depend on the lower bounds of `allocatable` input arguments.<br>
    !>  The following is an illustration of the bug,<br>
    !>  \code{.F90}
    !>      PURE module lower_bound
    !>          interface lower
    !>              PURE module subroutine test_lower(array, indices)
    !>                  integer, intent(inout), allocatable :: array(:)
    !>                  integer, intent(in), contiguous     :: indices(lbound(array,1):)
    !>              end subroutine test_lower
    !>          end interface
    !>      end module lower_bound
    !>  \endcode
    !>  \code{.F90}
    !>      submodule (lower_bound) smod
    !>      contains
    !>          PURE module procedure test_lower
    !>          end procedure
    !>      end submodule smod
    !>  \endcode
    !>  Compile with
    !>  \code{.sh}
    !>      ifort lower_bound.F90 lower_bound@smod.F90 -c
    !>  \endcode
    !>  to reproduce the compiler bug.
    !>  \remedy
    !>  For now, avoid such interfaces in the library.<br>
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to 1D input objects.
    !>
    !>  \todo
    !>  \pvhigh The gfortran bugs in the implementations of this generic interface must be resolved in the future.<br>
    !>
    !>  \final{setRemapped}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setRemapped

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRemappedForOld_D0_SK5(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(:,SKG)        , intent(inout) , allocatable   :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRemappedForOld_D0_SK4(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)        , intent(inout) , allocatable   :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRemappedForOld_D0_SK3(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)        , intent(inout) , allocatable   :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRemappedForOld_D0_SK2(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)        , intent(inout) , allocatable   :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRemappedForOld_D0_SK1(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)        , intent(inout) , allocatable   :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRemappedForOld_D1_SK5(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRemappedForOld_D1_SK4(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRemappedForOld_D1_SK3(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRemappedForOld_D1_SK2(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRemappedForOld_D1_SK1(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRemappedForOld_D1_IK5(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRemappedForOld_D1_IK4(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRemappedForOld_D1_IK3(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRemappedForOld_D1_IK2(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRemappedForOld_D1_IK1(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRemappedForOld_D1_LK5(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRemappedForOld_D1_LK4(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRemappedForOld_D1_LK3(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRemappedForOld_D1_LK2(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRemappedForOld_D1_LK1(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRemappedForOld_D1_CK5(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRemappedForOld_D1_CK4(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRemappedForOld_D1_CK3(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRemappedForOld_D1_CK2(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRemappedForOld_D1_CK1(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRemappedForOld_D1_RK5(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRemappedForOld_D1_RK4(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRemappedForOld_D1_RK3(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRemappedForOld_D1_RK2(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRemappedForOld_D1_RK1(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForOld_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRemappedRevOld_D0_SK5(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(:,SKG)        , intent(inout) , allocatable   :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRemappedRevOld_D0_SK4(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)        , intent(inout) , allocatable   :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRemappedRevOld_D0_SK3(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)        , intent(inout) , allocatable   :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRemappedRevOld_D0_SK2(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)        , intent(inout) , allocatable   :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRemappedRevOld_D0_SK1(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)        , intent(inout) , allocatable   :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRemappedRevOld_D1_SK5(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRemappedRevOld_D1_SK4(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRemappedRevOld_D1_SK3(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRemappedRevOld_D1_SK2(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRemappedRevOld_D1_SK1(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRemappedRevOld_D1_IK5(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRemappedRevOld_D1_IK4(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRemappedRevOld_D1_IK3(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRemappedRevOld_D1_IK2(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRemappedRevOld_D1_IK1(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRemappedRevOld_D1_LK5(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRemappedRevOld_D1_LK4(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRemappedRevOld_D1_LK3(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRemappedRevOld_D1_LK2(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRemappedRevOld_D1_LK1(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRemappedRevOld_D1_CK5(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRemappedRevOld_D1_CK4(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRemappedRevOld_D1_CK3(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRemappedRevOld_D1_CK2(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRemappedRevOld_D1_CK1(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRemappedRevOld_D1_RK5(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRemappedRevOld_D1_RK4(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRemappedRevOld_D1_RK3(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRemappedRevOld_D1_RK2(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRemappedRevOld_D1_RK1(array, index, action)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevOld_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRemappedForNew_D0_SK5(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        character(*,SKG)        , intent(out)                   :: arrayNew
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRemappedForNew_D0_SK4(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        character(*,SKG)        , intent(out)                   :: arrayNew
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRemappedForNew_D0_SK3(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        character(*,SKG)        , intent(out)                   :: arrayNew
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRemappedForNew_D0_SK2(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        character(*,SKG)        , intent(out)                   :: arrayNew
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRemappedForNew_D0_SK1(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        character(*,SKG)        , intent(out)                   :: arrayNew
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRemappedForNew_D1_SK5(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        character(*,SKG)        , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRemappedForNew_D1_SK4(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        character(*,SKG)        , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRemappedForNew_D1_SK3(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        character(*,SKG)        , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRemappedForNew_D1_SK2(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        character(*,SKG)        , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRemappedForNew_D1_SK1(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        character(*,SKG)        , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRemappedForNew_D1_IK5(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        integer(IKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRemappedForNew_D1_IK4(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        integer(IKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRemappedForNew_D1_IK3(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        integer(IKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRemappedForNew_D1_IK2(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        integer(IKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRemappedForNew_D1_IK1(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        integer(IKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRemappedForNew_D1_LK5(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        logical(LKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRemappedForNew_D1_LK4(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        logical(LKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRemappedForNew_D1_LK3(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        logical(LKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRemappedForNew_D1_LK2(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        logical(LKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRemappedForNew_D1_LK1(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        logical(LKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRemappedForNew_D1_CK5(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        complex(CKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRemappedForNew_D1_CK4(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        complex(CKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRemappedForNew_D1_CK3(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        complex(CKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRemappedForNew_D1_CK2(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        complex(CKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRemappedForNew_D1_CK1(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        complex(CKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRemappedForNew_D1_RK5(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        real(RKG)               , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRemappedForNew_D1_RK4(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        real(RKG)               , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRemappedForNew_D1_RK3(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        real(RKG)               , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRemappedForNew_D1_RK2(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        real(RKG)               , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRemappedForNew_D1_RK1(array, index, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedForNew_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        real(RKG)               , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRemappedRevNew_D0_SK5(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        character(*,SKG)        , intent(out)                   :: arrayNew
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRemappedRevNew_D0_SK4(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        character(*,SKG)        , intent(out)                   :: arrayNew
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRemappedRevNew_D0_SK3(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        character(*,SKG)        , intent(out)                   :: arrayNew
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRemappedRevNew_D0_SK2(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        character(*,SKG)        , intent(out)                   :: arrayNew
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRemappedRevNew_D0_SK1(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        character(*,SKG)        , intent(out)                   :: arrayNew
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRemappedRevNew_D1_SK5(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        character(*,SKG)        , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRemappedRevNew_D1_SK4(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        character(*,SKG)        , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRemappedRevNew_D1_SK3(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        character(*,SKG)        , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRemappedRevNew_D1_SK2(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        character(*,SKG)        , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRemappedRevNew_D1_SK1(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        character(*,SKG)        , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRemappedRevNew_D1_IK5(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        integer(IKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRemappedRevNew_D1_IK4(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        integer(IKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRemappedRevNew_D1_IK3(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        integer(IKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRemappedRevNew_D1_IK2(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        integer(IKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRemappedRevNew_D1_IK1(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        integer(IKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRemappedRevNew_D1_LK5(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        logical(LKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRemappedRevNew_D1_LK4(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        logical(LKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRemappedRevNew_D1_LK3(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        logical(LKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRemappedRevNew_D1_LK2(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        logical(LKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRemappedRevNew_D1_LK1(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        logical(LKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRemappedRevNew_D1_CK5(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        complex(CKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRemappedRevNew_D1_CK4(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        complex(CKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRemappedRevNew_D1_CK3(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        complex(CKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRemappedRevNew_D1_CK2(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        complex(CKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRemappedRevNew_D1_CK1(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        complex(CKG)            , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRemappedRevNew_D1_RK5(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        real(RKG)               , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRemappedRevNew_D1_RK4(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        real(RKG)               , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRemappedRevNew_D1_RK3(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        real(RKG)               , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRemappedRevNew_D1_RK2(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        real(RKG)               , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRemappedRevNew_D1_RK1(array, index, action, arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemappedRevNew_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        integer(IK)             , intent(in)    , contiguous    :: index(:)
        type(reverse_type)      , intent(in)                    :: action
        real(RKG)               , intent(out)   , contiguous    :: arrayNew(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayRemap ! LCOV_EXCL_LINE