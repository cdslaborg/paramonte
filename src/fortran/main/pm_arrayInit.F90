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
!>  This module contains procedures and generic interfaces for efficient initialization of arbitrary rectangular cores and
!>  surrounding halos of arrays of arbitrary size, shape, and rank of arbitrary intrinsic type and kind.<br>
!>
!>  \details
!>  The **halo of an array of arbitrary rank** is defined as the complement of rectangular (hollow) region of arbitrary size within a super-rectangle.<br>
!>  The following figures illustrate example scalar and array cores and halos of two 1D and 2D arrays.<br>
!>
!>  \htmlonly
!>      <img src="pm_arrayInit@scalar.png" style="width:50%;">
!>  \endhtmlonly
!>
!>  \htmlonly
!>      <img src="pm_arrayInit@array.png" style="width:50%;">
!>  \endhtmlonly
!>
!>  \note
!>  All procedures of this module initialize **all elements** of the output arrays.<br>
!>
!>  \see
!>  [pm_arrayFill](@ref pm_arrayFill)<br>
!>  [pm_matrixInit](@ref pm_matrixInit)<br>
!>
!>  \test
!>  [test_pm_arrayInit](@ref test_pm_arrayInit)
!>
!>  \todo
!>  \pmed
!>  A separate module `pm_arrayInitHalo` should be added in future to handle initialization of halos exclusively (without changing array cores).
!>
!>  \bug
!>  \status \unresolved
!>  \source \ifx{2023.0.0 20221201}
!>  \desc 
!>  \ifx yields an ICE while compiling this module.<br>
!>  By contrast, the Intel classic Fortran compiler `ifort` can successfully compile this file.<br>
!>  The origins of the ICE remain unexplored as of today.<br>
!>  \remedy
!>  Avoid compiling the library with \ifx.
!>  
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayInit

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_arrayInit"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return an array of specified rank and shape of arbitrary intrinsic type and kind
    !>  with its rectangular core and halo set to the corresponding user-specified values.<br>
    !>
    !>  \details
    !>  The **halo of an array of arbitrary rank** is defined as the complement of rectangular core region of arbitrary size within a super-rectangle array.<br>
    !>  The following figure illustrates example cores and halos of two 1D and 2D arrays.<br>
    !>
    !>  \htmlonly
    !>      <img src="pm_arrayInit@array.png" style="width:50%;">
    !>  \endhtmlonly
    !>
    !>  \param[in]      Size    :   The input non-negative scalar or vector of size `rank(array)` of type `integer` of default kind \IK,
    !>                              representing the shape of the output array, `array(@Size)`.
    !>                              <ul>
    !>                                  <li>    If the desired rank of the output `array` is `<= 1`, then `Size` must be a scalar.<br>
    !>                                  <li>    If the desired rank of the output `array` is `>= 2`, then `Size` must be a vector of size `rank(array)`.<br>
    !>                              </ul>
    !>  \param[in]      core    :   The input `contiguous` array of the same rank as the output `array`, of the same type and kind as the output `array`,
    !>                              containing the array value to paste into the rectangular core of `array`.<br>
    !>                              The extents of `core` (i.e., `shape(core)`) delineate the rectangular `core` of `array` into which `core` will be pasted.<br>
    !>                              <ul>
    !>                                  <li>    If `array` is of type `character`, then the condition `len(core) <= len(array)` must hold.<br>
    !>                              </ul>
    !>  \param[in]      halo    :   The input scalar of the same type and kind as the output `array`, containing the value to write to the halo surrounding the specified rectangular core of the `array`.<br>
    !>                              <ul>
    !>                                  <li>    If `array` is a scalar string, then `halo` must be a single `character` (with length type parameter `1`).<br>
    !>                                  <li>    If `array` is an array of type `character`, then `halo` can have arbitrary length type parameter satisfying the condition `len(halo) <= len(array)`.<br>
    !>                              </ul>
    !>  \param[in]      coffset :   The input non-negative scalar or vector of size `rank(array)` of type `integer` of default kind \IK representing the index offset from the first (top-left corner) element of `array`,
    !>                              such that `array(@coffset + 1)` represents the indices of first element of the rectangular core of `array`.<br>
    !>                              <ul>
    !>                                  <li>    If the output `array` is of rank `<= 1`, then `coffset` must be a scalar.<br>
    !>                                  <li>    If the output `array` is of rank `>= 2`, then `coffset` must be a vector of size `rank(array)`.<br>
    !>                              </ul>
    !>
    !>  \return
    !>  `array`                 :   The output scalar of,<br>
    !>                              <ul>
    !>                                  <li>    type `character` of kind \SKALL of the specified length type parameter `size` or,<br>
    !>                              </ul>
    !>                              the output array of rank `1`, `2`, or `3`, of the user-specified shape `size`, of either,<br>
    !>                              <ul>
    !>                                  <li>    type `character` of kind \SKALL, or<br>
    !>                                  <li>    type `integer` of kind \IKALL, or<br>
    !>                                  <li>    type `logical` of kind \LKALL, or<br>
    !>                                  <li>    type `complex` of kind \CKALL, or<br>
    !>                                  <li>    type `real` of kind \RKALL.<br>
    !>                              </ul>
    !>                              On output, the rectangular core and halo of the array will be initialized to the specified input arguments `core` and `halo` respectively.<br>
    !>
    !>  \interface{getCoreHalo}
    !>  \code{.F90}
    !>
    !>      use pm_arrayInit, only: getCoreHalo
    !>
    !>      call getCoreHalo(array, core, halo, coffset) ! `array` is scalar string. All other arguments are scalar.
    !>      call getCoreHalo(array(:), core(@shape(array)), halo, coffset) ! `array` and `core` are vectors. All other arguments are scalar.
    !>      call getCoreHalo(array(..), core(@shape(array)), halo, coffset(1:rank(array))) ! `array` and `core` are of arbitrary but identical ranks.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  All warnings and conditions associated with [setCoreHalo](@ref pm_arrayInit::setCoreHalo) also apply to this generic interface.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  This generic interface is merely a convenience wrapper around the generic interface [setCoreHalo](@ref pm_arrayInit::setCoreHalo).<br>
    !>
    !>  \see
    !>  [setCoreHalo](@ref pm_arrayInit::setCoreHalo)<br>
    !>  [setResized](@ref pm_arrayResize::setResized)<br>
    !>  [setRebound](@ref pm_arrayRebind::setRebound)<br>
    !>  [setRebilled](@ref pm_arrayRebill::setRebilled)<br>
    !>  [getCoreHalo](@ref pm_arrayInit::getCoreHalo)<br>
    !>  [setCoreHalo](@ref pm_arrayInit::setCoreHalo)<br>
    !>  [setRefilled](@ref pm_arrayRefill::setRefilled)<br>
    !>  [getCentered](@ref pm_arrayCenter::getCentered)<br>
    !>  [setCentered](@ref pm_arrayCenter::setCentered)<br>
    !>  [getPadded](@ref pm_arrayPad::getPadded)<br>
    !>  [setPadded](@ref pm_arrayPad::setPadded)<br>
    !>
    !>  \example{getCoreHalo}
    !>  \include{lineno} example/pm_arrayInit/getCoreHalo/main.F90
    !>  \compilef{getCoreHalo}
    !>  \output{getCoreHalo}
    !>  \include{lineno} example/pm_arrayInit/getCoreHalo/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayInit](@ref test_pm_arrayInit)
    !>
    !>  \todo
    !>  \plow
    !>  This generic interface can be extended to `array` arguments of higher ranks.<br>
    !>
    !>  \todo
    !>  \pmed
    !>  Unlike [setCoreHalo](@ref pm_arrayInit::setCoreHalo), this functional generic interface lacks the option for scalar input `core`.<br>
    !>  This is because of ambiguity created by scalar `core` for output `array` of ranks two and higher.<br>
    !>  This can be resolved once the Fortran standard allows deferred rank function output.<br>
    !>
    !>  \final{getCoreHalo}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getCoreHalo

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE elemental module function getCoreHaloArr_D0_SK5(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        character(1,SKG)            , intent(in)                        :: halo
        character(*,SKG)            , intent(in)                        :: core
        character(size,SKG)                                             :: array
    end function
#endif

#if SK4_ENABLED
    PURE elemental module function getCoreHaloArr_D0_SK4(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        character(1,SKG)            , intent(in)                        :: halo
        character(*,SKG)            , intent(in)                        :: core
        character(size,SKG)                                             :: array
    end function
#endif

#if SK3_ENABLED
    PURE elemental module function getCoreHaloArr_D0_SK3(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        character(1,SKG)            , intent(in)                        :: halo
        character(*,SKG)            , intent(in)                        :: core
        character(size,SKG)                                             :: array
    end function
#endif

#if SK2_ENABLED
    PURE elemental module function getCoreHaloArr_D0_SK2(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        character(1,SKG)            , intent(in)                        :: halo
        character(*,SKG)            , intent(in)                        :: core
        character(size,SKG)                                             :: array
    end function
#endif

#if SK1_ENABLED
    PURE elemental module function getCoreHaloArr_D0_SK1(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        character(1,SKG)            , intent(in)                        :: halo
        character(*,SKG)            , intent(in)                        :: core
        character(size,SKG)                                             :: array
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getCoreHaloArr_D1_SK5(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        character(*,SKG)            , intent(in)                        :: halo
        character(*,SKG)            , intent(in)    , contiguous        :: core(:)
        character(max(len(halo,IK),len(core,IK)),SKG)                   :: array(size)
    end function
#endif

#if SK4_ENABLED
    PURE module function getCoreHaloArr_D1_SK4(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        character(*,SKG)            , intent(in)                        :: halo
        character(*,SKG)            , intent(in)    , contiguous        :: core(:)
        character(max(len(halo,IK),len(core,IK)),SKG)                   :: array(size)
    end function
#endif

#if SK3_ENABLED
    PURE module function getCoreHaloArr_D1_SK3(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        character(*,SKG)            , intent(in)                        :: halo
        character(*,SKG)            , intent(in)    , contiguous        :: core(:)
        character(max(len(halo,IK),len(core,IK)),SKG)                   :: array(size)
    end function
#endif

#if SK2_ENABLED
    PURE module function getCoreHaloArr_D1_SK2(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        character(*,SKG)            , intent(in)                        :: halo
        character(*,SKG)            , intent(in)    , contiguous        :: core(:)
        character(max(len(halo,IK),len(core,IK)),SKG)                   :: array(size)
    end function
#endif

#if SK1_ENABLED
    PURE module function getCoreHaloArr_D1_SK1(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        character(*,SKG)            , intent(in)                        :: halo
        character(*,SKG)            , intent(in)    , contiguous        :: core(:)
        character(max(len(halo,IK),len(core,IK)),SKG)                   :: array(size)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getCoreHaloArr_D1_IK5(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        integer(IKG)                , intent(in)                        :: halo
        integer(IKG)                , intent(in)    , contiguous        :: core(:)
        integer(IKG)                                                    :: array(size)
    end function
#endif

#if IK4_ENABLED
    PURE module function getCoreHaloArr_D1_IK4(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        integer(IKG)                , intent(in)                        :: halo
        integer(IKG)                , intent(in)    , contiguous        :: core(:)
        integer(IKG)                                                    :: array(size)
    end function
#endif

#if IK3_ENABLED
    PURE module function getCoreHaloArr_D1_IK3(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        integer(IKG)                , intent(in)                        :: halo
        integer(IKG)                , intent(in)    , contiguous        :: core(:)
        integer(IKG)                                                    :: array(size)
    end function
#endif

#if IK2_ENABLED
    PURE module function getCoreHaloArr_D1_IK2(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        integer(IKG)                , intent(in)                        :: halo
        integer(IKG)                , intent(in)    , contiguous        :: core(:)
        integer(IKG)                                                    :: array(size)
    end function
#endif

#if IK1_ENABLED
    PURE module function getCoreHaloArr_D1_IK1(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        integer(IKG)                , intent(in)                        :: halo
        integer(IKG)                , intent(in)    , contiguous        :: core(:)
        integer(IKG)                                                    :: array(size)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getCoreHaloArr_D1_LK5(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        logical(LKG)                , intent(in)                        :: halo
        logical(LKG)                , intent(in)    , contiguous        :: core(:)
        logical(LKG)                                                    :: array(size)
    end function
#endif

#if LK4_ENABLED
    PURE module function getCoreHaloArr_D1_LK4(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        logical(LKG)                , intent(in)                        :: halo
        logical(LKG)                , intent(in)    , contiguous        :: core(:)
        logical(LKG)                                                    :: array(size)
    end function
#endif

#if LK3_ENABLED
    PURE module function getCoreHaloArr_D1_LK3(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        logical(LKG)                , intent(in)                        :: halo
        logical(LKG)                , intent(in)    , contiguous        :: core(:)
        logical(LKG)                                                    :: array(size)
    end function
#endif

#if LK2_ENABLED
    PURE module function getCoreHaloArr_D1_LK2(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        logical(LKG)                , intent(in)                        :: halo
        logical(LKG)                , intent(in)    , contiguous        :: core(:)
        logical(LKG)                                                    :: array(size)
    end function
#endif

#if LK1_ENABLED
    PURE module function getCoreHaloArr_D1_LK1(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        logical(LKG)                , intent(in)                        :: halo
        logical(LKG)                , intent(in)    , contiguous        :: core(:)
        logical(LKG)                                                    :: array(size)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCoreHaloArr_D1_CK5(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        complex(CKG)                , intent(in)                        :: halo
        complex(CKG)                , intent(in)    , contiguous        :: core(:)
        complex(CKG)                                                    :: array(size)
    end function
#endif

#if CK4_ENABLED
    PURE module function getCoreHaloArr_D1_CK4(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        complex(CKG)                , intent(in)                        :: halo
        complex(CKG)                , intent(in)    , contiguous        :: core(:)
        complex(CKG)                                                    :: array(size)
    end function
#endif

#if CK3_ENABLED
    PURE module function getCoreHaloArr_D1_CK3(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        complex(CKG)                , intent(in)                        :: halo
        complex(CKG)                , intent(in)    , contiguous        :: core(:)
        complex(CKG)                                                    :: array(size)
    end function
#endif

#if CK2_ENABLED
    PURE module function getCoreHaloArr_D1_CK2(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        complex(CKG)                , intent(in)                        :: halo
        complex(CKG)                , intent(in)    , contiguous        :: core(:)
        complex(CKG)                                                    :: array(size)
    end function
#endif

#if CK1_ENABLED
    PURE module function getCoreHaloArr_D1_CK1(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        complex(CKG)                , intent(in)                        :: halo
        complex(CKG)                , intent(in)    , contiguous        :: core(:)
        complex(CKG)                                                    :: array(size)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCoreHaloArr_D1_RK5(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        real(RKG)                   , intent(in)                        :: halo
        real(RKG)                   , intent(in)    , contiguous        :: core(:)
        real(RKG)                                                       :: array(size)
    end function
#endif

#if RK4_ENABLED
    PURE module function getCoreHaloArr_D1_RK4(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        real(RKG)                   , intent(in)                        :: halo
        real(RKG)                   , intent(in)    , contiguous        :: core(:)
        real(RKG)                                                       :: array(size)
    end function
#endif

#if RK3_ENABLED
    PURE module function getCoreHaloArr_D1_RK3(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        real(RKG)                   , intent(in)                        :: halo
        real(RKG)                   , intent(in)    , contiguous        :: core(:)
        real(RKG)                                                       :: array(size)
    end function
#endif

#if RK2_ENABLED
    PURE module function getCoreHaloArr_D1_RK2(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        real(RKG)                   , intent(in)                        :: halo
        real(RKG)                   , intent(in)    , contiguous        :: core(:)
        real(RKG)                                                       :: array(size)
    end function
#endif

#if RK1_ENABLED
    PURE module function getCoreHaloArr_D1_RK1(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)                 , intent(in)                        :: coffset
        integer(IK)                 , intent(in)                        :: size
        real(RKG)                   , intent(in)                        :: halo
        real(RKG)                   , intent(in)    , contiguous        :: core(:)
        real(RKG)                                                       :: array(size)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getCoreHaloArr_D2_SK5(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        character(*,SKG)            , intent(in)                        :: halo
        character(*,SKG)            , intent(in)    , contiguous        :: core(:,:)
        character(max(len(halo,IK),len(core,IK)),SKG)                   :: array(Size(1), Size(2))
    end function
#endif

#if SK4_ENABLED
    PURE module function getCoreHaloArr_D2_SK4(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        character(*,SKG)            , intent(in)                        :: halo
        character(*,SKG)            , intent(in)    , contiguous        :: core(:,:)
        character(max(len(halo,IK),len(core,IK)),SKG)                   :: array(Size(1), Size(2))
    end function
#endif

#if SK3_ENABLED
    PURE module function getCoreHaloArr_D2_SK3(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        character(*,SKG)            , intent(in)                        :: halo
        character(*,SKG)            , intent(in)    , contiguous        :: core(:,:)
        character(max(len(halo,IK),len(core,IK)),SKG)                   :: array(Size(1), Size(2))
    end function
#endif

#if SK2_ENABLED
    PURE module function getCoreHaloArr_D2_SK2(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        character(*,SKG)            , intent(in)                        :: halo
        character(*,SKG)            , intent(in)    , contiguous        :: core(:,:)
        character(max(len(halo,IK),len(core,IK)),SKG)                   :: array(Size(1), Size(2))
    end function
#endif

#if SK1_ENABLED
    PURE module function getCoreHaloArr_D2_SK1(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        character(*,SKG)            , intent(in)                        :: halo
        character(*,SKG)            , intent(in)    , contiguous        :: core(:,:)
        character(max(len(halo,IK),len(core,IK)),SKG)                   :: array(Size(1), Size(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getCoreHaloArr_D2_IK5(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        integer(IKG)                , intent(in)                        :: halo
        integer(IKG)                , intent(in)    , contiguous        :: core(:,:)
        integer(IKG)                                                    :: array(Size(1), Size(2))
    end function
#endif

#if IK4_ENABLED
    PURE module function getCoreHaloArr_D2_IK4(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        integer(IKG)                , intent(in)                        :: halo
        integer(IKG)                , intent(in)    , contiguous        :: core(:,:)
        integer(IKG)                                                    :: array(Size(1), Size(2))
    end function
#endif

#if IK3_ENABLED
    PURE module function getCoreHaloArr_D2_IK3(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        integer(IKG)                , intent(in)                        :: halo
        integer(IKG)                , intent(in)    , contiguous        :: core(:,:)
        integer(IKG)                                                    :: array(Size(1), Size(2))
    end function
#endif

#if IK2_ENABLED
    PURE module function getCoreHaloArr_D2_IK2(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        integer(IKG)                , intent(in)                        :: halo
        integer(IKG)                , intent(in)    , contiguous        :: core(:,:)
        integer(IKG)                                                    :: array(Size(1), Size(2))
    end function
#endif

#if IK1_ENABLED
    PURE module function getCoreHaloArr_D2_IK1(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        integer(IKG)                , intent(in)                        :: halo
        integer(IKG)                , intent(in)    , contiguous        :: core(:,:)
        integer(IKG)                                                    :: array(Size(1), Size(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getCoreHaloArr_D2_LK5(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        logical(LKG)                , intent(in)                        :: halo
        logical(LKG)                , intent(in)    , contiguous        :: core(:,:)
        logical(LKG)                                                    :: array(Size(1), Size(2))
    end function
#endif

#if LK4_ENABLED
    PURE module function getCoreHaloArr_D2_LK4(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        logical(LKG)                , intent(in)                        :: halo
        logical(LKG)                , intent(in)    , contiguous        :: core(:,:)
        logical(LKG)                                                    :: array(Size(1), Size(2))
    end function
#endif

#if LK3_ENABLED
    PURE module function getCoreHaloArr_D2_LK3(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        logical(LKG)                , intent(in)                        :: halo
        logical(LKG)                , intent(in)    , contiguous        :: core(:,:)
        logical(LKG)                                                    :: array(Size(1), Size(2))
    end function
#endif

#if LK2_ENABLED
    PURE module function getCoreHaloArr_D2_LK2(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        logical(LKG)                , intent(in)                        :: halo
        logical(LKG)                , intent(in)    , contiguous        :: core(:,:)
        logical(LKG)                                                    :: array(Size(1), Size(2))
    end function
#endif

#if LK1_ENABLED
    PURE module function getCoreHaloArr_D2_LK1(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        logical(LKG)                , intent(in)                        :: halo
        logical(LKG)                , intent(in)    , contiguous        :: core(:,:)
        logical(LKG)                                                    :: array(Size(1), Size(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCoreHaloArr_D2_CK5(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        complex(CKG)                , intent(in)                        :: halo
        complex(CKG)                , intent(in)    , contiguous        :: core(:,:)
        complex(CKG)                                                    :: array(Size(1), Size(2))
    end function
#endif

#if CK4_ENABLED
    PURE module function getCoreHaloArr_D2_CK4(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        complex(CKG)                , intent(in)                        :: halo
        complex(CKG)                , intent(in)    , contiguous        :: core(:,:)
        complex(CKG)                                                    :: array(Size(1), Size(2))
    end function
#endif

#if CK3_ENABLED
    PURE module function getCoreHaloArr_D2_CK3(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        complex(CKG)                , intent(in)                        :: halo
        complex(CKG)                , intent(in)    , contiguous        :: core(:,:)
        complex(CKG)                                                    :: array(Size(1), Size(2))
    end function
#endif

#if CK2_ENABLED
    PURE module function getCoreHaloArr_D2_CK2(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        complex(CKG)                , intent(in)                        :: halo
        complex(CKG)                , intent(in)    , contiguous        :: core(:,:)
        complex(CKG)                                                    :: array(Size(1), Size(2))
    end function
#endif

#if CK1_ENABLED
    PURE module function getCoreHaloArr_D2_CK1(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        complex(CKG)                , intent(in)                        :: halo
        complex(CKG)                , intent(in)    , contiguous        :: core(:,:)
        complex(CKG)                                                    :: array(Size(1), Size(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCoreHaloArr_D2_RK5(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        real(RKG)                   , intent(in)                        :: halo
        real(RKG)                   , intent(in)    , contiguous        :: core(:,:)
        real(RKG)                                                       :: array(Size(1), Size(2))
    end function
#endif

#if RK4_ENABLED
    PURE module function getCoreHaloArr_D2_RK4(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        real(RKG)                   , intent(in)                        :: halo
        real(RKG)                   , intent(in)    , contiguous        :: core(:,:)
        real(RKG)                                                       :: array(Size(1), Size(2))
    end function
#endif

#if RK3_ENABLED
    PURE module function getCoreHaloArr_D2_RK3(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        real(RKG)                   , intent(in)                        :: halo
        real(RKG)                   , intent(in)    , contiguous        :: core(:,:)
        real(RKG)                                                       :: array(Size(1), Size(2))
    end function
#endif

#if RK2_ENABLED
    PURE module function getCoreHaloArr_D2_RK2(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        real(RKG)                   , intent(in)                        :: halo
        real(RKG)                   , intent(in)    , contiguous        :: core(:,:)
        real(RKG)                                                       :: array(Size(1), Size(2))
    end function
#endif

#if RK1_ENABLED
    PURE module function getCoreHaloArr_D2_RK1(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)                 , intent(in)                        :: coffset(2)
        integer(IK)                 , intent(in)                        :: Size(2)
        real(RKG)                   , intent(in)                        :: halo
        real(RKG)                   , intent(in)    , contiguous        :: core(:,:)
        real(RKG)                                                       :: array(Size(1), Size(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getCoreHaloArr_D3_SK5(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_SK5
#endif
        use pm_kind, only: SKG => SK5
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        character(*,SKG)            , intent(in)                        :: halo
        character(*,SKG)            , intent(in)    , contiguous        :: core(:,:,:)
        character(max(len(halo,IK),len(core,IK)),SKG)                   :: array(Size(1), Size(2), Size(3))
    end function
#endif

#if SK4_ENABLED
    PURE module function getCoreHaloArr_D3_SK4(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_SK4
#endif
        use pm_kind, only: SKG => SK4
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        character(*,SKG)            , intent(in)                        :: halo
        character(*,SKG)            , intent(in)    , contiguous        :: core(:,:,:)
        character(max(len(halo,IK),len(core,IK)),SKG)                   :: array(Size(1), Size(2), Size(3))
    end function
#endif

#if SK3_ENABLED
    PURE module function getCoreHaloArr_D3_SK3(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_SK3
#endif
        use pm_kind, only: SKG => SK3
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        character(*,SKG)            , intent(in)                        :: halo
        character(*,SKG)            , intent(in)    , contiguous        :: core(:,:,:)
        character(max(len(halo,IK),len(core,IK)),SKG)                   :: array(Size(1), Size(2), Size(3))
    end function
#endif

#if SK2_ENABLED
    PURE module function getCoreHaloArr_D3_SK2(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_SK2
#endif
        use pm_kind, only: SKG => SK2
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        character(*,SKG)            , intent(in)                        :: halo
        character(*,SKG)            , intent(in)    , contiguous        :: core(:,:,:)
        character(max(len(halo,IK),len(core,IK)),SKG)                   :: array(Size(1), Size(2), Size(3))
    end function
#endif

#if SK1_ENABLED
    PURE module function getCoreHaloArr_D3_SK1(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_SK1
#endif
        use pm_kind, only: SKG => SK1
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        character(*,SKG)            , intent(in)                        :: halo
        character(*,SKG)            , intent(in)    , contiguous        :: core(:,:,:)
        character(max(len(halo,IK),len(core,IK)),SKG)                   :: array(Size(1), Size(2), Size(3))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getCoreHaloArr_D3_IK5(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        integer(IKG)                , intent(in)                        :: halo
        integer(IKG)                , intent(in)    , contiguous        :: core(:,:,:)
        integer(IKG)                                                    :: array(Size(1), Size(2), Size(3))
    end function
#endif

#if IK4_ENABLED
    PURE module function getCoreHaloArr_D3_IK4(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        integer(IKG)                , intent(in)                        :: halo
        integer(IKG)                , intent(in)    , contiguous        :: core(:,:,:)
        integer(IKG)                                                    :: array(Size(1), Size(2), Size(3))
    end function
#endif

#if IK3_ENABLED
    PURE module function getCoreHaloArr_D3_IK3(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        integer(IKG)                , intent(in)                        :: halo
        integer(IKG)                , intent(in)    , contiguous        :: core(:,:,:)
        integer(IKG)                                                    :: array(Size(1), Size(2), Size(3))
    end function
#endif

#if IK2_ENABLED
    PURE module function getCoreHaloArr_D3_IK2(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        integer(IKG)                , intent(in)                        :: halo
        integer(IKG)                , intent(in)    , contiguous        :: core(:,:,:)
        integer(IKG)                                                    :: array(Size(1), Size(2), Size(3))
    end function
#endif

#if IK1_ENABLED
    PURE module function getCoreHaloArr_D3_IK1(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        integer(IKG)                , intent(in)                        :: halo
        integer(IKG)                , intent(in)    , contiguous        :: core(:,:,:)
        integer(IKG)                                                    :: array(Size(1), Size(2), Size(3))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getCoreHaloArr_D3_LK5(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_LK5
#endif
        use pm_kind, only: LKG => LK5
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        logical(LKG)                , intent(in)                        :: halo
        logical(LKG)                , intent(in)    , contiguous        :: core(:,:,:)
        logical(LKG)                                                    :: array(Size(1), Size(2), Size(3))
    end function
#endif

#if LK4_ENABLED
    PURE module function getCoreHaloArr_D3_LK4(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_LK4
#endif
        use pm_kind, only: LKG => LK4
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        logical(LKG)                , intent(in)                        :: halo
        logical(LKG)                , intent(in)    , contiguous        :: core(:,:,:)
        logical(LKG)                                                    :: array(Size(1), Size(2), Size(3))
    end function
#endif

#if LK3_ENABLED
    PURE module function getCoreHaloArr_D3_LK3(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_LK3
#endif
        use pm_kind, only: LKG => LK3
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        logical(LKG)                , intent(in)                        :: halo
        logical(LKG)                , intent(in)    , contiguous        :: core(:,:,:)
        logical(LKG)                                                    :: array(Size(1), Size(2), Size(3))
    end function
#endif

#if LK2_ENABLED
    PURE module function getCoreHaloArr_D3_LK2(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_LK2
#endif
        use pm_kind, only: LKG => LK2
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        logical(LKG)                , intent(in)                        :: halo
        logical(LKG)                , intent(in)    , contiguous        :: core(:,:,:)
        logical(LKG)                                                    :: array(Size(1), Size(2), Size(3))
    end function
#endif

#if LK1_ENABLED
    PURE module function getCoreHaloArr_D3_LK1(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_LK1
#endif
        use pm_kind, only: LKG => LK1
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        logical(LKG)                , intent(in)                        :: halo
        logical(LKG)                , intent(in)    , contiguous        :: core(:,:,:)
        logical(LKG)                                                    :: array(Size(1), Size(2), Size(3))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCoreHaloArr_D3_CK5(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        complex(CKG)                , intent(in)                        :: halo
        complex(CKG)                , intent(in)    , contiguous        :: core(:,:,:)
        complex(CKG)                                                    :: array(Size(1), Size(2), Size(3))
    end function
#endif

#if CK4_ENABLED
    PURE module function getCoreHaloArr_D3_CK4(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        complex(CKG)                , intent(in)                        :: halo
        complex(CKG)                , intent(in)    , contiguous        :: core(:,:,:)
        complex(CKG)                                                    :: array(Size(1), Size(2), Size(3))
    end function
#endif

#if CK3_ENABLED
    PURE module function getCoreHaloArr_D3_CK3(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        complex(CKG)                , intent(in)                        :: halo
        complex(CKG)                , intent(in)    , contiguous        :: core(:,:,:)
        complex(CKG)                                                    :: array(Size(1), Size(2), Size(3))
    end function
#endif

#if CK2_ENABLED
    PURE module function getCoreHaloArr_D3_CK2(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        complex(CKG)                , intent(in)                        :: halo
        complex(CKG)                , intent(in)    , contiguous        :: core(:,:,:)
        complex(CKG)                                                    :: array(Size(1), Size(2), Size(3))
    end function
#endif

#if CK1_ENABLED
    PURE module function getCoreHaloArr_D3_CK1(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        complex(CKG)                , intent(in)                        :: halo
        complex(CKG)                , intent(in)    , contiguous        :: core(:,:,:)
        complex(CKG)                                                    :: array(Size(1), Size(2), Size(3))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCoreHaloArr_D3_RK5(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        real(RKG)                   , intent(in)                        :: halo
        real(RKG)                   , intent(in)    , contiguous        :: core(:,:,:)
        real(RKG)                                                       :: array(Size(1), Size(2), Size(3))
    end function
#endif

#if RK4_ENABLED
    PURE module function getCoreHaloArr_D3_RK4(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        real(RKG)                   , intent(in)                        :: halo
        real(RKG)                   , intent(in)    , contiguous        :: core(:,:,:)
        real(RKG)                                                       :: array(Size(1), Size(2), Size(3))
    end function
#endif

#if RK3_ENABLED
    PURE module function getCoreHaloArr_D3_RK3(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        real(RKG)                   , intent(in)                        :: halo
        real(RKG)                   , intent(in)    , contiguous        :: core(:,:,:)
        real(RKG)                                                       :: array(Size(1), Size(2), Size(3))
    end function
#endif

#if RK2_ENABLED
    PURE module function getCoreHaloArr_D3_RK2(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        real(RKG)                   , intent(in)                        :: halo
        real(RKG)                   , intent(in)    , contiguous        :: core(:,:,:)
        real(RKG)                                                       :: array(Size(1), Size(2), Size(3))
    end function
#endif

#if RK1_ENABLED
    PURE module function getCoreHaloArr_D3_RK1(size, core, halo, coffset) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCoreHaloArr_D3_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)                 , intent(in)                        :: coffset(3)
        integer(IK)                 , intent(in)                        :: Size(3)
        real(RKG)                   , intent(in)                        :: halo
        real(RKG)                   , intent(in)    , contiguous        :: core(:,:,:)
        real(RKG)                                                       :: array(Size(1), Size(2), Size(3))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Initialize the rectangular core and halo of an input array of arbitrary rank and shape of arbitrary intrinsic type and kind
    !>  with the corresponding user-specified values.<br>
    !>
    !>  \details
    !>  The **halo of an array of arbitrary rank** is defined as the complement of rectangular core region of arbitrary size within a super-rectangle array.<br>
    !>  The following figure illustrates example cores and halos of two 1D and 2D arrays.<br>
    !>
    !>  \htmlonly
    !>      <img src="pm_arrayInit@scalar.png" style="width:50%;">
    !>  \endhtmlonly
    !>
    !>  \htmlonly
    !>      <img src="pm_arrayInit@array.png" style="width:50%;">
    !>  \endhtmlonly
    !>
    !>  \param[out]     array   :   The output scalar of,<br>
    !>                              <ul>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary length type parameter or,<br>
    !>                              </ul>
    !>                              the output `contiguous` array of rank `1`, `2`, or `3`, of arbitrary shape, of either,<br>
    !>                              <ul>
    !>                                  <li>    type `character` of kind \SKALL, or<br>
    !>                                  <li>    type `integer` of kind \IKALL, or<br>
    !>                                  <li>    type `logical` of kind \LKALL, or<br>
    !>                                  <li>    type `complex` of kind \CKALL, or<br>
    !>                                  <li>    type `real` of kind \RKALL.<br>
    !>                              </ul>
    !>                              On output, the rectangular core and halo of the array will be initialized to the specified input arguments `core` and `halo` respectively.<br>
    !>  \param[in]      core    :   The input scalar, or `contiguous` array of the same rank as the output `array`, of the same type and kind as the output `array`,
    !>                              containing the value to write to the rectangular core of `array`.<br>
    !>                              <ul>
    !>                                  <li>    If `core` is scalar, then the extents along each dimension of `array` are determined by the input vector argument `csize` of size `rank(array)`.<br>
    !>                                          <ul>
    !>                                              <li>    If `array` is a scalar string and the input argument `csize` is present, then `core` must be of type `character` with length type parameter `1`.<br>
    !>                                              <li>    If `array` is a scalar string and the input argument `csize` is missing, then `core` must be of type `character` with length type parameter `len(core) <= len(array)`.<br>
    !>                                          </ul>
    !>                                  <li>    If `core` is an array, then the extents of it (`shape(core)`) delineate the rectangular `core` of `array` into which `core` will be pasted.<br>
    !>                                          <ul>
    !>                                              <li>    If `array` is of type `character`, then the condition `len(core) <= len(array)` must hold.<br>
    !>                                          </ul>
    !>                              </ul>
    !>  \param[in]      halo    :   The input scalar of the same type and kind as the output `array`, containing the value to write to the halo surrounding the specified rectangular core of the `array`.<br>
    !>                              <ul>
    !>                                  <li>    If `array` is a scalar string, then `halo` must be a single `character` (with length type parameter `1`).<br>
    !>                                  <li>    If `array` is an array of type `character`, then `halo` can have arbitrary length type parameter satisfying the condition `len(halo) <= len(array)`.<br>
    !>                              </ul>
    !>  \param[in]      coffset :   The input non-negative scalar or vector of size `rank(array)` of type `integer` of default kind \IK representing the index offset from the first (top-left corner) element of `array`,
    !>                              such that `array(@coffset + 1)` represents the indices of first element of the rectangular core of `array`.<br>
    !>                              <ul>
    !>                                  <li>    If the output `array` is of rank `<= 1`, then `coffset` must be a scalar.<br>
    !>                                  <li>    If the output `array` is of rank `>= 2`, then `coffset` must be a vector of size `rank(array)`.<br>
    !>                              </ul>
    !>  \param[in]      csize   :   The input non-negative scalar or vector of size `rank(array)` of type `integer` of default kind \IK representing the extent of the rectangular core of `array` along each dimension,
    !>                              such that `array(@coffset + @csize)` represents the indices of last element of the rectangular core of `array`.<br>
    !>                              <ul>
    !>                                  <li>    If the output `array` is of rank `<= 1`, then `csize` must be a scalar.<br>
    !>                                  <li>    If the output `array` is of rank `>= 2`, then `csize` must be a vector of size `rank(array)`.<br>
    !>                              </ul>
    !>                              (**optional**. default = `shape(core)`. It must be present **if and only if** the input argument `core` is scalar.)
    !>
    !>  \interface{setCoreHalo}
    !>  \code{.F90}
    !>
    !>      use pm_arrayInit, only: setCoreHalo
    !>
    !>      call setCoreHalo(array(..), core(@shape(array)), halo, coffset(1:rank(array)))
    !>      call setCoreHalo(array(..), core               , halo, coffset(1:rank(array)), csize(1:rank(array)))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `all(0 <= csize)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0 <= coffset)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(coffset + csize <= shape(array))` must hold for the corresponding input arguments.<br>
    !>  The condition `len(core) <= len(array)` must hold for the corresponding input arguments of type `character` of arbitrary rank.<br>
    !>  The condition `len(halo) <= len(array)` must hold for the corresponding input arguments of type `character` of non-zero rank.<br>
    !>  The condition `len(halo) == 1` must hold for the corresponding input arguments of type `character` of zero rank.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  Significant deliberation has gone into this generic interface to ensure its ease of use and future extensibility.<br>
    !>  In particular, a generic interface with arguments `coffset, csize` was preferred over input arguments `lbcore, ubcore`
    !>  representing the lower and upper bounds of the scalar core in the output `array`.<br>
    !>
    !>  \see
    !>  [getCoreHalo](@ref pm_arrayInit::getCoreHalo)<br>
    !>  [setCoreHalo](@ref pm_arrayInit::setCoreHalo)<br>
    !>  [setResized](@ref pm_arrayResize::setResized)<br>
    !>  [setRebound](@ref pm_arrayRebind::setRebound)<br>
    !>  [setRebilled](@ref pm_arrayRebill::setRebilled)<br>
    !>  [setRefilled](@ref pm_arrayRefill::setRefilled)<br>
    !>  [getCentered](@ref pm_arrayCenter::getCentered)<br>
    !>  [setCentered](@ref pm_arrayCenter::setCentered)<br>
    !>  [getPadded](@ref pm_arrayPad::getPadded)<br>
    !>  [setPadded](@ref pm_arrayPad::setPadded)<br>
    !>
    !>  \example{setCoreHalo}
    !>  \include{lineno} example/pm_arrayInit/setCoreHalo/main.F90
    !>  \compilef{setCoreHalo}
    !>  \output{setCoreHalo}
    !>  \include{lineno} example/pm_arrayInit/setCoreHalo/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayInit](@ref test_pm_arrayInit)
    !>
    !>  \todo
    !>  \plow
    !>  This generic interface can be extended to `array` arguments of higher ranks.<br>
    !>
    !>  \todo
    !>  \pvhigh
    !>  The examples of this generic interface incorrectly call the functional interface.<br>
    !>  This must be fixed before release of the library.<br>
    !>
    !>  \todo
    !>  \pmed
    !>  The current interface requires the input `core` to be `contiguous` when it is an array.<br>
    !>  This restriction can be lifted by either removing the `contiguous` attribute of the argument or
    !>  adding additional `lbcore, ubcore` subsetting arguments to the interface.<br>
    !>  The former offers a simple yet powerful solution.<br>
    !>  However, a benchmark must be done to investigate the performance penalty of switching the `contiguous` attribute off.<br>
    !>
    !>  \final{setCoreHalo}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setCoreHalo

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE elemental module subroutine setCoreHaloArr_D0_SK5(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(out)                       :: array
        character(*,SKG)            , intent(in)                        :: core
        character(1,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if SK4_ENABLED
    PURE elemental module subroutine setCoreHaloArr_D0_SK4(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(out)                       :: array
        character(*,SKG)            , intent(in)                        :: core
        character(1,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if SK3_ENABLED
    PURE elemental module subroutine setCoreHaloArr_D0_SK3(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(out)                       :: array
        character(*,SKG)            , intent(in)                        :: core
        character(1,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if SK2_ENABLED
    PURE elemental module subroutine setCoreHaloArr_D0_SK2(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(out)                       :: array
        character(*,SKG)            , intent(in)                        :: core
        character(1,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if SK1_ENABLED
    PURE elemental module subroutine setCoreHaloArr_D0_SK1(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(out)                       :: array
        character(*,SKG)            , intent(in)                        :: core
        character(1,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setCoreHaloArr_D1_SK5(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(out)   , contiguous        :: array(:)
        character(*,SKG)            , intent(in)    , contiguous        :: core(:)
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setCoreHaloArr_D1_SK4(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(out)   , contiguous        :: array(:)
        character(*,SKG)            , intent(in)    , contiguous        :: core(:)
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setCoreHaloArr_D1_SK3(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(out)   , contiguous        :: array(:)
        character(*,SKG)            , intent(in)    , contiguous        :: core(:)
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setCoreHaloArr_D1_SK2(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(out)   , contiguous        :: array(:)
        character(*,SKG)            , intent(in)    , contiguous        :: core(:)
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setCoreHaloArr_D1_SK1(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(out)   , contiguous        :: array(:)
        character(*,SKG)            , intent(in)    , contiguous        :: core(:)
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setCoreHaloArr_D1_IK5(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(out)   , contiguous        :: array(:)
        integer(IKG)                , intent(in)    , contiguous        :: core(:)
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setCoreHaloArr_D1_IK4(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(out)   , contiguous        :: array(:)
        integer(IKG)                , intent(in)    , contiguous        :: core(:)
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setCoreHaloArr_D1_IK3(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(out)   , contiguous        :: array(:)
        integer(IKG)                , intent(in)    , contiguous        :: core(:)
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setCoreHaloArr_D1_IK2(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(out)   , contiguous        :: array(:)
        integer(IKG)                , intent(in)    , contiguous        :: core(:)
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setCoreHaloArr_D1_IK1(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(out)   , contiguous        :: array(:)
        integer(IKG)                , intent(in)    , contiguous        :: core(:)
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setCoreHaloArr_D1_LK5(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(out)   , contiguous        :: array(:)
        logical(LKG)                , intent(in)    , contiguous        :: core(:)
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setCoreHaloArr_D1_LK4(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(out)   , contiguous        :: array(:)
        logical(LKG)                , intent(in)    , contiguous        :: core(:)
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setCoreHaloArr_D1_LK3(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(out)   , contiguous        :: array(:)
        logical(LKG)                , intent(in)    , contiguous        :: core(:)
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setCoreHaloArr_D1_LK2(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(out)   , contiguous        :: array(:)
        logical(LKG)                , intent(in)    , contiguous        :: core(:)
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setCoreHaloArr_D1_LK1(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(out)   , contiguous        :: array(:)
        logical(LKG)                , intent(in)    , contiguous        :: core(:)
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCoreHaloArr_D1_CK5(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(out)   , contiguous        :: array(:)
        complex(CKG)                , intent(in)    , contiguous        :: core(:)
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCoreHaloArr_D1_CK4(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(out)   , contiguous        :: array(:)
        complex(CKG)                , intent(in)    , contiguous        :: core(:)
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCoreHaloArr_D1_CK3(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(out)   , contiguous        :: array(:)
        complex(CKG)                , intent(in)    , contiguous        :: core(:)
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCoreHaloArr_D1_CK2(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(out)   , contiguous        :: array(:)
        complex(CKG)                , intent(in)    , contiguous        :: core(:)
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCoreHaloArr_D1_CK1(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(out)   , contiguous        :: array(:)
        complex(CKG)                , intent(in)    , contiguous        :: core(:)
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCoreHaloArr_D1_RK5(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(out)   , contiguous        :: array(:)
        real(RKG)                   , intent(in)    , contiguous        :: core(:)
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCoreHaloArr_D1_RK4(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(out)   , contiguous        :: array(:)
        real(RKG)                   , intent(in)    , contiguous        :: core(:)
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCoreHaloArr_D1_RK3(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(out)   , contiguous        :: array(:)
        real(RKG)                   , intent(in)    , contiguous        :: core(:)
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCoreHaloArr_D1_RK2(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(out)   , contiguous        :: array(:)
        real(RKG)                   , intent(in)    , contiguous        :: core(:)
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCoreHaloArr_D1_RK1(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(out)   , contiguous        :: array(:)
        real(RKG)                   , intent(in)    , contiguous        :: core(:)
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setCoreHaloArr_D2_SK5(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(out)   , contiguous        :: array(:,:)
        character(*,SKG)            , intent(in)    , contiguous        :: core(:,:)
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setCoreHaloArr_D2_SK4(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(out)   , contiguous        :: array(:,:)
        character(*,SKG)            , intent(in)    , contiguous        :: core(:,:)
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setCoreHaloArr_D2_SK3(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(out)   , contiguous        :: array(:,:)
        character(*,SKG)            , intent(in)    , contiguous        :: core(:,:)
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setCoreHaloArr_D2_SK2(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(out)   , contiguous        :: array(:,:)
        character(*,SKG)            , intent(in)    , contiguous        :: core(:,:)
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setCoreHaloArr_D2_SK1(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(out)   , contiguous        :: array(:,:)
        character(*,SKG)            , intent(in)    , contiguous        :: core(:,:)
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setCoreHaloArr_D2_IK5(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(out)   , contiguous        :: array(:,:)
        integer(IKG)                , intent(in)    , contiguous        :: core(:,:)
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setCoreHaloArr_D2_IK4(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(out)   , contiguous        :: array(:,:)
        integer(IKG)                , intent(in)    , contiguous        :: core(:,:)
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setCoreHaloArr_D2_IK3(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(out)   , contiguous        :: array(:,:)
        integer(IKG)                , intent(in)    , contiguous        :: core(:,:)
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setCoreHaloArr_D2_IK2(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(out)   , contiguous        :: array(:,:)
        integer(IKG)                , intent(in)    , contiguous        :: core(:,:)
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setCoreHaloArr_D2_IK1(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(out)   , contiguous        :: array(:,:)
        integer(IKG)                , intent(in)    , contiguous        :: core(:,:)
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setCoreHaloArr_D2_LK5(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(out)   , contiguous        :: array(:,:)
        logical(LKG)                , intent(in)    , contiguous        :: core(:,:)
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setCoreHaloArr_D2_LK4(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(out)   , contiguous        :: array(:,:)
        logical(LKG)                , intent(in)    , contiguous        :: core(:,:)
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setCoreHaloArr_D2_LK3(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(out)   , contiguous        :: array(:,:)
        logical(LKG)                , intent(in)    , contiguous        :: core(:,:)
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setCoreHaloArr_D2_LK2(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(out)   , contiguous        :: array(:,:)
        logical(LKG)                , intent(in)    , contiguous        :: core(:,:)
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setCoreHaloArr_D2_LK1(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(out)   , contiguous        :: array(:,:)
        logical(LKG)                , intent(in)    , contiguous        :: core(:,:)
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCoreHaloArr_D2_CK5(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(out)   , contiguous        :: array(:,:)
        complex(CKG)                , intent(in)    , contiguous        :: core(:,:)
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCoreHaloArr_D2_CK4(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(out)   , contiguous        :: array(:,:)
        complex(CKG)                , intent(in)    , contiguous        :: core(:,:)
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCoreHaloArr_D2_CK3(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(out)   , contiguous        :: array(:,:)
        complex(CKG)                , intent(in)    , contiguous        :: core(:,:)
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCoreHaloArr_D2_CK2(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(out)   , contiguous        :: array(:,:)
        complex(CKG)                , intent(in)    , contiguous        :: core(:,:)
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCoreHaloArr_D2_CK1(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(out)   , contiguous        :: array(:,:)
        complex(CKG)                , intent(in)    , contiguous        :: core(:,:)
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCoreHaloArr_D2_RK5(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(out)   , contiguous        :: array(:,:)
        real(RKG)                   , intent(in)    , contiguous        :: core(:,:)
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCoreHaloArr_D2_RK4(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(out)   , contiguous        :: array(:,:)
        real(RKG)                   , intent(in)    , contiguous        :: core(:,:)
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCoreHaloArr_D2_RK3(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(out)   , contiguous        :: array(:,:)
        real(RKG)                   , intent(in)    , contiguous        :: core(:,:)
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCoreHaloArr_D2_RK2(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(out)   , contiguous        :: array(:,:)
        real(RKG)                   , intent(in)    , contiguous        :: core(:,:)
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCoreHaloArr_D2_RK1(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(out)   , contiguous        :: array(:,:)
        real(RKG)                   , intent(in)    , contiguous        :: core(:,:)
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setCoreHaloArr_D3_SK5(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(out)   , contiguous        :: array(:,:,:)
        character(*,SKG)            , intent(in)    , contiguous        :: core(:,:,:)
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setCoreHaloArr_D3_SK4(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(out)   , contiguous        :: array(:,:,:)
        character(*,SKG)            , intent(in)    , contiguous        :: core(:,:,:)
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setCoreHaloArr_D3_SK3(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(out)   , contiguous        :: array(:,:,:)
        character(*,SKG)            , intent(in)    , contiguous        :: core(:,:,:)
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setCoreHaloArr_D3_SK2(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(out)   , contiguous        :: array(:,:,:)
        character(*,SKG)            , intent(in)    , contiguous        :: core(:,:,:)
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setCoreHaloArr_D3_SK1(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(out)   , contiguous        :: array(:,:,:)
        character(*,SKG)            , intent(in)    , contiguous        :: core(:,:,:)
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setCoreHaloArr_D3_IK5(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(out)   , contiguous        :: array(:,:,:)
        integer(IKG)                , intent(in)    , contiguous        :: core(:,:,:)
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setCoreHaloArr_D3_IK4(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(out)   , contiguous        :: array(:,:,:)
        integer(IKG)                , intent(in)    , contiguous        :: core(:,:,:)
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setCoreHaloArr_D3_IK3(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(out)   , contiguous        :: array(:,:,:)
        integer(IKG)                , intent(in)    , contiguous        :: core(:,:,:)
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setCoreHaloArr_D3_IK2(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(out)   , contiguous        :: array(:,:,:)
        integer(IKG)                , intent(in)    , contiguous        :: core(:,:,:)
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setCoreHaloArr_D3_IK1(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(out)   , contiguous        :: array(:,:,:)
        integer(IKG)                , intent(in)    , contiguous        :: core(:,:,:)
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setCoreHaloArr_D3_LK5(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(out)   , contiguous        :: array(:,:,:)
        logical(LKG)                , intent(in)    , contiguous        :: core(:,:,:)
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setCoreHaloArr_D3_LK4(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(out)   , contiguous        :: array(:,:,:)
        logical(LKG)                , intent(in)    , contiguous        :: core(:,:,:)
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setCoreHaloArr_D3_LK3(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(out)   , contiguous        :: array(:,:,:)
        logical(LKG)                , intent(in)    , contiguous        :: core(:,:,:)
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setCoreHaloArr_D3_LK2(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(out)   , contiguous        :: array(:,:,:)
        logical(LKG)                , intent(in)    , contiguous        :: core(:,:,:)
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setCoreHaloArr_D3_LK1(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(out)   , contiguous        :: array(:,:,:)
        logical(LKG)                , intent(in)    , contiguous        :: core(:,:,:)
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCoreHaloArr_D3_CK5(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(out)   , contiguous        :: array(:,:,:)
        complex(CKG)                , intent(in)    , contiguous        :: core(:,:,:)
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCoreHaloArr_D3_CK4(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(out)   , contiguous        :: array(:,:,:)
        complex(CKG)                , intent(in)    , contiguous        :: core(:,:,:)
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCoreHaloArr_D3_CK3(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(out)   , contiguous        :: array(:,:,:)
        complex(CKG)                , intent(in)    , contiguous        :: core(:,:,:)
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCoreHaloArr_D3_CK2(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(out)   , contiguous        :: array(:,:,:)
        complex(CKG)                , intent(in)    , contiguous        :: core(:,:,:)
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCoreHaloArr_D3_CK1(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(out)   , contiguous        :: array(:,:,:)
        complex(CKG)                , intent(in)    , contiguous        :: core(:,:,:)
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCoreHaloArr_D3_RK5(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(out)   , contiguous        :: array(:,:,:)
        real(RKG)                   , intent(in)    , contiguous        :: core(:,:,:)
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCoreHaloArr_D3_RK4(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(out)   , contiguous        :: array(:,:,:)
        real(RKG)                   , intent(in)    , contiguous        :: core(:,:,:)
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCoreHaloArr_D3_RK3(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(out)   , contiguous        :: array(:,:,:)
        real(RKG)                   , intent(in)    , contiguous        :: core(:,:,:)
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCoreHaloArr_D3_RK2(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(out)   , contiguous        :: array(:,:,:)
        real(RKG)                   , intent(in)    , contiguous        :: core(:,:,:)
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCoreHaloArr_D3_RK1(array, core, halo, coffset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloArr_D3_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(out)   , contiguous        :: array(:,:,:)
        real(RKG)                   , intent(in)    , contiguous        :: core(:,:,:)
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array))
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
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE elemental module subroutine setCoreHaloSca_D0_SK5(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(out)                       :: array
        character(1,SKG)            , intent(in)                        :: core
        character(1,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if SK4_ENABLED
    PURE elemental module subroutine setCoreHaloSca_D0_SK4(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(out)                       :: array
        character(1,SKG)            , intent(in)                        :: core
        character(1,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if SK3_ENABLED
    PURE elemental module subroutine setCoreHaloSca_D0_SK3(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(out)                       :: array
        character(1,SKG)            , intent(in)                        :: core
        character(1,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if SK2_ENABLED
    PURE elemental module subroutine setCoreHaloSca_D0_SK2(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(out)                       :: array
        character(1,SKG)            , intent(in)                        :: core
        character(1,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if SK1_ENABLED
    PURE elemental module subroutine setCoreHaloSca_D0_SK1(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(out)                       :: array
        character(1,SKG)            , intent(in)                        :: core
        character(1,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setCoreHaloSca_D1_SK5(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(out)   , contiguous        :: array(:)
        character(*,SKG)            , intent(in)                        :: core
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setCoreHaloSca_D1_SK4(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(out)   , contiguous        :: array(:)
        character(*,SKG)            , intent(in)                        :: core
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setCoreHaloSca_D1_SK3(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(out)   , contiguous        :: array(:)
        character(*,SKG)            , intent(in)                        :: core
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setCoreHaloSca_D1_SK2(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(out)   , contiguous        :: array(:)
        character(*,SKG)            , intent(in)                        :: core
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setCoreHaloSca_D1_SK1(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(out)   , contiguous        :: array(:)
        character(*,SKG)            , intent(in)                        :: core
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setCoreHaloSca_D1_IK5(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(out)   , contiguous        :: array(:)
        integer(IKG)                , intent(in)                        :: core
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setCoreHaloSca_D1_IK4(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(out)   , contiguous        :: array(:)
        integer(IKG)                , intent(in)                        :: core
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setCoreHaloSca_D1_IK3(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(out)   , contiguous        :: array(:)
        integer(IKG)                , intent(in)                        :: core
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setCoreHaloSca_D1_IK2(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(out)   , contiguous        :: array(:)
        integer(IKG)                , intent(in)                        :: core
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setCoreHaloSca_D1_IK1(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(out)   , contiguous        :: array(:)
        integer(IKG)                , intent(in)                        :: core
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setCoreHaloSca_D1_LK5(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(out)   , contiguous        :: array(:)
        logical(LKG)                , intent(in)                        :: core
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setCoreHaloSca_D1_LK4(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(out)   , contiguous        :: array(:)
        logical(LKG)                , intent(in)                        :: core
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setCoreHaloSca_D1_LK3(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(out)   , contiguous        :: array(:)
        logical(LKG)                , intent(in)                        :: core
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setCoreHaloSca_D1_LK2(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(out)   , contiguous        :: array(:)
        logical(LKG)                , intent(in)                        :: core
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setCoreHaloSca_D1_LK1(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(out)   , contiguous        :: array(:)
        logical(LKG)                , intent(in)                        :: core
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCoreHaloSca_D1_CK5(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(out)   , contiguous        :: array(:)
        complex(CKG)                , intent(in)                        :: core
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCoreHaloSca_D1_CK4(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(out)   , contiguous        :: array(:)
        complex(CKG)                , intent(in)                        :: core
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCoreHaloSca_D1_CK3(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(out)   , contiguous        :: array(:)
        complex(CKG)                , intent(in)                        :: core
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCoreHaloSca_D1_CK2(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(out)   , contiguous        :: array(:)
        complex(CKG)                , intent(in)                        :: core
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCoreHaloSca_D1_CK1(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(out)   , contiguous        :: array(:)
        complex(CKG)                , intent(in)                        :: core
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCoreHaloSca_D1_RK5(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(out)   , contiguous        :: array(:)
        real(RKG)                   , intent(in)                        :: core
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCoreHaloSca_D1_RK4(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(out)   , contiguous        :: array(:)
        real(RKG)                   , intent(in)                        :: core
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCoreHaloSca_D1_RK3(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(out)   , contiguous        :: array(:)
        real(RKG)                   , intent(in)                        :: core
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCoreHaloSca_D1_RK2(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(out)   , contiguous        :: array(:)
        real(RKG)                   , intent(in)                        :: core
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCoreHaloSca_D1_RK1(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(out)   , contiguous        :: array(:)
        real(RKG)                   , intent(in)                        :: core
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset, csize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setCoreHaloSca_D2_SK5(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(out)   , contiguous        :: array(:,:)
        character(*,SKG)            , intent(in)                        :: core
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setCoreHaloSca_D2_SK4(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(out)   , contiguous        :: array(:,:)
        character(*,SKG)            , intent(in)                        :: core
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setCoreHaloSca_D2_SK3(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(out)   , contiguous        :: array(:,:)
        character(*,SKG)            , intent(in)                        :: core
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setCoreHaloSca_D2_SK2(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(out)   , contiguous        :: array(:,:)
        character(*,SKG)            , intent(in)                        :: core
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setCoreHaloSca_D2_SK1(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(out)   , contiguous        :: array(:,:)
        character(*,SKG)            , intent(in)                        :: core
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setCoreHaloSca_D2_IK5(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(out)   , contiguous        :: array(:,:)
        integer(IKG)                , intent(in)                        :: core
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setCoreHaloSca_D2_IK4(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(out)   , contiguous        :: array(:,:)
        integer(IKG)                , intent(in)                        :: core
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setCoreHaloSca_D2_IK3(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(out)   , contiguous        :: array(:,:)
        integer(IKG)                , intent(in)                        :: core
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setCoreHaloSca_D2_IK2(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(out)   , contiguous        :: array(:,:)
        integer(IKG)                , intent(in)                        :: core
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setCoreHaloSca_D2_IK1(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(out)   , contiguous        :: array(:,:)
        integer(IKG)                , intent(in)                        :: core
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setCoreHaloSca_D2_LK5(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(out)   , contiguous        :: array(:,:)
        logical(LKG)                , intent(in)                        :: core
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setCoreHaloSca_D2_LK4(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(out)   , contiguous        :: array(:,:)
        logical(LKG)                , intent(in)                        :: core
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setCoreHaloSca_D2_LK3(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(out)   , contiguous        :: array(:,:)
        logical(LKG)                , intent(in)                        :: core
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setCoreHaloSca_D2_LK2(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(out)   , contiguous        :: array(:,:)
        logical(LKG)                , intent(in)                        :: core
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setCoreHaloSca_D2_LK1(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(out)   , contiguous        :: array(:,:)
        logical(LKG)                , intent(in)                        :: core
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCoreHaloSca_D2_CK5(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(out)   , contiguous        :: array(:,:)
        complex(CKG)                , intent(in)                        :: core
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCoreHaloSca_D2_CK4(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(out)   , contiguous        :: array(:,:)
        complex(CKG)                , intent(in)                        :: core
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCoreHaloSca_D2_CK3(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(out)   , contiguous        :: array(:,:)
        complex(CKG)                , intent(in)                        :: core
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCoreHaloSca_D2_CK2(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(out)   , contiguous        :: array(:,:)
        complex(CKG)                , intent(in)                        :: core
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCoreHaloSca_D2_CK1(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(out)   , contiguous        :: array(:,:)
        complex(CKG)                , intent(in)                        :: core
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCoreHaloSca_D2_RK5(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(out)   , contiguous        :: array(:,:)
        real(RKG)                   , intent(in)                        :: core
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCoreHaloSca_D2_RK4(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(out)   , contiguous        :: array(:,:)
        real(RKG)                   , intent(in)                        :: core
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCoreHaloSca_D2_RK3(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(out)   , contiguous        :: array(:,:)
        real(RKG)                   , intent(in)                        :: core
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCoreHaloSca_D2_RK2(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(out)   , contiguous        :: array(:,:)
        real(RKG)                   , intent(in)                        :: core
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCoreHaloSca_D2_RK1(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(out)   , contiguous        :: array(:,:)
        real(RKG)                   , intent(in)                        :: core
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setCoreHaloSca_D3_SK5(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(out)   , contiguous        :: array(:,:,:)
        character(*,SKG)            , intent(in)                        :: core
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setCoreHaloSca_D3_SK4(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(out)   , contiguous        :: array(:,:,:)
        character(*,SKG)            , intent(in)                        :: core
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setCoreHaloSca_D3_SK3(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(out)   , contiguous        :: array(:,:,:)
        character(*,SKG)            , intent(in)                        :: core
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setCoreHaloSca_D3_SK2(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(out)   , contiguous        :: array(:,:,:)
        character(*,SKG)            , intent(in)                        :: core
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setCoreHaloSca_D3_SK1(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(out)   , contiguous        :: array(:,:,:)
        character(*,SKG)            , intent(in)                        :: core
        character(*,SKG)            , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setCoreHaloSca_D3_IK5(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(out)   , contiguous        :: array(:,:,:)
        integer(IKG)                , intent(in)                        :: core
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setCoreHaloSca_D3_IK4(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(out)   , contiguous        :: array(:,:,:)
        integer(IKG)                , intent(in)                        :: core
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setCoreHaloSca_D3_IK3(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(out)   , contiguous        :: array(:,:,:)
        integer(IKG)                , intent(in)                        :: core
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setCoreHaloSca_D3_IK2(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(out)   , contiguous        :: array(:,:,:)
        integer(IKG)                , intent(in)                        :: core
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setCoreHaloSca_D3_IK1(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(out)   , contiguous        :: array(:,:,:)
        integer(IKG)                , intent(in)                        :: core
        integer(IKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setCoreHaloSca_D3_LK5(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(out)   , contiguous        :: array(:,:,:)
        logical(LKG)                , intent(in)                        :: core
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setCoreHaloSca_D3_LK4(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(out)   , contiguous        :: array(:,:,:)
        logical(LKG)                , intent(in)                        :: core
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setCoreHaloSca_D3_LK3(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(out)   , contiguous        :: array(:,:,:)
        logical(LKG)                , intent(in)                        :: core
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setCoreHaloSca_D3_LK2(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(out)   , contiguous        :: array(:,:,:)
        logical(LKG)                , intent(in)                        :: core
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setCoreHaloSca_D3_LK1(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(out)   , contiguous        :: array(:,:,:)
        logical(LKG)                , intent(in)                        :: core
        logical(LKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCoreHaloSca_D3_CK5(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(out)   , contiguous        :: array(:,:,:)
        complex(CKG)                , intent(in)                        :: core
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCoreHaloSca_D3_CK4(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(out)   , contiguous        :: array(:,:,:)
        complex(CKG)                , intent(in)                        :: core
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCoreHaloSca_D3_CK3(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(out)   , contiguous        :: array(:,:,:)
        complex(CKG)                , intent(in)                        :: core
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCoreHaloSca_D3_CK2(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(out)   , contiguous        :: array(:,:,:)
        complex(CKG)                , intent(in)                        :: core
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCoreHaloSca_D3_CK1(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(out)   , contiguous        :: array(:,:,:)
        complex(CKG)                , intent(in)                        :: core
        complex(CKG)                , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCoreHaloSca_D3_RK5(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(out)   , contiguous        :: array(:,:,:)
        real(RKG)                   , intent(in)                        :: core
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCoreHaloSca_D3_RK4(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(out)   , contiguous        :: array(:,:,:)
        real(RKG)                   , intent(in)                        :: core
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCoreHaloSca_D3_RK3(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(out)   , contiguous        :: array(:,:,:)
        real(RKG)                   , intent(in)                        :: core
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCoreHaloSca_D3_RK2(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(out)   , contiguous        :: array(:,:,:)
        real(RKG)                   , intent(in)                        :: core
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCoreHaloSca_D3_RK1(array, core, halo, coffset, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCoreHaloSca_D3_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(out)   , contiguous        :: array(:,:,:)
        real(RKG)                   , intent(in)                        :: core
        real(RKG)                   , intent(in)                        :: halo
        integer(IK)                 , intent(in)                        :: coffset(rank(array)), csize(rank(array))
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayInit ! LCOV_EXCL_LINE