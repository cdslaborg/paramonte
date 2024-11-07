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
!>  This module contains procedures and generic interfaces for resizing `allocatable` arrays of various types,
!>  relocating their contents **and filling** the newly added elements with specific values.
!>
!>  \details
!>  In brief,
!>  <ol>
!>      <li>    Use [setResized](@ref pm_arrayResize::setResized) to **resize** strings/arrays **without** changing the lower bound.<br>
!>              The following figure illustrates example resizing of a 1D array and transferal of its contents.<br>
!>              <br>
!>              \htmlonly
!>                  <img src="pm_arrayResize.png" style="width:40%;">
!>              \endhtmlonly
!>              <br>
!>      <li>    Use [setRefilled](@ref pm_arrayRefill::setRefilled) to **resize** strings/arrays **without** changing the lower bound **and** to **initialize** the new elements.<br>
!>              The following figure illustrates example refilling of a 1D array and transferal of its contents.<br>
!>              <br>
!>              \htmlonly
!>                  <img src="pm_arrayRefill.png" style="width:40%;">
!>              \endhtmlonly
!>              <br>
!>      <li>    Use [setRebound](@ref pm_arrayRebind::setRebound) to **resize arrays** by changing their lower and/or upper bounds.<br>
!>              The following figure illustrates example rebinding of a 1D array and transferal of its contents.<br>
!>              <br>
!>              \htmlonly
!>                  <img src="pm_arrayRebind.png" style="width:33%;">
!>              \endhtmlonly
!>              <br>
!>      <li>    Use [setRebilled](@ref pm_arrayRebill::setRebilled) to **resize arrays** by changing their lower and/or upper bounds **and** to **initialize** the new elements.<br>
!>              The following figure illustrates example rebinding and refilling of a 1D array and transferal of its contents.<br>
!>              <br>
!>              \htmlonly
!>                  <img src="pm_arrayRebill.png" style="width:33%;">
!>              \endhtmlonly
!>              <br>
!>  </ol>
!>
!>  \benchmarks
!>
!>  \benchmark{setRefilled_vs_direct, The runtime performance of [setRefilled](@ref pm_arrayRefill::setRefilled) vs. direct reversing}
!>  \include{lineno} benchmark/pm_arrayRefill/setRefilled_vs_direct/main.F90
!>  \compilefb{setRefilled_vs_direct}
!>  \postprocb{setRefilled_vs_direct}
!>  \include{lineno} benchmark/pm_arrayRefill/setRefilled_vs_direct/main.py
!>  \visb{setRefilled_vs_direct}
!>  \image html benchmark/pm_arrayRefill/setRefilled_vs_direct/benchmark.setRefilled_vs_direct.runtime.png width=1000
!>  \image html benchmark/pm_arrayRefill/setRefilled_vs_direct/benchmark.setRefilled_vs_direct.runtime.ratio.png width=1000
!>  \moralb{setRefilled_vs_direct}
!>      -#  The procedures under the generic interface [setRefilled](@ref pm_arrayRefill::setRefilled) are about \f$\ms{10-20%}\f$ faster than
!>          naively copying and resizing via temporary arrays and initializing the values, simply because
!>          [setRefilled](@ref pm_arrayRefill::setRefilled) avoids an extra unnecessary reallocation and data copy.<br>
!>
!>  \test
!>  [test_pm_arrayRefill](@ref test_pm_arrayRefill)
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayRefill

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_arrayRefill"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Allocate or resize (shrink or expand) and refill an input `allocatable` scalar string or
    !>  array of rank `1..3` to an arbitrary size while preserving the original contents or a subset of it.<br>
    !>
    !>  \details
    !>  The new `array` size is set to twice its current size or, to the requested input `size(..)`.<br>
    !>  The array contents or a requested subset of it are kept in the original indices in the output resized array or
    !>  shifted to a new starting location `lbc` in the output `array`.<br>
    !>  The rest of the elements (including the newly-added elements) are filled with the user-specified `fill`.<br>
    !>
    !>  The following figure illustrates example resizing and refilling of a 1D array and transferal of its contents.<br>
    !>  <br>
    !>  \htmlonly
    !>      <img src="pm_arrayRefill.png" style="width:40%;">
    !>  \endhtmlonly
    !>  <br>
    !>
    !>  \param[inout]   array   :   The input/output `allocatable` scalar of<br>
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL type `character` of kind \SKALL
    !>                              </ol>
    !>                              or array of rank `1..3` of either<br>
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL
    !>                                  <li>    type `integer` of kind \IKALL
    !>                                  <li>    type `logical` of kind \LKALL
    !>                                  <li>    type `complex` of kind \CKALL
    !>                                  <li>    type `real` of kind \RKALL
    !>                              </ol>
    !>                              On output, the array will be (re)allocated to the requested new `size` with the **same lower bound** as before (or `1` if unallocated).
    !>  \param[in]      fill    :   The input scalar of the same type and kind as the input `array` containing the value to fill the new
    !>                              elements (if any) of `array`. If `array` is of type `character`, then <br>
    !>                              <ol>
    !>                                  <li>    The equality `len(fill) == 1` must also hold if `array` is of zero rank (i.e., a scalar string).
    !>                                  <li>    The equality `len(fill) <= len(array)` must also hold if `array` is of non-zero rank.
    !>                              </ol>
    !>  \param[in]      size    :   The input non-negative scalar or array of type `integer` of default kind \IK representing the new size of the output array.<br>
    !>                              <ol>
    !>                                  <li>    If `array` is a scalar or array of rank `1`, then `size` must be a scalar.<br>
    !>                                  <li>    If `array` is an array of rank `> 1`, then `size` must be a vector of the same length as `rank(array)`.<br>
    !>                              </ol>
    !>                              (**optional**, default = `2 * len/shape(array)` where the condition `0 < len/shape(array)` must hold, otherwise infinite loops within the program can occur.)
    !>  \param[in]      lbc     :   The input scalar or array of type `integer` of default kind \IK,
    !>                              representing the <b>L</b>ower <b>B</b>ound(s) of the <b>C</b>ontents in the newly resized output `array`.<br>
    !>                              <ol>
    !>                                  <li>    If `array` is a scalar or array of rank `1`, then `lbc` must be a scalar.<br>
    !>                                  <li>    If `array` is an array of rank `> 1`, then `lbc` must be a vector of the same length as `rank(array)`.<br>
    !>                              </ol>
    !>                              (**optional**, default = `lbcold`. It can be present **only if** the `size` argument is also present.)
    !>  \param[in]      lbcold  :   The input scalar or array of type `integer` of default kind \IK, representing the <b>L</b>ower <b>B</b>ound(s) of the <b>C</b>ontents
    !>                              in the original (<b>old</b>) input `array` that is to be copied to the newly allocated output `array` starting at the new lower bound(s) `lbc`.<br>
    !>                              (**optional**, default = `ubound(array)`. If `array` is a scalar string, then default = `1`.
    !>                              It can be present **only if** the `size`, `lbc`, and `ubcold` input arguments are also present.)
    !>  \param[in]      ubcold  :   The input scalar or array of type `integer` of default kind \IK, representing the <b>U</b>pper <b>B</b>ound(s) of the <b>C</b>ontents
    !>                              in the original (<b>old</b>) input `array` that is to be copied to the newly allocated output `array` starting at the new lower bound(s) `lbc`.<br>
    !>                              (**optional**, default = `ubound(array)`. If `array` is a scalar string, then default = `len(array)`
    !>                              It can be present **only if** the `size` and `lbc` and `lbcold` input arguments are also present.)
    !>  \param[out]     failed  :   The output scalar `logical` of default kind \LK that is `.false.` if and only if the requested array resizing is successful,
    !>                              otherwise it is set to `.true.` to signal the occurrence of an allocation error.<br>
    !>                              The value of `failed` is `.true.` only if the `stat` argument returned by the Fortran intrinsic `allocate()` statement is non-zero.<br>
    !>                              (**optional**, if missing and an allocation error occurs, the processor dictates the program behavior (normally execution stops).)
    !>  \param[out]     errmsg  :   The output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If the `optional` output argument `failed` is present and an error occurs, `errmsg` will be set to a message describing the nature of the error.<br>
    !>                              This behavior conforms with the standard Fortran behavior for the intrinsic `allocate()` statement.<br>
    !>                              A length type parameter of `127` or more for `errmsg` should be sufficient for capturing most if not all error messages in entirety.<br>
    !>                              (**optional**. Its presence is relevant **if and only if** the `optional` output argument `failed` is also present.)
    !>
    !>  \interface{setRefilled}
    !>  \code{.F90}
    !>
    !>      use pm_arrayRefill, only: setRefilled
    !>
    !>      call setRefilled(array, fill, failed = failed, errmsg = errmsg)                             ! Refill to twice as large with the same lower bound as before and optionally gracefully return if reallocation fails.
    !>      call setRefilled(array, fill, size, failed = failed, errmsg = errmsg)                       ! Refill to the requested `size` with the same lower bound as before and optionally gracefully return if reallocation fails.
    !>      call setRefilled(array, fill, size, lbc, failed = failed, errmsg = errmsg)                  ! Refill to the requested `size`, write `array` to the output `array(\@lbc:)`, optionally gracefully return upon failure.
    !>      call setRefilled(array, fill, size, lbc, lbcold, ubcold, failed = failed, errmsg = errmsg)  ! Refill to the requested `size` with the same lower bound, write `array(\@lbcold:ubcold)` to `array(\@lbc:\@lbc-lbcold+ubcold)`, optionally gracefully return upon failure.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `allocated(array)` must hold (the input `array` must be preallocated) when any or all of the optional input arguments `lbc, lbcold, ubcold` are present or when the `size` argument is missing.<br>
    !>  The condition `all(0_IK <= size)` must hold for the corresponding input argument.<br>
    !>  The condition `all(lbound(array) <= lbcold .and. lbcold <= ubound(array))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(lbound(array) <= ubcold .and. ubcold <= ubound(array))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(lbound(array) <= lbc)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(lbc + ubcold - lbcold <= lbound(array) + size - 1)` must hold for the corresponding input arguments (i.e., the upper bound(s) of contents cannot overflow the upper bound(s) of the new array).<br>
    !>  The equality `len(fill) == 1` must also hold if `array` is of zero rank (i.e., a scalar string).<br>
    !>  The equality `len(fill) <= len(array)` must also hold if `array` is of non-zero rank.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  If the input `array` is unallocated, it will be allocated to the requested shape and filled entirely with `fill`, equivalent to `allocate(array(@size), source = fill)`.<br>
    !>
    !>  \note
    !>  If the initialization of the new elements with `fill` is not necessary,
    !>  use [setResized](@ref pm_arrayResize::setResized) to resize arrays without initialization.<br>
    !>
    !>  \note
    !>  While shrinking an array using this generic interface is pointless dues to the lack of the use of `fill` argument for any purpose, shrinking (without actually doing any refills) is allowed.<br>
    !>  This behavior is quite useful for graceful handling of situations where the task (shrinking or expanding) is only determined at runtime or can change from one call to another.<br>
    !>  Otherwise, use [setResized](@ref pm_arrayResize::setResized) to shrink arrays or resize arrays without initialization.<br>
    !>
    !>  \note
    !>  The sole purpose of this generic interface is to provide a **convenient** but **fast** method of resizing allocatable
    !>  arrays without losing the contents of the array and refilling the new elements with requested `fill`.<br>
    !>
    !>  \devnote
    !>  An optional dummy argument `stat` (instead of `failed`) for the procedures of this generic interface are impossible as it creates ambiguous interfaces.<br>
    !>
    !>  \see
    !>  [setResized](@ref pm_arrayResize::setResized)<br>
    !>  [setRebound](@ref pm_arrayRebind::setRebound)<br>
    !>  [setRefilled](@ref pm_arrayRefill::setRefilled)<br>
    !>  [setRebilled](@ref pm_arrayRebill::setRebilled)<br>
    !>  [getCoreHalo](@ref pm_arrayInit::getCoreHalo)<br>
    !>  [setCoreHalo](@ref pm_arrayInit::setCoreHalo)<br>
    !>  [getCentered](@ref pm_arrayCenter::getCentered)<br>
    !>  [setCentered](@ref pm_arrayCenter::setCentered)<br>
    !>  [getPadded](@ref pm_arrayPad::getPadded)<br>
    !>  [setPadded](@ref pm_arrayPad::setPadded)<br>
    !>
    !>  \example{setRefilled}
    !>  \include{lineno} example/pm_arrayRefill/setRefilled/main.F90
    !>  \compilef{setRefilled}
    !>  \output{setRefilled}
    !>  \include{lineno} example/pm_arrayRefill/setRefilled/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayRefill](@ref test_pm_arrayRefill)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{10-12}
    !>  \desc
    !>  There is an annoying gfortran bug concerning allocation of allocatable arrays of strings with assumed length type parameter.<br>
    !>  The typical compiler error message is `around line 230: Error allocating 283223642230368 bytes: Cannot allocate memory`.<br>
    !>  This requires the allocation statement be explicit for `character` arrays of non-zero rank.<br>
    !>  This makes the already complex code superbly more complex and messy.<br>
    !>  \remedy
    !>  For now, the `allocatable` arrays of type `character` are allocated with explicit shape in the allocation statement.<br>
    !>  This explicit allocation for `character` types must be removed and replaced with the generic allocation once the bug is resolved.<br>
    !>
    !>  \todo
    !>  \pvlow
    !>  This generic interface can be extended to arrays of higher ranks than currently supported.<br>
    !>
    !>  \todo
    !>  \pmed
    !>  This generic interface should be extended to arrays of container type as done in [pm_arrayResize](@ref pm_arrayResize).<br>
    !>
    !>  \final{setRefilled}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setRefilled

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRefilledDDDD_D0_SK5(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(1,SKG)            , intent(in)                        :: fill
        character(:,SKG)            , intent(inout) , allocatable       :: array
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRefilledDDDD_D0_SK4(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(1,SKG)            , intent(in)                        :: fill
        character(:,SKG)            , intent(inout) , allocatable       :: array
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRefilledDDDD_D0_SK3(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(1,SKG)            , intent(in)                        :: fill
        character(:,SKG)            , intent(inout) , allocatable       :: array
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRefilledDDDD_D0_SK2(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(1,SKG)            , intent(in)                        :: fill
        character(:,SKG)            , intent(inout) , allocatable       :: array
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRefilledDDDD_D0_SK1(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(1,SKG)            , intent(in)                        :: fill
        character(:,SKG)            , intent(inout) , allocatable       :: array
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRefilledDDDD_D1_SK5(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRefilledDDDD_D1_SK4(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRefilledDDDD_D1_SK3(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRefilledDDDD_D1_SK2(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRefilledDDDD_D1_SK1(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRefilledDDDD_D1_IK5(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRefilledDDDD_D1_IK4(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRefilledDDDD_D1_IK3(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRefilledDDDD_D1_IK2(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRefilledDDDD_D1_IK1(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRefilledDDDD_D1_LK5(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRefilledDDDD_D1_LK4(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRefilledDDDD_D1_LK3(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRefilledDDDD_D1_LK2(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRefilledDDDD_D1_LK1(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRefilledDDDD_D1_CK5(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRefilledDDDD_D1_CK4(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRefilledDDDD_D1_CK3(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRefilledDDDD_D1_CK2(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRefilledDDDD_D1_CK1(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRefilledDDDD_D1_RK5(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRefilledDDDD_D1_RK4(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRefilledDDDD_D1_RK3(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRefilledDDDD_D1_RK2(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRefilledDDDD_D1_RK1(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRefilledDDDD_D2_SK5(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRefilledDDDD_D2_SK4(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRefilledDDDD_D2_SK3(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRefilledDDDD_D2_SK2(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRefilledDDDD_D2_SK1(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRefilledDDDD_D2_IK5(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRefilledDDDD_D2_IK4(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRefilledDDDD_D2_IK3(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRefilledDDDD_D2_IK2(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRefilledDDDD_D2_IK1(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRefilledDDDD_D2_LK5(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRefilledDDDD_D2_LK4(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRefilledDDDD_D2_LK3(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRefilledDDDD_D2_LK2(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRefilledDDDD_D2_LK1(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRefilledDDDD_D2_CK5(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRefilledDDDD_D2_CK4(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRefilledDDDD_D2_CK3(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRefilledDDDD_D2_CK2(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRefilledDDDD_D2_CK1(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRefilledDDDD_D2_RK5(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRefilledDDDD_D2_RK4(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRefilledDDDD_D2_RK3(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRefilledDDDD_D2_RK2(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRefilledDDDD_D2_RK1(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRefilledDDDD_D3_SK5(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRefilledDDDD_D3_SK4(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRefilledDDDD_D3_SK3(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRefilledDDDD_D3_SK2(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRefilledDDDD_D3_SK1(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRefilledDDDD_D3_IK5(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRefilledDDDD_D3_IK4(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRefilledDDDD_D3_IK3(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRefilledDDDD_D3_IK2(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRefilledDDDD_D3_IK1(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRefilledDDDD_D3_LK5(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRefilledDDDD_D3_LK4(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRefilledDDDD_D3_LK3(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRefilledDDDD_D3_LK2(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRefilledDDDD_D3_LK1(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRefilledDDDD_D3_CK5(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRefilledDDDD_D3_CK4(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRefilledDDDD_D3_CK3(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRefilledDDDD_D3_CK2(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRefilledDDDD_D3_CK1(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRefilledDDDD_D3_RK5(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRefilledDDDD_D3_RK4(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRefilledDDDD_D3_RK3(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRefilledDDDD_D3_RK2(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRefilledDDDD_D3_RK1(array, fill, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledDDDD_D3_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
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

#if SK5_ENABLED
    PURE module subroutine setRefilledSDDD_D0_SK5(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(1,SKG)            , intent(in)                        :: fill
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRefilledSDDD_D0_SK4(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(1,SKG)            , intent(in)                        :: fill
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRefilledSDDD_D0_SK3(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(1,SKG)            , intent(in)                        :: fill
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRefilledSDDD_D0_SK2(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(1,SKG)            , intent(in)                        :: fill
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRefilledSDDD_D0_SK1(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(1,SKG)            , intent(in)                        :: fill
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRefilledSDDD_D1_SK5(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRefilledSDDD_D1_SK4(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRefilledSDDD_D1_SK3(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRefilledSDDD_D1_SK2(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRefilledSDDD_D1_SK1(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRefilledSDDD_D1_IK5(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRefilledSDDD_D1_IK4(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRefilledSDDD_D1_IK3(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRefilledSDDD_D1_IK2(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRefilledSDDD_D1_IK1(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRefilledSDDD_D1_LK5(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRefilledSDDD_D1_LK4(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRefilledSDDD_D1_LK3(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRefilledSDDD_D1_LK2(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRefilledSDDD_D1_LK1(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRefilledSDDD_D1_CK5(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRefilledSDDD_D1_CK4(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRefilledSDDD_D1_CK3(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRefilledSDDD_D1_CK2(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRefilledSDDD_D1_CK1(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRefilledSDDD_D1_RK5(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRefilledSDDD_D1_RK4(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRefilledSDDD_D1_RK3(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRefilledSDDD_D1_RK2(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRefilledSDDD_D1_RK1(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRefilledSDDD_D2_SK5(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRefilledSDDD_D2_SK4(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRefilledSDDD_D2_SK3(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRefilledSDDD_D2_SK2(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRefilledSDDD_D2_SK1(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRefilledSDDD_D2_IK5(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRefilledSDDD_D2_IK4(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRefilledSDDD_D2_IK3(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRefilledSDDD_D2_IK2(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRefilledSDDD_D2_IK1(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRefilledSDDD_D2_LK5(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRefilledSDDD_D2_LK4(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRefilledSDDD_D2_LK3(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRefilledSDDD_D2_LK2(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRefilledSDDD_D2_LK1(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRefilledSDDD_D2_CK5(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRefilledSDDD_D2_CK4(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRefilledSDDD_D2_CK3(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRefilledSDDD_D2_CK2(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRefilledSDDD_D2_CK1(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRefilledSDDD_D2_RK5(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRefilledSDDD_D2_RK4(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRefilledSDDD_D2_RK3(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRefilledSDDD_D2_RK2(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRefilledSDDD_D2_RK1(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRefilledSDDD_D3_SK5(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRefilledSDDD_D3_SK4(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRefilledSDDD_D3_SK3(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRefilledSDDD_D3_SK2(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRefilledSDDD_D3_SK1(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRefilledSDDD_D3_IK5(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRefilledSDDD_D3_IK4(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRefilledSDDD_D3_IK3(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRefilledSDDD_D3_IK2(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRefilledSDDD_D3_IK1(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRefilledSDDD_D3_LK5(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRefilledSDDD_D3_LK4(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRefilledSDDD_D3_LK3(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRefilledSDDD_D3_LK2(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRefilledSDDD_D3_LK1(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRefilledSDDD_D3_CK5(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRefilledSDDD_D3_CK4(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRefilledSDDD_D3_CK3(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRefilledSDDD_D3_CK2(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRefilledSDDD_D3_CK1(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRefilledSDDD_D3_RK5(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRefilledSDDD_D3_RK4(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRefilledSDDD_D3_RK3(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRefilledSDDD_D3_RK2(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRefilledSDDD_D3_RK1(array, fill, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSDDD_D3_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
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

#if SK5_ENABLED
    PURE module subroutine setRefilledSLDD_D0_SK5(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(1,SKG)            , intent(in)                        :: fill
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRefilledSLDD_D0_SK4(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(1,SKG)            , intent(in)                        :: fill
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRefilledSLDD_D0_SK3(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(1,SKG)            , intent(in)                        :: fill
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRefilledSLDD_D0_SK2(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(1,SKG)            , intent(in)                        :: fill
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRefilledSLDD_D0_SK1(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(1,SKG)            , intent(in)                        :: fill
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRefilledSLDD_D1_SK5(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRefilledSLDD_D1_SK4(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRefilledSLDD_D1_SK3(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRefilledSLDD_D1_SK2(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRefilledSLDD_D1_SK1(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRefilledSLDD_D1_IK5(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRefilledSLDD_D1_IK4(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRefilledSLDD_D1_IK3(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRefilledSLDD_D1_IK2(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRefilledSLDD_D1_IK1(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRefilledSLDD_D1_LK5(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRefilledSLDD_D1_LK4(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRefilledSLDD_D1_LK3(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRefilledSLDD_D1_LK2(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRefilledSLDD_D1_LK1(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRefilledSLDD_D1_CK5(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRefilledSLDD_D1_CK4(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRefilledSLDD_D1_CK3(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRefilledSLDD_D1_CK2(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRefilledSLDD_D1_CK1(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRefilledSLDD_D1_RK5(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRefilledSLDD_D1_RK4(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRefilledSLDD_D1_RK3(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRefilledSLDD_D1_RK2(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRefilledSLDD_D1_RK1(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRefilledSLDD_D2_SK5(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRefilledSLDD_D2_SK4(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRefilledSLDD_D2_SK3(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRefilledSLDD_D2_SK2(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRefilledSLDD_D2_SK1(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRefilledSLDD_D2_IK5(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRefilledSLDD_D2_IK4(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRefilledSLDD_D2_IK3(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRefilledSLDD_D2_IK2(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRefilledSLDD_D2_IK1(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRefilledSLDD_D2_LK5(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRefilledSLDD_D2_LK4(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRefilledSLDD_D2_LK3(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRefilledSLDD_D2_LK2(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRefilledSLDD_D2_LK1(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRefilledSLDD_D2_CK5(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRefilledSLDD_D2_CK4(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRefilledSLDD_D2_CK3(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRefilledSLDD_D2_CK2(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRefilledSLDD_D2_CK1(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRefilledSLDD_D2_RK5(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRefilledSLDD_D2_RK4(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRefilledSLDD_D2_RK3(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRefilledSLDD_D2_RK2(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRefilledSLDD_D2_RK1(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRefilledSLDD_D3_SK5(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRefilledSLDD_D3_SK4(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRefilledSLDD_D3_SK3(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRefilledSLDD_D3_SK2(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRefilledSLDD_D3_SK1(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRefilledSLDD_D3_IK5(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRefilledSLDD_D3_IK4(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRefilledSLDD_D3_IK3(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRefilledSLDD_D3_IK2(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRefilledSLDD_D3_IK1(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRefilledSLDD_D3_LK5(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRefilledSLDD_D3_LK4(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRefilledSLDD_D3_LK3(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRefilledSLDD_D3_LK2(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRefilledSLDD_D3_LK1(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRefilledSLDD_D3_CK5(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRefilledSLDD_D3_CK4(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRefilledSLDD_D3_CK3(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRefilledSLDD_D3_CK2(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRefilledSLDD_D3_CK1(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRefilledSLDD_D3_RK5(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRefilledSLDD_D3_RK4(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRefilledSLDD_D3_RK3(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRefilledSLDD_D3_RK2(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRefilledSLDD_D3_RK1(array, fill, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLDD_D3_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
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

#if SK5_ENABLED
    PURE module subroutine setRefilledSLLU_D0_SK5(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(1,SKG)            , intent(in)                        :: fill
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRefilledSLLU_D0_SK4(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(1,SKG)            , intent(in)                        :: fill
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRefilledSLLU_D0_SK3(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(1,SKG)            , intent(in)                        :: fill
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRefilledSLLU_D0_SK2(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(1,SKG)            , intent(in)                        :: fill
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRefilledSLLU_D0_SK1(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(1,SKG)            , intent(in)                        :: fill
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRefilledSLLU_D1_SK5(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRefilledSLLU_D1_SK4(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRefilledSLLU_D1_SK3(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRefilledSLLU_D1_SK2(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRefilledSLLU_D1_SK1(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRefilledSLLU_D1_IK5(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRefilledSLLU_D1_IK4(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRefilledSLLU_D1_IK3(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRefilledSLLU_D1_IK2(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRefilledSLLU_D1_IK1(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRefilledSLLU_D1_LK5(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRefilledSLLU_D1_LK4(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRefilledSLLU_D1_LK3(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRefilledSLLU_D1_LK2(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRefilledSLLU_D1_LK1(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRefilledSLLU_D1_CK5(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRefilledSLLU_D1_CK4(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRefilledSLLU_D1_CK3(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRefilledSLLU_D1_CK2(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRefilledSLLU_D1_CK1(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRefilledSLLU_D1_RK5(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRefilledSLLU_D1_RK4(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRefilledSLLU_D1_RK3(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRefilledSLLU_D1_RK2(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRefilledSLLU_D1_RK1(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRefilledSLLU_D2_SK5(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRefilledSLLU_D2_SK4(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRefilledSLLU_D2_SK3(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRefilledSLLU_D2_SK2(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRefilledSLLU_D2_SK1(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRefilledSLLU_D2_IK5(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRefilledSLLU_D2_IK4(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRefilledSLLU_D2_IK3(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRefilledSLLU_D2_IK2(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRefilledSLLU_D2_IK1(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRefilledSLLU_D2_LK5(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRefilledSLLU_D2_LK4(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRefilledSLLU_D2_LK3(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRefilledSLLU_D2_LK2(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRefilledSLLU_D2_LK1(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRefilledSLLU_D2_CK5(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRefilledSLLU_D2_CK4(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRefilledSLLU_D2_CK3(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRefilledSLLU_D2_CK2(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRefilledSLLU_D2_CK1(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRefilledSLLU_D2_RK5(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRefilledSLLU_D2_RK4(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRefilledSLLU_D2_RK3(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRefilledSLLU_D2_RK2(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRefilledSLLU_D2_RK1(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRefilledSLLU_D3_SK5(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRefilledSLLU_D3_SK4(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRefilledSLLU_D3_SK3(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRefilledSLLU_D3_SK2(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRefilledSLLU_D3_SK1(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                        :: fill
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRefilledSLLU_D3_IK5(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRefilledSLLU_D3_IK4(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRefilledSLLU_D3_IK3(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRefilledSLLU_D3_IK2(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRefilledSLLU_D3_IK1(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)                        :: fill
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRefilledSLLU_D3_LK5(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRefilledSLLU_D3_LK4(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRefilledSLLU_D3_LK3(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRefilledSLLU_D3_LK2(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRefilledSLLU_D3_LK1(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)                        :: fill
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRefilledSLLU_D3_CK5(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRefilledSLLU_D3_CK4(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRefilledSLLU_D3_CK3(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRefilledSLLU_D3_CK2(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRefilledSLLU_D3_CK1(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)                        :: fill
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRefilledSLLU_D3_RK5(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRefilledSLLU_D3_RK4(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRefilledSLLU_D3_RK3(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRefilledSLLU_D3_RK2(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRefilledSLLU_D3_RK1(array, fill, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefilledSLLU_D3_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)                        :: fill
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
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

end module pm_arrayRefill ! LCOV_EXCL_LINE