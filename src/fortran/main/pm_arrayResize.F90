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
!>  This module contains procedures and generic interfaces for resizing `allocatable` arrays of various types and relocating their contents,
!>  **without** initializing or filling the newly added elements with specific values.<br>
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
!>  \benchmark{setResized_vs_direct, The runtime performance of [setResized](@ref pm_arrayResize::setResized) vs. direct reversing}
!>  \include{lineno} benchmark/pm_arrayResize/setResized_vs_direct/main.F90
!>  \compilefb{setResized_vs_direct}
!>  \postprocb{setResized_vs_direct}
!>  \include{lineno} benchmark/pm_arrayResize/setResized_vs_direct/main.py
!>  \visb{setResized_vs_direct}
!>  \image html benchmark/pm_arrayResize/setResized_vs_direct/benchmark.setResized_vs_direct.runtime.png width=1000
!>  \image html benchmark/pm_arrayResize/setResized_vs_direct/benchmark.setResized_vs_direct.runtime.ratio.png width=1000
!>  \moralb{setResized_vs_direct}
!>  <ol>
!>      <li>    The procedures under the generic interface [setResized](@ref pm_arrayResize::setResized) are about \f$\ms{10-20%}\f$ faster than
!>              naively copying and resizing via temporary arrays, simply because [setResized](@ref pm_arrayResize::setResized)
!>              avoids an extra unnecessary reallocation and data copy.<br>
!>  </ol>
!>
!>  \test
!>  [test_pm_arrayResize](@ref test_pm_arrayResize)
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayResize

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_arrayResize"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Allocate or resize (shrink or expand) an input `allocatable` scalar string or
    !>  array of rank `1..3` to an arbitrary size while preserving the original contents or a subset of it.<br>
    !>
    !>  \details
    !>  The new `array` size is set to twice its current size or, to the requested input `size(..)`.<br>
    !>  The array contents or a requested subset of it are kept in the original indices in the output resized array or
    !>  shifted to a new starting location `lbc` in the output `array`.<br>
    !>
    !>  The following figure illustrates example resizing of a 1D array and transferal of its contents.<br>
    !>  <br>
    !>  \htmlonly
    !>      <img src="pm_arrayResize.png" style="width:40%;">
    !>  \endhtmlonly
    !>  <br>
    !>
    !>  \param[inout]   array   :   The input/output `allocatable` scalar of<br>
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL type `character` of kind \SKALL
    !>                              </ol>
    !>                              or array of rank `1..3` of either<br>
    !>                              <ol>
    !>                                  <li>    type [css_pdt](@ref pm_container::css_pdt) (parameterized container of string of kind \SKALL) or,<br>
    !>                                  <li>    type [css_type](@ref pm_container::css_type) (container of string of default kind \SK) or,<br>
    !>                                  <li>    type `character` of kind \SKALL
    !>                                  <li>    type `integer` of kind \IKALL
    !>                                  <li>    type `logical` of kind \LKALL
    !>                                  <li>    type `complex` of kind \CKALL
    !>                                  <li>    type `real` of kind \RKALL
    !>                              </ol>
    !>                              On output, the array will be (re)allocated to the requested new `size` with the **same lower bound** as before (or `1` if unallocated).
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
    !>  \interface{setResized}
    !>  \code{.F90}
    !>
    !>      use pm_arrayResize, only: setResized
    !>
    !>      call setResized(array, failed = failed, errmsg = errmsg)                            ! Resize to twice as large with the same lower bound as before and optionally gracefully return if reallocation fails.
    !>      call setResized(array, size, failed = failed, errmsg = errmsg)                      ! Resize to the requested `size` with the same lower bound as before and optionally gracefully return if reallocation fails.
    !>      call setResized(array, size, lbc, failed = failed, errmsg = errmsg)                 ! Resize to the requested `size`, write `array` to the output `array(\@lbc:)`, optionally gracefully return upon failure.
    !>      call setResized(array, size, lbc, lbcold, ubcold, failed = failed, errmsg = errmsg) ! Resize to the requested `size` with the same lower bound, write `array(\@lbcold:ubcold)` to `array(\@lbc:\@lbc-lbcold+ubcold)`, optionally gracefully return upon failure.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  Note that the new elements of the newly allocated `array` are **not** initialized to any particular value on output.<br>
    !>  if `array` is a container of an `allocatable` component, the new elements remain unallocated upon return.<br>
    !>  In such a case, the contents of the new elements of the output `array` is processor dependent, frequently meaningless, and should not be relied upon, even if they seem to have been initialized.<br>
    !>  If the initialization of the new elements with a specific `fill` is necessary, use [setRefilled](@ref pm_arrayRefill::setRefilled) to resize arrays and filling the new elements.<br>
    !>
    !>  \warning
    !>  The condition `all(0 < len/shape(array))` must hold for the corresponding input arguments when the input `size` argument is missing.<br>
    !>  The condition `allocated(array)` must hold (the input `array` must be preallocated) when any or all of the optional input arguments `lbc, lbcold, ubcold` are present or when the `size` argument is missing.<br>
    !>  The condition `all(0 <= size)` must hold for the corresponding input argument.<br>
    !>  The condition `all(lbound(array) <= lbcold .and. lbcold <= ubound(array))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(lbound(array) <= ubcold .and. ubcold <= ubound(array))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(lbound(array) <= lbc)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(lbc + ubcold - lbcold <= lbound(array) + size - 1)` must hold for the corresponding input arguments (i.e., the upper bound(s) of contents cannot overflow the upper bound(s) of the new array).<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  If the input `array` is unallocated, it will be allocated to the requested shape, equivalent to `allocate(array(@size))`.<br>
    !>
    !>  \note
    !>  The sole purpose of this generic interface is to provide a **convenient** but **fast**
    !>  method of resizing allocatable arrays without losing the contents of the array.<br>
    !>  See [pm_arrayResize](@ref pm_arrayResize) for the relevant benchmarks.<br>
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
    !>  \example{setResized}
    !>  \include{lineno} example/pm_arrayResize/setResized/main.F90
    !>  \compilef{setResized}
    !>  \output{setResized}
    !>  \include{lineno} example/pm_arrayResize/setResized/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayResize](@ref test_pm_arrayResize)
    !>
    !>  \todo
    !>  \pvlow
    !>  This generic interface can be extended to arrays of higher ranks than currently supported.<br>
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
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021-2022}
    !>  \desc
    !>  There is an Intel compiler 2020-2022 bug in processing multiple `CHECK_ASSERTION` macros in individual routines of this module in `debug` compile mode.<br>
    !>  The problem does not appear to exist in other compilation modes.<br>
    !>  However, this bug does seem to be related to other similar instances where \ifort cannot tolerate frequent appearance of `use` statements within a single `submodule`.<br>
    !>  \remedy
    !>  For now, the resolution was to remove and replace the checking macros with explicit merged `block` statements.<br>
    !>  A similar problem also was present in the implementation of [pm_quadPack](@ref pm_quadPack).<br>
    !>
    !>  \final{setResized}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! DDDD_D0

    interface setResized

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setResizedDDDD_D0_SK5(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(:,SKG)            , intent(inout) , allocatable       :: array
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setResizedDDDD_D0_SK4(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)            , intent(inout) , allocatable       :: array
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setResizedDDDD_D0_SK3(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)            , intent(inout) , allocatable       :: array
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setResizedDDDD_D0_SK2(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)            , intent(inout) , allocatable       :: array
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setResizedDDDD_D0_SK1(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)            , intent(inout) , allocatable       :: array
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! DDDD_D1

    interface setResized

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setResizedDDDD_D1_SK5(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setResizedDDDD_D1_SK4(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setResizedDDDD_D1_SK3(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setResizedDDDD_D1_SK2(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setResizedDDDD_D1_SK1(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setResizedDDDD_D1_IK5(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setResizedDDDD_D1_IK4(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setResizedDDDD_D1_IK3(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setResizedDDDD_D1_IK2(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setResizedDDDD_D1_IK1(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setResizedDDDD_D1_LK5(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setResizedDDDD_D1_LK4(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setResizedDDDD_D1_LK3(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setResizedDDDD_D1_LK2(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setResizedDDDD_D1_LK1(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setResizedDDDD_D1_CK5(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setResizedDDDD_D1_CK4(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setResizedDDDD_D1_CK3(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setResizedDDDD_D1_CK2(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setResizedDDDD_D1_CK1(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setResizedDDDD_D1_RK5(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setResizedDDDD_D1_RK4(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setResizedDDDD_D1_RK3(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setResizedDDDD_D1_RK2(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setResizedDDDD_D1_RK1(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setResizedDDDD_D1_PSSK5(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setResizedDDDD_D1_PSSK4(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setResizedDDDD_D1_PSSK3(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setResizedDDDD_D1_PSSK2(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setResizedDDDD_D1_PSSK1(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setResizedDDDD_D1_BSSK(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(inout) , allocatable       :: array(:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! DDDD_D2

    interface setResized

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setResizedDDDD_D2_SK5(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setResizedDDDD_D2_SK4(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setResizedDDDD_D2_SK3(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setResizedDDDD_D2_SK2(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setResizedDDDD_D2_SK1(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setResizedDDDD_D2_IK5(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setResizedDDDD_D2_IK4(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setResizedDDDD_D2_IK3(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setResizedDDDD_D2_IK2(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setResizedDDDD_D2_IK1(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setResizedDDDD_D2_LK5(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setResizedDDDD_D2_LK4(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setResizedDDDD_D2_LK3(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setResizedDDDD_D2_LK2(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setResizedDDDD_D2_LK1(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setResizedDDDD_D2_CK5(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setResizedDDDD_D2_CK4(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setResizedDDDD_D2_CK3(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setResizedDDDD_D2_CK2(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setResizedDDDD_D2_CK1(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setResizedDDDD_D2_RK5(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setResizedDDDD_D2_RK4(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setResizedDDDD_D2_RK3(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setResizedDDDD_D2_RK2(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setResizedDDDD_D2_RK1(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setResizedDDDD_D2_PSSK5(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setResizedDDDD_D2_PSSK4(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setResizedDDDD_D2_PSSK3(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setResizedDDDD_D2_PSSK2(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setResizedDDDD_D2_PSSK1(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setResizedDDDD_D2_BSSK(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D2_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(inout) , allocatable       :: array(:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! DDDD_D3

    interface setResized

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setResizedDDDD_D3_SK5(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setResizedDDDD_D3_SK4(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setResizedDDDD_D3_SK3(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setResizedDDDD_D3_SK2(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setResizedDDDD_D3_SK1(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setResizedDDDD_D3_IK5(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setResizedDDDD_D3_IK4(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setResizedDDDD_D3_IK3(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setResizedDDDD_D3_IK2(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setResizedDDDD_D3_IK1(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setResizedDDDD_D3_LK5(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setResizedDDDD_D3_LK4(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setResizedDDDD_D3_LK3(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setResizedDDDD_D3_LK2(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setResizedDDDD_D3_LK1(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setResizedDDDD_D3_CK5(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setResizedDDDD_D3_CK4(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setResizedDDDD_D3_CK3(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setResizedDDDD_D3_CK2(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setResizedDDDD_D3_CK1(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setResizedDDDD_D3_RK5(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setResizedDDDD_D3_RK4(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setResizedDDDD_D3_RK3(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setResizedDDDD_D3_RK2(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setResizedDDDD_D3_RK1(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setResizedDDDD_D3_PSSK5(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setResizedDDDD_D3_PSSK4(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setResizedDDDD_D3_PSSK3(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setResizedDDDD_D3_PSSK2(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setResizedDDDD_D3_PSSK1(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setResizedDDDD_D3_BSSK(array, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedDDDD_D3_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(inout) , allocatable       :: array(:,:,:)
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! SDDD_D0

    interface setResized

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setResizedSDDD_D0_SK5(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setResizedSDDD_D0_SK4(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setResizedSDDD_D0_SK3(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setResizedSDDD_D0_SK2(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setResizedSDDD_D0_SK1(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SDDD_D1

    interface setResized

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setResizedSDDD_D1_SK5(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setResizedSDDD_D1_SK4(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setResizedSDDD_D1_SK3(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setResizedSDDD_D1_SK2(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setResizedSDDD_D1_SK1(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setResizedSDDD_D1_IK5(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setResizedSDDD_D1_IK4(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setResizedSDDD_D1_IK3(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setResizedSDDD_D1_IK2(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setResizedSDDD_D1_IK1(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setResizedSDDD_D1_LK5(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setResizedSDDD_D1_LK4(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setResizedSDDD_D1_LK3(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setResizedSDDD_D1_LK2(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setResizedSDDD_D1_LK1(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setResizedSDDD_D1_CK5(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setResizedSDDD_D1_CK4(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setResizedSDDD_D1_CK3(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setResizedSDDD_D1_CK2(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setResizedSDDD_D1_CK1(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setResizedSDDD_D1_RK5(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setResizedSDDD_D1_RK4(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setResizedSDDD_D1_RK3(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setResizedSDDD_D1_RK2(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setResizedSDDD_D1_RK1(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setResizedSDDD_D1_PSSK5(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setResizedSDDD_D1_PSSK4(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setResizedSDDD_D1_PSSK3(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setResizedSDDD_D1_PSSK2(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setResizedSDDD_D1_PSSK1(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setResizedSDDD_D1_BSSK(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SDDD_D2

    interface setResized

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setResizedSDDD_D2_SK5(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setResizedSDDD_D2_SK4(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setResizedSDDD_D2_SK3(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setResizedSDDD_D2_SK2(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setResizedSDDD_D2_SK1(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setResizedSDDD_D2_IK5(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setResizedSDDD_D2_IK4(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setResizedSDDD_D2_IK3(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setResizedSDDD_D2_IK2(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setResizedSDDD_D2_IK1(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setResizedSDDD_D2_LK5(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setResizedSDDD_D2_LK4(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setResizedSDDD_D2_LK3(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setResizedSDDD_D2_LK2(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setResizedSDDD_D2_LK1(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setResizedSDDD_D2_CK5(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setResizedSDDD_D2_CK4(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setResizedSDDD_D2_CK3(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setResizedSDDD_D2_CK2(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setResizedSDDD_D2_CK1(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setResizedSDDD_D2_RK5(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setResizedSDDD_D2_RK4(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setResizedSDDD_D2_RK3(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setResizedSDDD_D2_RK2(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setResizedSDDD_D2_RK1(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setResizedSDDD_D2_PSSK5(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setResizedSDDD_D2_PSSK4(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setResizedSDDD_D2_PSSK3(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setResizedSDDD_D2_PSSK2(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setResizedSDDD_D2_PSSK1(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setResizedSDDD_D2_BSSK(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D2_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SDDD_D3

    interface setResized

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setResizedSDDD_D3_SK5(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setResizedSDDD_D3_SK4(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setResizedSDDD_D3_SK3(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setResizedSDDD_D3_SK2(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setResizedSDDD_D3_SK1(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setResizedSDDD_D3_IK5(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setResizedSDDD_D3_IK4(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setResizedSDDD_D3_IK3(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setResizedSDDD_D3_IK2(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setResizedSDDD_D3_IK1(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setResizedSDDD_D3_LK5(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setResizedSDDD_D3_LK4(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setResizedSDDD_D3_LK3(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setResizedSDDD_D3_LK2(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setResizedSDDD_D3_LK1(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setResizedSDDD_D3_CK5(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setResizedSDDD_D3_CK4(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setResizedSDDD_D3_CK3(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setResizedSDDD_D3_CK2(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setResizedSDDD_D3_CK1(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setResizedSDDD_D3_RK5(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setResizedSDDD_D3_RK4(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setResizedSDDD_D3_RK3(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setResizedSDDD_D3_RK2(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setResizedSDDD_D3_RK1(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setResizedSDDD_D3_PSSK5(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setResizedSDDD_D3_PSSK4(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setResizedSDDD_D3_PSSK3(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setResizedSDDD_D3_PSSK2(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setResizedSDDD_D3_PSSK1(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setResizedSDDD_D3_BSSK(array, size, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSDDD_D3_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)                        :: size
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! SLDD_D0

    interface setResized

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setResizedSLDD_D0_SK5(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setResizedSLDD_D0_SK4(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setResizedSLDD_D0_SK3(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setResizedSLDD_D0_SK2(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setResizedSLDD_D0_SK1(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SLDD_D1

    interface setResized

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setResizedSLDD_D1_SK5(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setResizedSLDD_D1_SK4(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setResizedSLDD_D1_SK3(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setResizedSLDD_D1_SK2(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setResizedSLDD_D1_SK1(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setResizedSLDD_D1_IK5(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setResizedSLDD_D1_IK4(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setResizedSLDD_D1_IK3(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setResizedSLDD_D1_IK2(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setResizedSLDD_D1_IK1(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setResizedSLDD_D1_LK5(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setResizedSLDD_D1_LK4(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setResizedSLDD_D1_LK3(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setResizedSLDD_D1_LK2(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setResizedSLDD_D1_LK1(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setResizedSLDD_D1_CK5(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setResizedSLDD_D1_CK4(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setResizedSLDD_D1_CK3(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setResizedSLDD_D1_CK2(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setResizedSLDD_D1_CK1(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setResizedSLDD_D1_RK5(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setResizedSLDD_D1_RK4(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setResizedSLDD_D1_RK3(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setResizedSLDD_D1_RK2(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setResizedSLDD_D1_RK1(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setResizedSLDD_D1_PSSK5(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setResizedSLDD_D1_PSSK4(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setResizedSLDD_D1_PSSK3(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setResizedSLDD_D1_PSSK2(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setResizedSLDD_D1_PSSK1(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setResizedSLDD_D1_BSSK(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SLDD_D2

    interface setResized

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setResizedSLDD_D2_SK5(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setResizedSLDD_D2_SK4(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setResizedSLDD_D2_SK3(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setResizedSLDD_D2_SK2(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setResizedSLDD_D2_SK1(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setResizedSLDD_D2_IK5(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setResizedSLDD_D2_IK4(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setResizedSLDD_D2_IK3(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setResizedSLDD_D2_IK2(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setResizedSLDD_D2_IK1(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setResizedSLDD_D2_LK5(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setResizedSLDD_D2_LK4(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setResizedSLDD_D2_LK3(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setResizedSLDD_D2_LK2(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setResizedSLDD_D2_LK1(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setResizedSLDD_D2_CK5(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setResizedSLDD_D2_CK4(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setResizedSLDD_D2_CK3(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setResizedSLDD_D2_CK2(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setResizedSLDD_D2_CK1(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setResizedSLDD_D2_RK5(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setResizedSLDD_D2_RK4(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setResizedSLDD_D2_RK3(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setResizedSLDD_D2_RK2(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setResizedSLDD_D2_RK1(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setResizedSLDD_D2_PSSK5(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setResizedSLDD_D2_PSSK4(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setResizedSLDD_D2_PSSK3(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setResizedSLDD_D2_PSSK2(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setResizedSLDD_D2_PSSK1(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setResizedSLDD_D2_BSSK(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D2_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SLDD_D3

    interface setResized

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setResizedSLDD_D3_SK5(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setResizedSLDD_D3_SK4(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setResizedSLDD_D3_SK3(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setResizedSLDD_D3_SK2(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setResizedSLDD_D3_SK1(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setResizedSLDD_D3_IK5(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setResizedSLDD_D3_IK4(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setResizedSLDD_D3_IK3(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setResizedSLDD_D3_IK2(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setResizedSLDD_D3_IK1(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setResizedSLDD_D3_LK5(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setResizedSLDD_D3_LK4(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setResizedSLDD_D3_LK3(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setResizedSLDD_D3_LK2(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setResizedSLDD_D3_LK1(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setResizedSLDD_D3_CK5(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setResizedSLDD_D3_CK4(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setResizedSLDD_D3_CK3(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setResizedSLDD_D3_CK2(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setResizedSLDD_D3_CK1(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setResizedSLDD_D3_RK5(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setResizedSLDD_D3_RK4(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setResizedSLDD_D3_RK3(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setResizedSLDD_D3_RK2(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setResizedSLDD_D3_RK1(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setResizedSLDD_D3_PSSK5(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setResizedSLDD_D3_PSSK4(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setResizedSLDD_D3_PSSK3(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setResizedSLDD_D3_PSSK2(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setResizedSLDD_D3_PSSK1(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setResizedSLDD_D3_BSSK(array, size, lbc, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLDD_D3_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! SLLU_D0

    interface setResized

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setResizedSLLU_D0_SK5(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setResizedSLLU_D0_SK4(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setResizedSLLU_D0_SK3(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setResizedSLLU_D0_SK2(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setResizedSLLU_D0_SK1(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)            , intent(inout) , allocatable       :: array
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SLLU_D1

    interface setResized

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setResizedSLLU_D1_SK5(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setResizedSLLU_D1_SK4(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setResizedSLLU_D1_SK3(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setResizedSLLU_D1_SK2(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setResizedSLLU_D1_SK1(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setResizedSLLU_D1_IK5(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setResizedSLLU_D1_IK4(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setResizedSLLU_D1_IK3(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setResizedSLLU_D1_IK2(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setResizedSLLU_D1_IK1(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setResizedSLLU_D1_LK5(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setResizedSLLU_D1_LK4(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setResizedSLLU_D1_LK3(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setResizedSLLU_D1_LK2(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setResizedSLLU_D1_LK1(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setResizedSLLU_D1_CK5(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setResizedSLLU_D1_CK4(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setResizedSLLU_D1_CK3(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setResizedSLLU_D1_CK2(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setResizedSLLU_D1_CK1(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setResizedSLLU_D1_RK5(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setResizedSLLU_D1_RK4(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setResizedSLLU_D1_RK3(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setResizedSLLU_D1_RK2(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setResizedSLLU_D1_RK1(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setResizedSLLU_D1_PSSK5(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setResizedSLLU_D1_PSSK4(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setResizedSLLU_D1_PSSK3(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setResizedSLLU_D1_PSSK2(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setResizedSLLU_D1_PSSK1(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setResizedSLLU_D1_BSSK(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(inout) , allocatable       :: array(:)
        integer(IK)                 , intent(in)                        :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SLLU_D2

    interface setResized

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setResizedSLLU_D2_SK5(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setResizedSLLU_D2_SK4(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setResizedSLLU_D2_SK3(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setResizedSLLU_D2_SK2(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setResizedSLLU_D2_SK1(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setResizedSLLU_D2_IK5(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setResizedSLLU_D2_IK4(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setResizedSLLU_D2_IK3(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setResizedSLLU_D2_IK2(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setResizedSLLU_D2_IK1(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setResizedSLLU_D2_LK5(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setResizedSLLU_D2_LK4(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setResizedSLLU_D2_LK3(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setResizedSLLU_D2_LK2(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setResizedSLLU_D2_LK1(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setResizedSLLU_D2_CK5(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setResizedSLLU_D2_CK4(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setResizedSLLU_D2_CK3(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setResizedSLLU_D2_CK2(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setResizedSLLU_D2_CK1(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setResizedSLLU_D2_RK5(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setResizedSLLU_D2_RK4(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setResizedSLLU_D2_RK3(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setResizedSLLU_D2_RK2(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setResizedSLLU_D2_RK1(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setResizedSLLU_D2_PSSK5(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setResizedSLLU_D2_PSSK4(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setResizedSLLU_D2_PSSK3(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setResizedSLLU_D2_PSSK2(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setResizedSLLU_D2_PSSK1(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setResizedSLLU_D2_BSSK(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D2_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(inout) , allocatable       :: array(:,:)
        integer(IK)                 , intent(in)    , dimension(2)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! SLLU_D3

    interface setResized

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setResizedSLLU_D3_SK5(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setResizedSLLU_D3_SK4(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setResizedSLLU_D3_SK3(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setResizedSLLU_D3_SK2(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setResizedSLLU_D3_SK1(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setResizedSLLU_D3_IK5(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setResizedSLLU_D3_IK4(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setResizedSLLU_D3_IK3(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setResizedSLLU_D3_IK2(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setResizedSLLU_D3_IK1(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setResizedSLLU_D3_LK5(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setResizedSLLU_D3_LK4(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setResizedSLLU_D3_LK3(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setResizedSLLU_D3_LK2(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setResizedSLLU_D3_LK1(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setResizedSLLU_D3_CK5(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setResizedSLLU_D3_CK4(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setResizedSLLU_D3_CK3(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setResizedSLLU_D3_CK2(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setResizedSLLU_D3_CK1(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setResizedSLLU_D3_RK5(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setResizedSLLU_D3_RK4(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setResizedSLLU_D3_RK3(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setResizedSLLU_D3_RK2(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setResizedSLLU_D3_RK1(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setResizedSLLU_D3_PSSK5(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setResizedSLLU_D3_PSSK4(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setResizedSLLU_D3_PSSK3(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setResizedSLLU_D3_PSSK2(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setResizedSLLU_D3_PSSK1(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setResizedSLLU_D3_BSSK(array, size, lbc, lbcold, ubcold, failed, errmsg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setResizedSLLU_D3_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(inout) , allocatable       :: array(:,:,:)
        integer(IK)                 , intent(in)    , dimension(3)      :: size, lbc, lbcold, ubcold
        character(*, SK)            , intent(out)   , optional          :: errmsg
        logical(LK)                 , intent(out)   , optional          :: failed
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayResize ! LCOV_EXCL_LINE