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
!>  This module contains procedures and generic interfaces for resizing an input array and centering the original contents of the array in a new array.
!>
!>  \note
!>  The array centering would be perfect if the new and old sizes of the array are both even or both odd.<br>
!>  For example, the contents of an array `"aba"` of size `3` centered in a new size of `5` will be perfectly symmetric
!>  while the same content of length `3` will be slightly asymmetric in a new array with even size, like `4`.
!>
!>  \test
!>  [test_pm_arrayCenter](@ref test_pm_arrayCenter)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayCenter

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_arrayCenter"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a **resized and centered copy** of the input `array` within the newly allocated `arrayCentered`
    !>  and filling the newly added elements (if any) with the user-specified `fill`.<br>
    !>
    !>  \details
    !>  Additionally, if left and right margins `lmsize, rmsize` are specified, resize the input allocatable `array(lb:ub)`
    !>  to `arrayCentered(lb : lb + size + lmsize + rmsize - 1)` while centering the original `array` contents within
    !>  the range `(lb + lmsize : lb + lmsize + size)` in the newly allocated array.
    !>
    !>  \param[inout]   array   :   The input `allocatable` array of shape `(:)` of either <br>
    !>                              <ul>
    !>                                  <li>    type `character` of kind \SKALL, or <br>
    !>                                  <li>    type `integer` of kind \IKALL, or <br>
    !>                                  <li>    type `logical` of kind \LKALL, or <br>
    !>                                  <li>    type `complex` of kind \CKALL, or <br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ul>
    !>                              or,
    !>                              <ul>
    !>                                  <li>    a **scalar assumed-length `character`** of kind \SKALL.<br>
    !>                              </ul>
    !>                              On output, the array will be reallocated to the requested new size, with the same lower bound as before.
    !>  \param[in]      size    :   The input scalar `integer` of default kind \IK representing the length of
    !>                              the section of the array within which the array contents will be centered.<br>
    !>                              If the `lmsize` and `rmsize` input arguments are missing, then `size` represents
    !>                              the total length of the output array, otherwise the length of the output array is `lmsize + size + rmsize`.
    !>  \param[in]      lmsize  :   The input scalar `integer` of default kind \IK representing the size of the left-margin of the output `array`<br>
    !>                              (**optional**, default = `0`. It can be present **only if** the input argument `rmsize` is also present.)
    !>  \param[in]      lmsize  :   The input scalar `integer` of default kind \IK representing the size of the right-margin of the output `array`<br>
    !>                              (**optional**, default = `0`. It can be present **only if** the input argument `lmsize` is also present.)
    !>  \param[in]      fill    :   The input scalar of the same type and kind as the input `array` containing the value to fill the new elements (if any)
    !>                              of `array` surrounding the original `array` contents (excluding the margins). If `array` is of type `character`, then <br>
    !>                              <ol>
    !>                                  <li> the equality `len(fill) == 1` must also hold if `array` is a scalar string.<br>
    !>                                  <li> the equality `len(fill) == len(array)` must also hold if `array` is an array of strings.<br>
    !>                              </ol>
    !>                              (**optional**, if missing, the new elements will remain uninitialized).
    !>  \param[in]      lmfill  :   The input scalar of the same type and kind as the input `array` containing the value to fill the left margin (if any)
    !>                              of newly allocated `array`. If `array` is of type `character`, then <br>
    !>                              <ol>
    !>                                  <li> the equality `len(lmfill) == 1` must also hold if `array` is a scalar string.<br>
    !>                                  <li> the equality `len(lmfill) == len(array)` must also hold if `array` is an array of strings.<br>
    !>                              </ol>
    !>                              (**optional**, if missing, the left-margin will remain uninitialized. It can be present **only if** `lmsize` and `rmsize` are also present.)
    !>  \param[in]      rmfill  :   The input scalar of the same type and kind as the input `array` containing the value to fill the right margin (if any)
    !>                              of newly allocated `array`. If `array` is of type `character`, then <br>
    !>                              <ol>
    !>                                  <li> the equality `len(rmfill) == 1` must also hold if `array` is a scalar string.<br>
    !>                                  <li> the equality `len(rmfill) == len(array)` must also hold if `array` is an array of strings.<br>
    !>                              </ol>
    !>                              (**optional**, if missing, the right-margin will remain uninitialized. It can be present **only if** `lmsize` and `rmsize` are also present.)
    !>
    !>  \return
    !>  `arrayCentered`         :   The output object of the same type, kind, and rank as the input `array`
    !>                              whose size is `size + lmsize + rmsize`, whose contents are the same as the contents of `array`
    !>                              but centered in the range `(1 + lmsize : lmsize + size)` while any new elements are optionally filled
    !>                              with `fill` and the sections `(1 : lmsize)` and `(lmsize + size + 1 : lmsize + size + rmsize)`
    !>                              corresponding to the array margins are filled with `lmfill` and `rmfill` respectively.
    !>
    !>  \interface{getCentered}
    !>  \code{.F90}
    !>
    !>      use pm_arrayCenter, only: getCentered
    !>
    !>      arrayCentered = getCentered(array, size, fill = fill)                                                   ! resize to `(lb:lb+size-1)`, center contents within the new limits, and optionally fill the new elements with `fill`.
    !>      arrayCentered = getCentered(array, size, lmsize, rmsize, fill = fill, lmfill = lmfill, rmfill = rmfill) ! resize to `(lb:lb+size+lmsize+rmsize-1)`, center contents, fill the new elements, left and right margins.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  Note that the new elements of the newly allocated output `arrayCentered` are **not** initialized to any particular value if `fill` is missing.<br>
    !>  Therefore, the contents of the new elements of the output `arrayCentered` is processor dependent, frequently meaningless and should not be relied upon.<br>
    !>  However, if `fill` is specified, then all new elements will be filled with the specified `fill` value.
    !>
    !>  \warning
    !>  Similarly, the margin elements will not be initialized to any particular
    !>  values unless the corresponding `lmfill` or `rmfill` arguments are present.<br>
    !>
    !>  \warning
    !>  When the specified new array size is smaller than the original,
    !>  the corresponding elements of the old array will be also trimmed from the resized shrunk array.<br>
    !>  In such a case, there will be no new elements to initialize their values and `fill` will have no effects on the output.<br>
    !>
    !>  \warning
    !>  The input `size`, `lmsize`, `rmsize` arguments must be non-negative input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  The sole purpose of this generic interface is to provide a **convenient** but **fast** method of generating a resized and centered copy of an array.<br>
    !>  See [pm_arrayResize](@ref pm_arrayResize) and [pm_arrayCenter](@ref pm_arrayCenter) for the relevant benchmarks.<br>
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
    !>  [isFailedGetShellWidth](@ref pm_sysShell::isFailedGetShellWidth)<br>
    !>  [isFailedGetShellHeight](@ref pm_sysShell::isFailedGetShellHeight)<br>
    !>
    !>  \example{getCentered}
    !>  \include{lineno} example/pm_arrayCenter/getCentered/main.F90
    !>  \compilef{getCentered}
    !>  \output{getCentered}
    !>  \include{lineno} example/pm_arrayCenter/getCentered/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayCenter](@ref test_pm_arrayCenter)
    !>
    !>  \todo
    !>  \pmed Two new optional input scalar `lbcold` and `ubcold` arguments can be added to procedures to specify
    !>  a subset of the contents of the original array that has to be kept in the newly allocated centered array.<br>
    !>
    !>  \final{getCentered}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getCentered

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getCenteredAsis_D0_SK5(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        character(size,SKG)                                         :: arrayCentered
    end function
#endif

#if SK4_ENABLED
    PURE module function getCenteredAsis_D0_SK4(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        character(size,SKG)                                         :: arrayCentered
    end function
#endif

#if SK3_ENABLED
    PURE module function getCenteredAsis_D0_SK3(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        character(size,SKG)                                         :: arrayCentered
    end function
#endif

#if SK2_ENABLED
    PURE module function getCenteredAsis_D0_SK2(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        character(size,SKG)                                         :: arrayCentered
    end function
#endif

#if SK1_ENABLED
    PURE module function getCenteredAsis_D0_SK1(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        character(size,SKG)                                         :: arrayCentered
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getCenteredAsis_D1_SK5(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: size
        character(len(array,IK),SKG), intent(in)    , optional      :: fill
        character(len(array,IK),SKG)                                :: arrayCentered(size)
    end function
#endif

#if SK4_ENABLED
    PURE module function getCenteredAsis_D1_SK4(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: size
        character(len(array,IK),SKG), intent(in)    , optional      :: fill
        character(len(array,IK),SKG)                                :: arrayCentered(size)
    end function
#endif

#if SK3_ENABLED
    PURE module function getCenteredAsis_D1_SK3(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: size
        character(len(array,IK),SKG), intent(in)    , optional      :: fill
        character(len(array,IK),SKG)                                :: arrayCentered(size)
    end function
#endif

#if SK2_ENABLED
    PURE module function getCenteredAsis_D1_SK2(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: size
        character(len(array,IK),SKG), intent(in)    , optional      :: fill
        character(len(array,IK),SKG)                                :: arrayCentered(size)
    end function
#endif

#if SK1_ENABLED
    PURE module function getCenteredAsis_D1_SK1(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: size
        character(len(array,IK),SKG), intent(in)    , optional      :: fill
        character(len(array,IK),SKG)                                :: arrayCentered(size)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getCenteredAsis_D1_IK5(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        integer(IKG)                                                :: arrayCentered(size)
    end function
#endif

#if IK4_ENABLED
    PURE module function getCenteredAsis_D1_IK4(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        integer(IKG)                                                :: arrayCentered(size)
    end function
#endif

#if IK3_ENABLED
    PURE module function getCenteredAsis_D1_IK3(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        integer(IKG)                                                :: arrayCentered(size)
    end function
#endif

#if IK2_ENABLED
    PURE module function getCenteredAsis_D1_IK2(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        integer(IKG)                                                :: arrayCentered(size)
    end function
#endif

#if IK1_ENABLED
    PURE module function getCenteredAsis_D1_IK1(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        integer(IKG)                                                :: arrayCentered(size)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getCenteredAsis_D1_LK5(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        logical(LKG)                                                :: arrayCentered(size)
    end function
#endif

#if LK4_ENABLED
    PURE module function getCenteredAsis_D1_LK4(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        logical(LKG)                                                :: arrayCentered(size)
    end function
#endif

#if LK3_ENABLED
    PURE module function getCenteredAsis_D1_LK3(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        logical(LKG)                                                :: arrayCentered(size)
    end function
#endif

#if LK2_ENABLED
    PURE module function getCenteredAsis_D1_LK2(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        logical(LKG)                                                :: arrayCentered(size)
    end function
#endif

#if LK1_ENABLED
    PURE module function getCenteredAsis_D1_LK1(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        logical(LKG)                                                :: arrayCentered(size)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCenteredAsis_D1_CK5(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        complex(CKG)                                                :: arrayCentered(size)
    end function
#endif

#if CK4_ENABLED
    PURE module function getCenteredAsis_D1_CK4(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        complex(CKG)                                                :: arrayCentered(size)
    end function
#endif

#if CK3_ENABLED
    PURE module function getCenteredAsis_D1_CK3(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        complex(CKG)                                                :: arrayCentered(size)
    end function
#endif

#if CK2_ENABLED
    PURE module function getCenteredAsis_D1_CK2(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        complex(CKG)                                                :: arrayCentered(size)
    end function
#endif

#if CK1_ENABLED
    PURE module function getCenteredAsis_D1_CK1(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        complex(CKG)                                                :: arrayCentered(size)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCenteredAsis_D1_RK5(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        real(RKG)                                                   :: arrayCentered(size)
    end function
#endif

#if RK4_ENABLED
    PURE module function getCenteredAsis_D1_RK4(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        real(RKG)                                                   :: arrayCentered(size)
    end function
#endif

#if RK3_ENABLED
    PURE module function getCenteredAsis_D1_RK3(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        real(RKG)                                                   :: arrayCentered(size)
    end function
#endif

#if RK2_ENABLED
    PURE module function getCenteredAsis_D1_RK2(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        real(RKG)                                                   :: arrayCentered(size)
    end function
#endif

#if RK1_ENABLED
    PURE module function getCenteredAsis_D1_RK1(array, size, fill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredAsis_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , optional      :: fill
        integer(IK)                 , intent(in)                    :: size
        real(RKG)                                                   :: arrayCentered(size)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getCenteredMarg_D0_SK5(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        character(size+lmsize+rmsize,SKG)                           :: arrayCentered
    end function
#endif

#if SK4_ENABLED
    PURE module function getCenteredMarg_D0_SK4(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        character(size+lmsize+rmsize,SKG)                           :: arrayCentered
    end function
#endif

#if SK3_ENABLED
    PURE module function getCenteredMarg_D0_SK3(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        character(size+lmsize+rmsize,SKG)                           :: arrayCentered
    end function
#endif

#if SK2_ENABLED
    PURE module function getCenteredMarg_D0_SK2(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        character(size+lmsize+rmsize,SKG)                           :: arrayCentered
    end function
#endif

#if SK1_ENABLED
    PURE module function getCenteredMarg_D0_SK1(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        character(size+lmsize+rmsize,SKG)                           :: arrayCentered
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getCenteredMarg_D1_SK5(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        character(len(array,IK),SKG), intent(in)    , optional      :: fill, lmfill, rmfill
        character(len(array,IK),SKG)                                :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

#if SK4_ENABLED
    PURE module function getCenteredMarg_D1_SK4(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        character(len(array,IK),SKG), intent(in)    , optional      :: fill, lmfill, rmfill
        character(len(array,IK),SKG)                                :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

#if SK3_ENABLED
    PURE module function getCenteredMarg_D1_SK3(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        character(len(array,IK),SKG), intent(in)    , optional      :: fill, lmfill, rmfill
        character(len(array,IK),SKG)                                :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

#if SK2_ENABLED
    PURE module function getCenteredMarg_D1_SK2(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        character(len(array,IK),SKG), intent(in)    , optional      :: fill, lmfill, rmfill
        character(len(array,IK),SKG)                                :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

#if SK1_ENABLED
    PURE module function getCenteredMarg_D1_SK1(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        character(len(array,IK),SKG), intent(in)    , optional      :: fill, lmfill, rmfill
        character(len(array,IK),SKG)                                :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getCenteredMarg_D1_IK5(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        integer(IKG)                                                :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

#if IK4_ENABLED
    PURE module function getCenteredMarg_D1_IK4(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        integer(IKG)                                                :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

#if IK3_ENABLED
    PURE module function getCenteredMarg_D1_IK3(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        integer(IKG)                                                :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

#if IK2_ENABLED
    PURE module function getCenteredMarg_D1_IK2(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        integer(IKG)                                                :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

#if IK1_ENABLED
    PURE module function getCenteredMarg_D1_IK1(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        integer(IKG)                                                :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getCenteredMarg_D1_LK5(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        logical(LKG)                                                :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

#if LK4_ENABLED
    PURE module function getCenteredMarg_D1_LK4(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        logical(LKG)                                                :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

#if LK3_ENABLED
    PURE module function getCenteredMarg_D1_LK3(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        logical(LKG)                                                :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

#if LK2_ENABLED
    PURE module function getCenteredMarg_D1_LK2(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        logical(LKG)                                                :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

#if LK1_ENABLED
    PURE module function getCenteredMarg_D1_LK1(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        logical(LKG)                                                :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#if CK5_ENABLED
    PURE module function getCenteredMarg_D1_CK5(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        complex(CKG)                                                :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

#if CK4_ENABLED
    PURE module function getCenteredMarg_D1_CK4(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        complex(CKG)                                                :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

#if CK3_ENABLED
    PURE module function getCenteredMarg_D1_CK3(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        complex(CKG)                                                :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

#if CK2_ENABLED
    PURE module function getCenteredMarg_D1_CK2(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        complex(CKG)                                                :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

#if CK1_ENABLED
    PURE module function getCenteredMarg_D1_CK1(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        complex(CKG)                                                :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCenteredMarg_D1_RK5(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        real(RKG)                                                   :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

#if RK4_ENABLED
    PURE module function getCenteredMarg_D1_RK4(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        real(RKG)                                                   :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

#if RK3_ENABLED
    PURE module function getCenteredMarg_D1_RK3(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        real(RKG)                                                   :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

#if RK2_ENABLED
    PURE module function getCenteredMarg_D1_RK2(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        real(RKG)                                                   :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

#if RK1_ENABLED
    PURE module function getCenteredMarg_D1_RK1(array, size, lmsize, rmsize, fill, lmfill, rmfill) result(arrayCentered)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCenteredMarg_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                 , intent(in)                    :: size, lmsize, rmsize
        real(RKG)                                                   :: arrayCentered(size+lmsize+rmsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Center the contents of the input `array` within the output `arrayCentered`
    !>  while filling the newly added elements (if any) with the user-specified `fill`.
    !>
    !>  \details
    !>  Additionally, if the left and right margins `lmsize, rmsize` are specified, center the
    !>  input `array` within the section `arrayCentered(lmsize + 1 : size(arrayCentered) - rmsize)`
    !>  and optionally fill the margins with the user-specified `lmfill` and `rmfill`.
    !>
    !>  \param[inout]   arrayCentered   :   The output `contiguous` array-like object of shape `(:)` of the same type and kind as the input `array`,
    !>                                      containing the contents of `array` centered within it, with the requested margins and fillings.<br>
    !>  \param[inout]   array           :   The input `contiguous` array of shape `(:)` of either <br>
    !>                                      <ul>
    !>                                          <li>    type `character` of kind \SKALL of the same length type-parameter as that of `array`, or
    !>                                          <li>    type `integer` of kind \IKALL, or
    !>                                          <li>    type `logical` of kind \LKALL, or
    !>                                          <li>    type `complex` of kind \CKALL, or
    !>                                          <li>    type `real` of kind \RKALL,
    !>                                      </ul>
    !>                                      or,
    !>                                      <ul>
    !>                                          <li>    a **scalar `character`** of kind \SKALL,
    !>                                      </ul>
    !>                                      whose contents will be copied to the center of `arrayCentered` after subtracting
    !>                                      the `lmsize` and `rmsize` elements from the left and the right of `arrayCentered`.
    !>  \param[in]      lmsize          :   The input scalar `integer` of default kind \IK representing the width of the left-margin of the output `arrayCentered`<br>
    !>                                      (**optional**, default = `0`. It must be present **if and only if** the input argument `rmsize` is also present.)
    !>  \param[in]      rmsize          :   The input scalar `integer` of default kind \IK representing the width of the right-margin of the output `arrayCentered`<br>
    !>                                      (**optional**, default = `0`. It must be present **if and only if** the input argument `lmsize` is also present.)
    !>  \param[in]      fill            :   The input scalar of the same type and kind as the input `array` containing the value to fill the new elements (if any)
    !>                                      of `arrayCentered` surrounding the original `array` contents (excluding the margins). If `array` is of type `character`, then <br>
    !>                                      <ol>
    !>                                          <li> the equality `len(fill) == 1` must also hold if `array` is a scalar string.<br>
    !>                                          <li> the equality `len(fill) == len(array)` must also hold if `array` is an array of strings.<br>
    !>                                      </ol>
    !>                                      (**optional**, if missing, the new elements in `arrayCentered` will remain uninitialized.)
    !>  \param[in]      lmfill          :   The input scalar of the same type and kind as the input `array` containing the value to fill the left margin (if any)
    !>                                      of output `arrayCentered`. If `array` is of type `character`, then <br>
    !>                                      <ol>
    !>                                          <li> the equality `len(lmfill) == 1` must also hold if `array` is a scalar string.<br>
    !>                                          <li> the equality `len(lmfill) == len(array)` must also hold if `array` is an array of strings.<br>
    !>                                      </ol>
    !>                                      (**optional**, if missing, the left-margin will remain uninitialized. It can be present **only if** `lmsize` and `rmsize` are also present.)
    !>  \param[in]      rmfill          :   The input scalar of the same type and kind as the input `array` containing the value to fill the right margin (if any)
    !>                                      of output `arrayCentered`. If `array` is of type `character`, then <br>
    !>                                      <ol>
    !>                                          <li> the equality `len(rmfill) == 1` must also hold if `array` is a scalar string.<br>
    !>                                          <li> the equality `len(rmfill) == len(array)` must also hold if `array` is an array of strings.<br>
    !>                                      </ol>
    !>                                      (**optional**, if missing, the right-margin will remain uninitialized. It can be present **only if** `lmsize` and `rmsize` are also present)
    !>
    !>  \interface{setCentered}
    !>  \code{.F90}
    !>
    !>      use pm_arrayCenter, only: setCentered
    !>
    !>      call setCentered(arrayCentered, array, fill = fill)
    !>      call setCentered(arrayCentered, array, lmsize, rmsize, fill = fill, lmfill = lmfill, rmfill = rmfill)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Note that the new elements of the output `arrayCentered` are **not** initialized to any particular value if `fill` is missing.<br>
    !>  Therefore, the contents of the new elements of the output `array` is processor dependent, frequently meaningless and should not be relied upon.<br>
    !>  However, if `fill` is specified, then all new elements will be filled with the specified `fill` value.<br>
    !>
    !>  \warning
    !>  Similarly, the margin elements will not be initialized to any particular values unless the corresponding `lmfill` or `rmfill` arguments are present.<br>
    !>
    !>  \warning
    !>  When the specified new array size is smaller than the original,
    !>  the corresponding elements of the old array will be also trimmed from the resized shrunk array (after taking into account the left and the right margins).<br>
    !>  In such a case, there will be no new elements to initialize their values and `fill` will have no effects on the output.<br>
    !>
    !>  \warning
    !>  The sum of the input `lmsize` and `rmsize` arguments must not be larger than the size of the input `arrayCentered`.<br>
    !>  \vericons
    !>
    !>  \warnpure
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
    !>  \example{setCentered}
    !>  \include{lineno} example/pm_arrayCenter/setCentered/main.F90
    !>  \compilef{setCentered}
    !>  \output{setCentered}
    !>  \include{lineno} example/pm_arrayCenter/setCentered/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayCenter](@ref test_pm_arrayCenter)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{10.3}
    !>  \desc
    !>  \gfortran cannot handle regular allocation for assumed-length allocatable character types and returns the following error.<br>
    !>  \code{.sh}
    !>      Fortran runtime error: Integer overflow when calculating the amount of memory to allocate
    !>  \endcode
    !>  \remedy{1.5}
    !>  The preprocessor conditions bypass this bug.<br>
    !>  \remedy{2.0.0}
    !>  The new interface allocatable output arguments entirely, thus obviating the need to handle this gracefully.<br>
    !>
    !>  \todo
    !>  \pmed
    !>  Two new optional input scalar `lbcold` and `ubcold` arguments can be added to procedures to specify
    !>  a subset of the contents of the original array that has to be kept in the newly allocated centered array.
    !>
    !>  \final{setCentered}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setCentered

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setCenteredAsis_D0_SK5(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)                            , intent(out)                   :: arrayCentered
        character(*,SKG)                            , intent(in)                    :: array
        character(1,SKG)                            , intent(in)    , optional      :: fill
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setCenteredAsis_D0_SK4(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)                            , intent(out)                   :: arrayCentered
        character(*,SKG)                            , intent(in)                    :: array
        character(1,SKG)                            , intent(in)    , optional      :: fill
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setCenteredAsis_D0_SK3(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)                            , intent(out)                   :: arrayCentered
        character(*,SKG)                            , intent(in)                    :: array
        character(1,SKG)                            , intent(in)    , optional      :: fill
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setCenteredAsis_D0_SK2(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)                            , intent(out)                   :: arrayCentered
        character(*,SKG)                            , intent(in)                    :: array
        character(1,SKG)                            , intent(in)    , optional      :: fill
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setCenteredAsis_D0_SK1(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)                            , intent(out)                   :: arrayCentered
        character(*,SKG)                            , intent(in)                    :: array
        character(1,SKG)                            , intent(in)    , optional      :: fill
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setCenteredAsis_D1_SK5(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)                            , intent(out)   , contiguous    :: arrayCentered(:)
        character(*,SKG)                            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG)                , intent(in)    , optional      :: fill
        character(*,SKG)                            , intent(out)   , contiguous    :: array(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setCenteredAsis_D1_SK4(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)                            , intent(out)   , contiguous    :: arrayCentered(:)
        character(*,SKG)                            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG)                , intent(in)    , optional      :: fill
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setCenteredAsis_D1_SK3(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)                            , intent(out)   , contiguous    :: arrayCentered(:)
        character(*,SKG)                            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG)                , intent(in)    , optional      :: fill
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setCenteredAsis_D1_SK2(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)                            , intent(out)   , contiguous    :: arrayCentered(:)
        character(*,SKG)                            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG)                , intent(in)    , optional      :: fill
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setCenteredAsis_D1_SK1(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)                            , intent(out)   , contiguous    :: arrayCentered(:)
        character(*,SKG)                            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG)                , intent(in)    , optional      :: fill
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setCenteredAsis_D1_IK5(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        integer(IKG)                                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                , intent(in)    , optional      :: fill
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setCenteredAsis_D1_IK4(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        integer(IKG)                                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                , intent(in)    , optional      :: fill
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setCenteredAsis_D1_IK3(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        integer(IKG)                                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                , intent(in)    , optional      :: fill
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setCenteredAsis_D1_IK2(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        integer(IKG)                                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                , intent(in)    , optional      :: fill
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setCenteredAsis_D1_IK1(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        integer(IKG)                                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                , intent(in)    , optional      :: fill
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setCenteredAsis_D1_LK5(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        logical(LKG)                                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                , intent(in)    , optional      :: fill
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setCenteredAsis_D1_LK4(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        logical(LKG)                                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                , intent(in)    , optional      :: fill
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setCenteredAsis_D1_LK3(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        logical(LKG)                                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                , intent(in)    , optional      :: fill
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setCenteredAsis_D1_LK2(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        logical(LKG)                                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                , intent(in)    , optional      :: fill
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setCenteredAsis_D1_LK1(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        logical(LKG)                                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                , intent(in)    , optional      :: fill
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCenteredAsis_D1_CK5(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        complex(CKG)                                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                , intent(in)    , optional      :: fill
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCenteredAsis_D1_CK4(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        complex(CKG)                                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                , intent(in)    , optional      :: fill
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCenteredAsis_D1_CK3(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        complex(CKG)                                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                , intent(in)    , optional      :: fill
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCenteredAsis_D1_CK2(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        complex(CKG)                                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                , intent(in)    , optional      :: fill
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCenteredAsis_D1_CK1(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        complex(CKG)                                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                , intent(in)    , optional      :: fill
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCenteredAsis_D1_RK5(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                                   , intent(out)   , contiguous    :: arrayCentered(:)
        real(RKG)                                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                                   , intent(in)    , optional      :: fill
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCenteredAsis_D1_RK4(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                                   , intent(out)   , contiguous    :: arrayCentered(:)
        real(RKG)                                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                                   , intent(in)    , optional      :: fill
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCenteredAsis_D1_RK3(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                                   , intent(out)   , contiguous    :: arrayCentered(:)
        real(RKG)                                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                                   , intent(in)    , optional      :: fill
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCenteredAsis_D1_RK2(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                                   , intent(out)   , contiguous    :: arrayCentered(:)
        real(RKG)                                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                                   , intent(in)    , optional      :: fill
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCenteredAsis_D1_RK1(arrayCentered, array, fill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredAsis_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                                   , intent(out)   , contiguous    :: arrayCentered(:)
        real(RKG)                                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                                   , intent(in)    , optional      :: fill
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setCenteredMarg_D0_SK5(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)                            , intent(out)                   :: arrayCentered
        character(*,SKG)                            , intent(in)                    :: array
        character(1,SKG)                            , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setCenteredMarg_D0_SK4(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)                            , intent(out)                   :: arrayCentered
        character(*,SKG)                            , intent(in)                    :: array
        character(1,SKG)                            , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setCenteredMarg_D0_SK3(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)                            , intent(out)                   :: arrayCentered
        character(*,SKG)                            , intent(in)                    :: array
        character(1,SKG)                            , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setCenteredMarg_D0_SK2(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)                            , intent(out)                   :: arrayCentered
        character(*,SKG)                            , intent(in)                    :: array
        character(1,SKG)                            , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setCenteredMarg_D0_SK1(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)                            , intent(out)                   :: arrayCentered
        character(*,SKG)                            , intent(in)                    :: array
        character(1,SKG)                            , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setCenteredMarg_D1_SK5(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)                            , intent(out)   , contiguous    :: arrayCentered(:)
        character(*,SKG)                            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG)                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setCenteredMarg_D1_SK4(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)                            , intent(out)   , contiguous    :: arrayCentered(:)
        character(*,SKG)                            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG)                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setCenteredMarg_D1_SK3(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)                            , intent(out)   , contiguous    :: arrayCentered(:)
        character(*,SKG)                            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG)                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setCenteredMarg_D1_SK2(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)                            , intent(out)   , contiguous    :: arrayCentered(:)
        character(*,SKG)                            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG)                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setCenteredMarg_D1_SK1(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)                            , intent(out)   , contiguous    :: arrayCentered(:)
        character(*,SKG)                            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG)                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setCenteredMarg_D1_IK5(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        integer(IKG)                                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setCenteredMarg_D1_IK4(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        integer(IKG)                                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setCenteredMarg_D1_IK3(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        integer(IKG)                                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setCenteredMarg_D1_IK2(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        integer(IKG)                                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setCenteredMarg_D1_IK1(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        integer(IKG)                                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setCenteredMarg_D1_LK5(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        logical(LKG)                                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setCenteredMarg_D1_LK4(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        logical(LKG)                                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setCenteredMarg_D1_LK3(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        logical(LKG)                                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setCenteredMarg_D1_LK2(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        logical(LKG)                                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setCenteredMarg_D1_LK1(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        logical(LKG)                                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCenteredMarg_D1_CK5(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        complex(CKG)                                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCenteredMarg_D1_CK4(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        complex(CKG)                                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCenteredMarg_D1_CK3(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        complex(CKG)                                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCenteredMarg_D1_CK2(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        complex(CKG)                                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCenteredMarg_D1_CK1(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                                , intent(out)   , contiguous    :: arrayCentered(:)
        complex(CKG)                                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCenteredMarg_D1_RK5(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                                   , intent(out)   , contiguous    :: arrayCentered(:)
        real(RKG)                                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                                   , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCenteredMarg_D1_RK4(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                                   , intent(out)   , contiguous    :: arrayCentered(:)
        real(RKG)                                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                                   , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCenteredMarg_D1_RK3(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                                   , intent(out)   , contiguous    :: arrayCentered(:)
        real(RKG)                                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                                   , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCenteredMarg_D1_RK2(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                                   , intent(out)   , contiguous    :: arrayCentered(:)
        real(RKG)                                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                                   , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCenteredMarg_D1_RK1(arrayCentered, array, lmsize, rmsize, fill, lmfill, rmfill)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenteredMarg_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                                   , intent(out)   , contiguous    :: arrayCentered(:)
        real(RKG)                                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                                   , intent(in)    , optional      :: fill, lmfill, rmfill
        integer(IK)                                 , intent(in)                    :: lmsize, rmsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayCenter ! LCOV_EXCL_LINE