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
!>  This module contains procedures and generic interfaces for resizing an input array and padding them with symbols on the left or right.
!>
!>  \test
!>  [test_pm_arrayPad](@ref test_pm_arrayPad)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayPad

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_arrayPad"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type :: padb_type; end type ! double padding.
    type :: padl_type; end type
    type :: padr_type; end type
    type(padb_type), parameter :: padb = padb_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: padb
#endif
    type(padl_type), parameter :: padl = padl_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: padl
#endif
    type(padr_type), parameter :: padr = padr_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: padr
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate a resized copy of the input `array` by padding it to the left and right with the
    !>  requested paddings and optionally adding margins to the left and the right of the padded `array`
    !>  optionally filled with the corresponding fills.
    !>
    !>  \param[in]  array   :   The input `contiguous` array of shape `(:)` of either <br>
    !>                          <ul>
    !>                              <li>    type `character` of kind \SKALL, or <br>
    !>                              <li>    type `integer` of kind \IKALL, or <br>
    !>                              <li>    type `logical` of kind \LKALL, or <br>
    !>                              <li>    type `complex` of kind \CKALL, or <br>
    !>                              <li>    type `real` of kind \RKALL,<br>
    !>                          </ul>
    !>                          or,
    !>                          <ul>
    !>                              <li>    a **scalar assumed-length `character`** of kind \SKALL.<br>
    !>                          </ul>
    !>                          If the `lmsize` and `rmsize` input arguments are missing, then the size of the output `arrayPadded`
    !>                          is `lenArray + lpsize + rpsize`, otherwise the length of the output `arrayPadded` is `lenArray + lmsize + lpsize + rpsize + rmsize`.
    !>  \param[in]  lpsize  :   The input scalar of type `integer` of default kind \IK representing the number of `lpfill` to add to the left of the array.
    !>  \param[in]  rpsize  :   The input scalar of type `integer` of default kind \IK representing the number of `rpfill` to add to the right of the array.
    !>  \param[in]  lpfill  :   The input scalar of the same type and kind as the input `array` containing the value to fill the left padding (if any)
    !>                          of the output `arrayPadded`. If `array` is of type `character`, then <br>
    !>                              -# the equality `len(lpfill) == 1` must also hold if `array` is a scalar string.<br>
    !>                              -# the equality `len(lpfill) == len(array)` must also hold if `array` is an array of strings.<br>
    !>  \param[in]  rpfill  :   The input scalar of the same type and kind as the input `array` containing the value to fill the right padding (if any)
    !>                          of the output `arrayPadded`. If `array` is of type `character`, then <br>
    !>                              -# the equality `len(rpfill) == 1` must also hold if `array` is a scalar string.<br>
    !>                              -# the equality `len(rpfill) == len(array)` must also hold if `array` is an array of strings.<br>
    !>  \param[in]  lmsize  :   The input scalar `integer` of default kind \IK representing the size of the left-margin of the output `arrayPadded`<br>
    !>                          (**optional**, default = `0`. It can be present **only if** the input argument `rmsize` is also present.)
    !>  \param[in]  rmsize  :   The input scalar `integer` of default kind \IK representing the size of the right-margin of the output `arrayPadded`<br>
    !>                          (**optional**, default = `0`. It can be present **only if** the input argument `lmsize` is also present.)
    !>  \param[in]  lmfill  :   The input scalar of the same type and kind as the input `array` containing the value to fill the left margin (if any)
    !>                          of newly allocated `array`. If `array` is of type `character`, then <br>
    !>                          <ol>
    !>                              <li> the equality `len(fill) == 1` must also hold if `array` is a scalar string.<br>
    !>                              <li> the equality `len(fill) == len(array)` must also hold if `array` is an array of strings.<br>
    !>                          </ol>
    !>                          (**optional**, if missing, the left-margin will remain uninitialized. It can be present **only if** `lmsize` and `rmsize` are also present.)
    !>  \param[in]  rmfill  :   The input scalar of the same type and kind as the input `array` containing the value to fill the right margin (if any)
    !>                          of newly allocated `array`. If `array` is of type `character`, then <br>
    !>                          <ol>
    !>                              <li> the equality `len(fill) == 1` must also hold if `array` is a scalar string.
    !>                              <li> the equality `len(fill) == len(array)` must also hold if `array` is an array of strings.
    !>                          </ol>
    !>                          (**optional**, if missing, the right-margin will remain uninitialized. It can be present **only if** `lmsize` and `rmsize` are also present.)
    !>
    !>  \return
    !>  `arrayPadded`       :   The output object of the same type, kind, and rank as the input `array` with the same lower bound,
    !>                          whose size is `lenArray + lpsize + rpsize + lmsize + rmsize`, whose contents are the same as the contents of `array`
    !>                          but padded to the left and right, optionally with the specified left and right margins.
    !>
    !>  \interface{getPadded}
    !>  \code{.F90}
    !>
    !>      use pm_arrayPad, only: getPadded
    !>
    !>      arrayPadded = getPadded(array, lpsize, rpsize, lpfill, rpfill)
    !>      arrayPadded = getPadded(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill = lmfill, rmfill = rmfill)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The margin elements will not be initialized to any particular values
    !>  unless the corresponding `lmfill` or `rmfill` arguments are present.
    !>
    !>  \warning
    !>  The input `lpsize`, `rpsize`, `lmsize`, `rmsize` arguments must be non-negative input arguments.<br>
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
    !>  \example{getPadded}
    !>  \include{lineno} example/pm_arrayPad/getPadded/main.F90
    !>  \compilef{getPadded}
    !>  \output{getPadded}
    !>  \include{lineno} example/pm_arrayPad/getPadded/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayPad](@ref test_pm_arrayPad)
    !>
    !>  \todo
    !>  \pmed
    !>  Two new optional input scalar `lbcold` and `ubcold` arguments can be added to procedures to specify
    !>  a subset of the contents of the original array that has to be kept in the newly allocated padded array.
    !>
    !>  \final{getPadded}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getPadded

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getPaddedAsisSB_D0_SK5(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        character(1,SKG)            , intent(in)                    :: lpfill, rpfill
        character(len(array,IK)+lpsize+rpsize,SKG)                  :: arrayPadded
    end function
#endif

#if SK4_ENABLED
    PURE module function getPaddedAsisSB_D0_SK4(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        character(1,SKG)            , intent(in)                    :: lpfill, rpfill
        character(len(array,IK)+lpsize+rpsize,SKG)                  :: arrayPadded
    end function
#endif

#if SK3_ENABLED
    PURE module function getPaddedAsisSB_D0_SK3(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        character(1,SKG)            , intent(in)                    :: lpfill, rpfill
        character(len(array,IK)+lpsize+rpsize,SKG)                  :: arrayPadded
    end function
#endif

#if SK2_ENABLED
    PURE module function getPaddedAsisSB_D0_SK2(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        character(1,SKG)            , intent(in)                    :: lpfill, rpfill
        character(len(array,IK)+lpsize+rpsize,SKG)                  :: arrayPadded
    end function
#endif

#if SK1_ENABLED
    PURE module function getPaddedAsisSB_D0_SK1(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        character(1,SKG)            , intent(in)                    :: lpfill, rpfill
        character(len(array,IK)+lpsize+rpsize,SKG)                  :: arrayPadded
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getPaddedAsisSB_D1_SK5(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        character(len(array,IK),SKG), intent(in)                    :: lpfill, rpfill
        character(len(array),SKG)                                   :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

#if SK4_ENABLED
    PURE module function getPaddedAsisSB_D1_SK4(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        character(len(array,IK),SKG), intent(in)                    :: lpfill, rpfill
        character(len(array),SKG)                                   :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

#if SK3_ENABLED
    PURE module function getPaddedAsisSB_D1_SK3(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        character(len(array,IK),SKG), intent(in)                    :: lpfill, rpfill
        character(len(array),SKG)                                   :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

#if SK2_ENABLED
    PURE module function getPaddedAsisSB_D1_SK2(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        character(len(array,IK),SKG), intent(in)                    :: lpfill, rpfill
        character(len(array),SKG)                                   :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

#if SK1_ENABLED
    PURE module function getPaddedAsisSB_D1_SK1(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        character(len(array,IK),SKG), intent(in)                    :: lpfill, rpfill
        character(len(array),SKG)                                   :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getPaddedAsisSB_D1_IK5(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        integer(IKG)                , intent(in)                    :: lpfill, rpfill
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

#if IK4_ENABLED
    PURE module function getPaddedAsisSB_D1_IK4(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        integer(IKG)                , intent(in)                    :: lpfill, rpfill
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

#if IK3_ENABLED
    PURE module function getPaddedAsisSB_D1_IK3(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        integer(IKG)                , intent(in)                    :: lpfill, rpfill
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

#if IK2_ENABLED
    PURE module function getPaddedAsisSB_D1_IK2(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        integer(IKG)                , intent(in)                    :: lpfill, rpfill
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

#if IK1_ENABLED
    PURE module function getPaddedAsisSB_D1_IK1(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        integer(IKG)                , intent(in)                    :: lpfill, rpfill
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getPaddedAsisSB_D1_LK5(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        logical(LKG)                , intent(in)                    :: lpfill, rpfill
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

#if LK4_ENABLED
    PURE module function getPaddedAsisSB_D1_LK4(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        logical(LKG)                , intent(in)                    :: lpfill, rpfill
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

#if LK3_ENABLED
    PURE module function getPaddedAsisSB_D1_LK3(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        logical(LKG)                , intent(in)                    :: lpfill, rpfill
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

#if LK2_ENABLED
    PURE module function getPaddedAsisSB_D1_LK2(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        logical(LKG)                , intent(in)                    :: lpfill, rpfill
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

#if LK1_ENABLED
    PURE module function getPaddedAsisSB_D1_LK1(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        logical(LKG)                , intent(in)                    :: lpfill, rpfill
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getPaddedAsisSB_D1_CK5(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        complex(CKG)                , intent(in)                    :: lpfill, rpfill
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

#if CK4_ENABLED
    PURE module function getPaddedAsisSB_D1_CK4(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        complex(CKG)                , intent(in)                    :: lpfill, rpfill
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

#if CK3_ENABLED
    PURE module function getPaddedAsisSB_D1_CK3(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        complex(CKG)                , intent(in)                    :: lpfill, rpfill
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

#if CK2_ENABLED
    PURE module function getPaddedAsisSB_D1_CK2(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        complex(CKG)                , intent(in)                    :: lpfill, rpfill
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

#if CK1_ENABLED
    PURE module function getPaddedAsisSB_D1_CK1(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        complex(CKG)                , intent(in)                    :: lpfill, rpfill
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPaddedAsisSB_D1_RK5(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        real(RKG)                   , intent(in)                    :: lpfill, rpfill
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

#if RK4_ENABLED
    PURE module function getPaddedAsisSB_D1_RK4(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        real(RKG)                   , intent(in)                    :: lpfill, rpfill
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

#if RK3_ENABLED
    PURE module function getPaddedAsisSB_D1_RK3(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        real(RKG)                   , intent(in)                    :: lpfill, rpfill
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

#if RK2_ENABLED
    PURE module function getPaddedAsisSB_D1_RK2(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        real(RKG)                   , intent(in)                    :: lpfill, rpfill
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

#if RK1_ENABLED
    PURE module function getPaddedAsisSB_D1_RK1(array, lpsize, rpsize, lpfill, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSB_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        real(RKG)                   , intent(in)                    :: lpfill, rpfill
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+lpsize+rpsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getPaddedMargSB_D0_SK5(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)                    :: lpfill, rpfill
        character(1,SKG)            , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        character(len(array,IK)+lpsize+rpsize+lmsize+rmsize,SKG)    :: arrayPadded
    end function
#endif

#if SK4_ENABLED
    PURE module function getPaddedMargSB_D0_SK4(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)                    :: lpfill, rpfill
        character(1,SKG)            , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        character(len(array,IK)+lpsize+rpsize+lmsize+rmsize,SKG)    :: arrayPadded
    end function
#endif

#if SK3_ENABLED
    PURE module function getPaddedMargSB_D0_SK3(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)                    :: lpfill, rpfill
        character(1,SKG)            , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        character(len(array,IK)+lpsize+rpsize+lmsize+rmsize,SKG)    :: arrayPadded
    end function
#endif

#if SK2_ENABLED
    PURE module function getPaddedMargSB_D0_SK2(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)                    :: lpfill, rpfill
        character(1,SKG)            , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        character(len(array,IK)+lpsize+rpsize+lmsize+rmsize,SKG)    :: arrayPadded
    end function
#endif

#if SK1_ENABLED
    PURE module function getPaddedMargSB_D0_SK1(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)                    :: lpfill, rpfill
        character(1,SKG)            , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        character(len(array,IK)+lpsize+rpsize+lmsize+rmsize,SKG)    :: arrayPadded
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getPaddedMargSB_D1_SK5(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: lpfill, rpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        character(len(array),SKG)                                   :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

#if SK4_ENABLED
    PURE module function getPaddedMargSB_D1_SK4(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: lpfill, rpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        character(len(array),SKG)                                   :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

#if SK3_ENABLED
    PURE module function getPaddedMargSB_D1_SK3(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: lpfill, rpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        character(len(array),SKG)                                   :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

#if SK2_ENABLED
    PURE module function getPaddedMargSB_D1_SK2(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: lpfill, rpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        character(len(array),SKG)                                   :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

#if SK1_ENABLED
    PURE module function getPaddedMargSB_D1_SK1(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: lpfill, rpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        character(len(array),SKG)                                   :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getPaddedMargSB_D1_IK5(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: lpfill, rpfill
        integer(IKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

#if IK4_ENABLED
    PURE module function getPaddedMargSB_D1_IK4(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: lpfill, rpfill
        integer(IKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

#if IK3_ENABLED
    PURE module function getPaddedMargSB_D1_IK3(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: lpfill, rpfill
        integer(IKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

#if IK2_ENABLED
    PURE module function getPaddedMargSB_D1_IK2(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: lpfill, rpfill
        integer(IKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

#if IK1_ENABLED
    PURE module function getPaddedMargSB_D1_IK1(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: lpfill, rpfill
        integer(IKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getPaddedMargSB_D1_LK5(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: lpfill, rpfill
        logical(LKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

#if LK4_ENABLED
    PURE module function getPaddedMargSB_D1_LK4(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: lpfill, rpfill
        logical(LKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

#if LK3_ENABLED
    PURE module function getPaddedMargSB_D1_LK3(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: lpfill, rpfill
        logical(LKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

#if LK2_ENABLED
    PURE module function getPaddedMargSB_D1_LK2(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: lpfill, rpfill
        logical(LKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

#if LK1_ENABLED
    PURE module function getPaddedMargSB_D1_LK1(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: lpfill, rpfill
        logical(LKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getPaddedMargSB_D1_CK5(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: lpfill, rpfill
        complex(CKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

#if CK4_ENABLED
    PURE module function getPaddedMargSB_D1_CK4(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: lpfill, rpfill
        complex(CKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

#if CK3_ENABLED
    PURE module function getPaddedMargSB_D1_CK3(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: lpfill, rpfill
        complex(CKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

#if CK2_ENABLED
    PURE module function getPaddedMargSB_D1_CK2(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: lpfill, rpfill
        complex(CKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

#if CK1_ENABLED
    PURE module function getPaddedMargSB_D1_CK1(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: lpfill, rpfill
        complex(CKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPaddedMargSB_D1_RK5(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: lpfill, rpfill
        real(RKG)                   , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

#if RK4_ENABLED
    PURE module function getPaddedMargSB_D1_RK4(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: lpfill, rpfill
        real(RKG)                   , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

#if RK3_ENABLED
    PURE module function getPaddedMargSB_D1_RK3(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: lpfill, rpfill
        real(RKG)                   , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

#if RK2_ENABLED
    PURE module function getPaddedMargSB_D1_RK2(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: lpfill, rpfill
        real(RKG)                   , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

#if RK1_ENABLED
    PURE module function getPaddedMargSB_D1_RK1(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSB_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: lpfill, rpfill
        real(RKG)                   , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+lpsize+rpsize+lmsize+rmsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Resize the input `array` by padding it to the left and right with the
    !>  requested paddings and optionally adding margins to the left and the right of the padded `array`
    !>  optionally filled with the corresponding fills.
    !>
    !>  \param[inout]   array   :   The input/output `allocatable` array of shape `(:)` of either <br>
    !>                              <ul>
    !>                                  <li>    type `character` of kind \SKALL, or <br>
    !>                                  <li>    type `integer` of kind \IKALL, or <br>
    !>                                  <li>    type `logical` of kind \LKALL, or <br>
    !>                                  <li>    type `complex` of kind \CKALL, or <br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ul>
    !>                              or,
    !>                              <ul>
    !>                                  <li>    a **scalar allocatable `character`** of default kind \SK.<br>
    !>                              </ul>
    !>                              If the `lmsize` and `rmsize` input arguments are missing, then the size of the output `array`
    !>                              is `lenArray + lpsize + rpsize`, otherwise the length of the output `array` is `lenArray + lmsize + lpsize + rpsize + rmsize`.<br>
    !>                              On output, the array will be resized (with the same lower bound as before) and padded and marginalized accordingly.
    !>  \param[in]      lpsize  :   The input scalar of type `integer` of default kind \IK representing the number of `lpfill` to add to the left of the array.
    !>  \param[in]      rpsize  :   The input scalar of type `integer` of default kind \IK representing the number of `rpfill` to add to the right of the array.
    !>  \param[in]      lpfill  :   The input scalar of the same type and kind as the input `array` containing the value to fill the left padding (if any)
    !>                              of the output `array`. If `array` is of type `character`, then <br>
    !>                                  -# the equality `len(lpfill) == 1` must also hold if `array` is a scalar string.<br>
    !>                                  -# the equality `len(lpfill) == len(array)` must also hold if `array` is an array of strings.<br>
    !>  \param[in]      rpfill  :   The input scalar of the same type and kind as the input `array` containing the value to fill the right padding (if any)
    !>                              of the output `array`. If `array` is of type `character`, then <br>
    !>                                  -# the equality `len(rpfill) == 1` must also hold if `array` is a scalar string.<br>
    !>                                  -# the equality `len(rpfill) == len(array)` must also hold if `array` is an array of strings.<br>
    !>  \param[in]      lmsize  :   The input scalar `integer` of default kind \IK representing the size of the left-margin of the output `array`<br>
    !>                              (**optional**, default = `0`. It can be present **only if** the input argument `rmsize` is also present.)
    !>  \param[in]      rmsize  :   The input scalar `integer` of default kind \IK representing the size of the right-margin of the output `array`<br>
    !>                              (**optional**, default = `0`. It can be present **only if** the input argument `lmsize` is also present.)
    !>  \param[in]      lmfill  :   The input scalar of the same type and kind as the input `array` containing the value to fill the left margin (if any)
    !>                              of newly allocated `array`. If `array` is of type `character`, then <br>
    !>                              <ol>
    !>                                  <li> the equality `len(fill) == 1` must also hold if `array` is a scalar string.<br>
    !>                                  <li> the equality `len(fill) == len(array)` must also hold if `array` is an array of strings.<br>
    !>                              </ol>
    !>                              (**optional**, if missing, the left-margin will remain uninitialized. It can be present **only if** `lmsize` and `rmsize` are also present.)
    !>  \param[in]      rmfill  :   The input scalar of the same type and kind as the input `array` containing the value to fill the right margin (if any)
    !>                              of newly allocated `array`. If `array` is of type `character`, then <br>
    !>                              <ol>
    !>                                  <li> the equality `len(fill) == 1` must also hold if `array` is a scalar string.
    !>                                  <li> the equality `len(fill) == len(array)` must also hold if `array` is an array of strings.
    !>                              </ol>
    !>                              (**optional**, if missing, the right-margin will remain uninitialized. It can be present **only if** `lmsize` and `rmsize` are also present.)
    !>  \param[in]      failed  :   The input scalar `logical` of default kind \LK that is `.true.` if the requested array resizing and padding is successful<br>
    !>                              (**optional**, if missing and an error occurs, the processor dictates the program behavior).
    !>
    !>  \interface{setPadded}
    !>  \code{.F90}
    !>
    !>      use pm_arrayPad, only: setPadded
    !>
    !>      call setPadded(array, lpsize, rpsize, lpfill, rpfill, failed = failed)
    !>      call setPadded(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill = lmfill, rmfill = rmfill, failed = failed)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The margin elements will not be initialized to any particular values
    !>  unless the corresponding `lmfill` or `rmfill` arguments are present.
    !>
    !>  \warning
    !>  The input `lpsize`, `rpsize`, `lmsize`, `rmsize` arguments must be non-negative input arguments.<br>
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
    !>  \example{setPadded}
    !>  \include{lineno} example/pm_arrayPad/setPadded/main.F90
    !>  \compilef{setPadded}
    !>  \output{setPadded}
    !>  \include{lineno} example/pm_arrayPad/setPadded/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayPad](@ref test_pm_arrayPad)
    !>
    !>  \todo
    !>  \pmed Two new optional input scalar `lbcold` and `ubcold` arguments can be added to procedures to specify
    !>  a subset of the contents of the original array that has to be kept in the newly allocated padded array.
    !>
    !>  \final{setPadded}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setPadded

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setPaddedAsisSB_D0_SK5(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(:,SKG)            , intent(inout) , allocatable   :: array
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        character(1,SKG)            , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setPaddedAsisSB_D0_SK4(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)            , intent(inout) , allocatable   :: array
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        character(1,SKG)            , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setPaddedAsisSB_D0_SK3(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)            , intent(inout) , allocatable   :: array
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        character(1,SKG)            , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setPaddedAsisSB_D0_SK2(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)            , intent(inout) , allocatable   :: array
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        character(1,SKG)            , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setPaddedAsisSB_D0_SK1(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)            , intent(inout) , allocatable   :: array
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        character(1,SKG)            , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_SK5(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        character(len(array,IK),SKG), intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_SK4(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        character(len(array,IK),SKG), intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_SK3(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        character(len(array,IK),SKG), intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_SK2(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        character(len(array,IK),SKG), intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_SK1(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        character(len(array,IK),SKG), intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_IK5(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        integer(IKG)                , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_IK4(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        integer(IKG)                , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_IK3(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        integer(IKG)                , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_IK2(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        integer(IKG)                , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_IK1(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        integer(IKG)                , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_LK5(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        logical(LKG)                , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_LK4(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        logical(LKG)                , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_LK3(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        logical(LKG)                , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_LK2(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        logical(LKG)                , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_LK1(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        logical(LKG)                , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_CK5(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        complex(CKG)                , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_CK4(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        complex(CKG)                , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_CK3(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        complex(CKG)                , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_CK2(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        complex(CKG)                , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_CK1(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        complex(CKG)                , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_RK5(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        real(RKG)                   , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_RK4(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        real(RKG)                   , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_RK3(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        real(RKG)                   , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_RK2(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        real(RKG)                   , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPaddedAsisSB_D1_RK1(array, lpsize, rpsize, lpfill, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSB_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize, rpsize
        real(RKG)                   , intent(in)                    :: lpfill, rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setPaddedMargSB_D0_SK5(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(:,SKG)            , intent(inout) , allocatable   :: array
        character(1,SKG)            , intent(in)                    :: lpfill, rpfill
        character(1,SKG)            , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setPaddedMargSB_D0_SK4(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)            , intent(inout) , allocatable   :: array
        character(1,SKG)            , intent(in)                    :: lpfill, rpfill
        character(1,SKG)            , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setPaddedMargSB_D0_SK3(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)            , intent(inout) , allocatable   :: array
        character(1,SKG)            , intent(in)                    :: lpfill, rpfill
        character(1,SKG)            , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setPaddedMargSB_D0_SK2(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)            , intent(inout) , allocatable   :: array
        character(1,SKG)            , intent(in)                    :: lpfill, rpfill
        character(1,SKG)            , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setPaddedMargSB_D0_SK1(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)            , intent(inout) , allocatable   :: array
        character(1,SKG)            , intent(in)                    :: lpfill, rpfill
        character(1,SKG)            , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setPaddedMargSB_D1_SK5(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: lpfill, rpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setPaddedMargSB_D1_SK4(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: lpfill, rpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setPaddedMargSB_D1_SK3(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: lpfill, rpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setPaddedMargSB_D1_SK2(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: lpfill, rpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setPaddedMargSB_D1_SK1(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: lpfill, rpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setPaddedMargSB_D1_IK5(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IKG)                , intent(in)                    :: lpfill, rpfill
        integer(IKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setPaddedMargSB_D1_IK4(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IKG)                , intent(in)                    :: lpfill, rpfill
        integer(IKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setPaddedMargSB_D1_IK3(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IKG)                , intent(in)                    :: lpfill, rpfill
        integer(IKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setPaddedMargSB_D1_IK2(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IKG)                , intent(in)                    :: lpfill, rpfill
        integer(IKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setPaddedMargSB_D1_IK1(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IKG)                , intent(in)                    :: lpfill, rpfill
        integer(IKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setPaddedMargSB_D1_LK5(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        logical(LKG)                , intent(in)                    :: lpfill, rpfill
        logical(LKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setPaddedMargSB_D1_LK4(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        logical(LKG)                , intent(in)                    :: lpfill, rpfill
        logical(LKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setPaddedMargSB_D1_LK3(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        logical(LKG)                , intent(in)                    :: lpfill, rpfill
        logical(LKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setPaddedMargSB_D1_LK2(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        logical(LKG)                , intent(in)                    :: lpfill, rpfill
        logical(LKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setPaddedMargSB_D1_LK1(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        logical(LKG)                , intent(in)                    :: lpfill, rpfill
        logical(LKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPaddedMargSB_D1_CK5(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        complex(CKG)                , intent(in)                    :: lpfill, rpfill
        complex(CKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPaddedMargSB_D1_CK4(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        complex(CKG)                , intent(in)                    :: lpfill, rpfill
        complex(CKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPaddedMargSB_D1_CK3(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        complex(CKG)                , intent(in)                    :: lpfill, rpfill
        complex(CKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPaddedMargSB_D1_CK2(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        complex(CKG)                , intent(in)                    :: lpfill, rpfill
        complex(CKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPaddedMargSB_D1_CK1(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        complex(CKG)                , intent(in)                    :: lpfill, rpfill
        complex(CKG)                , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPaddedMargSB_D1_RK5(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        real(RKG)                   , intent(in)                    :: lpfill, rpfill
        real(RKG)                   , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPaddedMargSB_D1_RK4(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        real(RKG)                   , intent(in)                    :: lpfill, rpfill
        real(RKG)                   , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPaddedMargSB_D1_RK3(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        real(RKG)                   , intent(in)                    :: lpfill, rpfill
        real(RKG)                   , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPaddedMargSB_D1_RK2(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        real(RKG)                   , intent(in)                    :: lpfill, rpfill
        real(RKG)                   , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPaddedMargSB_D1_RK1(array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSB_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        real(RKG)                   , intent(in)                    :: lpfill, rpfill
        real(RKG)                   , intent(in)    , optional      :: lmfill, rmfill
        integer(IK)                 , intent(in)                    :: lpsize, rpsize, lmsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate a resized copy of the input `array` by padding it to the left with the
    !>  requested paddings and optionally adding margins to the left of the padded `array`
    !>  optionally filled with the corresponding fills.
    !>
    !>  \param[in]  array   :   The input `contiguous` array of shape `(:)` of either <br>
    !>                          <ul>
    !>                              <li>    type `character` of kind \SKALL, or <br>
    !>                              <li>    type `integer` of kind \IKALL, or <br>
    !>                              <li>    type `logical` of kind \LKALL, or <br>
    !>                              <li>    type `complex` of kind \CKALL, or <br>
    !>                              <li>    type `real` of kind \RKALL,<br>
    !>                          </ul>
    !>                          or,
    !>                          <ul>
    !>                              <li>    a **scalar assumed-length `character`** of kind \SKALL.<br>
    !>                          </ul>
    !>                          If the `lmsize` input argument is missing, then the size of the output `arrayPadded`
    !>                          is `lenArray + lpsize`, otherwise the length of the output `arrayPadded` is `lenArray + lmsize + lpsize`.
    !>  \param[in]  lpsize  :   The input scalar of type `integer` of default kind \IK representing the number of `lpfill` to add to the left of the array.
    !>  \param[in]  lpfill  :   The input scalar of the same type and kind as the input `array` containing the value to fill the left padding (if any)
    !>                          of the output `arrayPadded`. If `array` is of type `character`, then <br>
    !>                          <ol>
    !>                              <li>    the equality `len(lpfill) == 1` must also hold if `array` is a scalar string.<br>
    !>                              <li>    the equality `len(lpfill) == len(array)` must also hold if `array` is an array of strings.<br>
    !>                          </ol>
    !>  \param[in]  lmsize  :   The input scalar `integer` of default kind \IK representing the size of the left-margin of the output `arrayPadded`<br>
    !>                          (**optional**, default = `0`).
    !>  \param[in]  lmfill  :   The input scalar of the same type and kind as the input `array` containing the value to fill the left margin (if any)
    !>                          of newly allocated `array`. If `array` is of type `character`, then <br>
    !>                          <ol>
    !>                              <li> the equality `len(fill) == 1` must also hold if `array` is a scalar string.<br>
    !>                              <li> the equality `len(fill) == len(array)` must also hold if `array` is an array of strings.<br>
    !>                          </ol>
    !>                          (**optional**, if missing, the left-margin will remain uninitialized. It can be present **only if** `lmsize` is also present).
    !>
    !>  \return
    !>  `arrayPadded`      :   The output object of the same type, kind, and rank as the input `array` with the same lower bound,
    !>                          whose size is `lenArray + lpsize + lmsize`, whose contents are the same as the contents of `array`
    !>                          but padded to the left, optionally with the specified left margin.
    !>
    !>  \interface{getPaddedl}
    !>  \code{.F90}
    !>
    !>      use pm_arrayPad, only: getPaddedl
    !>
    !>      arrayPadded = getPaddedl(array, lpsize, lpfill)
    !>      arrayPadded = getPaddedl(array, lpsize, lpfill, lmsize, lmfill = lmfill)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The left margin elements will not be initialized to any particular
    !>  values unless the corresponding `lmfill` argument is present.
    !>
    !>  \warning
    !>  The input `lpsize` and `lmsize` arguments must be non-negative input arguments.<br>
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
    !>  \example{getPaddedl}
    !>  \include{lineno} example/pm_arrayPad/getPaddedl/main.F90
    !>  \compilef{getPaddedl}
    !>  \output{getPaddedl}
    !>  \include{lineno} example/pm_arrayPad/getPaddedl/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayPad](@ref test_pm_arrayPad)
    !>
    !>  \todo
    !>  \pmed Two new optional input scalar `lbcold` and `ubcold` arguments can be added to procedures to specify
    !>  a subset of the contents of the original array that has to be kept in the newly allocated padded array.
    !>
    !>  \final{getPaddedl}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getPaddedl

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getPaddedAsisSL_D0_SK5(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: lpsize
        character(1,SKG)            , intent(in)                    :: lpfill
        character(len(array,IK)+lpsize,SKG)                         :: arrayPadded
    end function
#endif

#if SK4_ENABLED
    PURE module function getPaddedAsisSL_D0_SK4(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: lpsize
        character(1,SKG)            , intent(in)                    :: lpfill
        character(len(array,IK)+lpsize,SKG)                         :: arrayPadded
    end function
#endif

#if SK3_ENABLED
    PURE module function getPaddedAsisSL_D0_SK3(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: lpsize
        character(1,SKG)            , intent(in)                    :: lpfill
        character(len(array,IK)+lpsize,SKG)                         :: arrayPadded
    end function
#endif

#if SK2_ENABLED
    PURE module function getPaddedAsisSL_D0_SK2(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: lpsize
        character(1,SKG)            , intent(in)                    :: lpfill
        character(len(array,IK)+lpsize,SKG)                         :: arrayPadded
    end function
#endif

#if SK1_ENABLED
    PURE module function getPaddedAsisSL_D0_SK1(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: lpsize
        character(1,SKG)            , intent(in)                    :: lpfill
        character(len(array,IK)+lpsize,SKG)                         :: arrayPadded
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getPaddedAsisSL_D1_SK5(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        character(len(array,IK),SKG), intent(in)                    :: lpfill
        character(len(array,IK),SKG)                                :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

#if SK4_ENABLED
    PURE module function getPaddedAsisSL_D1_SK4(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        character(len(array,IK),SKG), intent(in)                    :: lpfill
        character(len(array,IK),SKG)                                :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

#if SK3_ENABLED
    PURE module function getPaddedAsisSL_D1_SK3(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        character(len(array,IK),SKG), intent(in)                    :: lpfill
        character(len(array,IK),SKG)                                :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

#if SK2_ENABLED
    PURE module function getPaddedAsisSL_D1_SK2(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        character(len(array,IK),SKG), intent(in)                    :: lpfill
        character(len(array,IK),SKG)                                :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

#if SK1_ENABLED
    PURE module function getPaddedAsisSL_D1_SK1(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        character(len(array,IK),SKG), intent(in)                    :: lpfill
        character(len(array,IK),SKG)                                :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getPaddedAsisSL_D1_IK5(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        integer(IKG)                , intent(in)                    :: lpfill
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

#if IK4_ENABLED
    PURE module function getPaddedAsisSL_D1_IK4(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        integer(IKG)                , intent(in)                    :: lpfill
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

#if IK3_ENABLED
    PURE module function getPaddedAsisSL_D1_IK3(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        integer(IKG)                , intent(in)                    :: lpfill
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

#if IK2_ENABLED
    PURE module function getPaddedAsisSL_D1_IK2(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        integer(IKG)                , intent(in)                    :: lpfill
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

#if IK1_ENABLED
    PURE module function getPaddedAsisSL_D1_IK1(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        integer(IKG)                , intent(in)                    :: lpfill
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getPaddedAsisSL_D1_LK5(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        logical(LKG)                , intent(in)                    :: lpfill
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

#if LK4_ENABLED
    PURE module function getPaddedAsisSL_D1_LK4(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        logical(LKG)                , intent(in)                    :: lpfill
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

#if LK3_ENABLED
    PURE module function getPaddedAsisSL_D1_LK3(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        logical(LKG)                , intent(in)                    :: lpfill
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

#if LK2_ENABLED
    PURE module function getPaddedAsisSL_D1_LK2(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        logical(LKG)                , intent(in)                    :: lpfill
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

#if LK1_ENABLED
    PURE module function getPaddedAsisSL_D1_LK1(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        logical(LKG)                , intent(in)                    :: lpfill
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getPaddedAsisSL_D1_CK5(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        complex(CKG)                , intent(in)                    :: lpfill
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

#if CK4_ENABLED
    PURE module function getPaddedAsisSL_D1_CK4(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        complex(CKG)                , intent(in)                    :: lpfill
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

#if CK3_ENABLED
    PURE module function getPaddedAsisSL_D1_CK3(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        complex(CKG)                , intent(in)                    :: lpfill
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

#if CK2_ENABLED
    PURE module function getPaddedAsisSL_D1_CK2(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        complex(CKG)                , intent(in)                    :: lpfill
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

#if CK1_ENABLED
    PURE module function getPaddedAsisSL_D1_CK1(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        complex(CKG)                , intent(in)                    :: lpfill
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPaddedAsisSL_D1_RK5(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        real(RKG)                   , intent(in)                    :: lpfill
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

#if RK4_ENABLED
    PURE module function getPaddedAsisSL_D1_RK4(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        real(RKG)                   , intent(in)                    :: lpfill
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

#if RK3_ENABLED
    PURE module function getPaddedAsisSL_D1_RK3(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        real(RKG)                   , intent(in)                    :: lpfill
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

#if RK2_ENABLED
    PURE module function getPaddedAsisSL_D1_RK2(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        real(RKG)                   , intent(in)                    :: lpfill
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

#if RK1_ENABLED
    PURE module function getPaddedAsisSL_D1_RK1(array, lpsize, lpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSL_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        real(RKG)                   , intent(in)                    :: lpfill
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+lpsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getPaddedMargSL_D0_SK5(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)                    :: lpfill
        character(1,SKG)            , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        character(len(array,IK)+lpsize+lmsize,SKG)                  :: arrayPadded
    end function
#endif

#if SK4_ENABLED
    PURE module function getPaddedMargSL_D0_SK4(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)                    :: lpfill
        character(1,SKG)            , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        character(len(array,IK)+lpsize+lmsize,SKG)                  :: arrayPadded
    end function
#endif

#if SK3_ENABLED
    PURE module function getPaddedMargSL_D0_SK3(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)                    :: lpfill
        character(1,SKG)            , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        character(len(array,IK)+lpsize+lmsize,SKG)                  :: arrayPadded
    end function
#endif

#if SK2_ENABLED
    PURE module function getPaddedMargSL_D0_SK2(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)                    :: lpfill
        character(1,SKG)            , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        character(len(array,IK)+lpsize+lmsize,SKG)                  :: arrayPadded
    end function
#endif

#if SK1_ENABLED
    PURE module function getPaddedMargSL_D0_SK1(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)                    :: lpfill
        character(1,SKG)            , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        character(len(array,IK)+lpsize+lmsize,SKG)                  :: arrayPadded
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getPaddedMargSL_D1_SK5(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: lpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        character(len(array,IK),SKG)                                :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

#if SK4_ENABLED
    PURE module function getPaddedMargSL_D1_SK4(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: lpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        character(len(array,IK),SKG)                                :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

#if SK3_ENABLED
    PURE module function getPaddedMargSL_D1_SK3(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: lpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        character(len(array,IK),SKG)                                :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

#if SK2_ENABLED
    PURE module function getPaddedMargSL_D1_SK2(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: lpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        character(len(array,IK),SKG)                                :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

#if SK1_ENABLED
    PURE module function getPaddedMargSL_D1_SK1(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: lpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        character(len(array,IK),SKG)                                :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getPaddedMargSL_D1_IK5(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: lpfill
        integer(IKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

#if IK4_ENABLED
    PURE module function getPaddedMargSL_D1_IK4(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: lpfill
        integer(IKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

#if IK3_ENABLED
    PURE module function getPaddedMargSL_D1_IK3(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: lpfill
        integer(IKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

#if IK2_ENABLED
    PURE module function getPaddedMargSL_D1_IK2(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: lpfill
        integer(IKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

#if IK1_ENABLED
    PURE module function getPaddedMargSL_D1_IK1(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: lpfill
        integer(IKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getPaddedMargSL_D1_LK5(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: lpfill
        logical(LKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

#if LK4_ENABLED
    PURE module function getPaddedMargSL_D1_LK4(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: lpfill
        logical(LKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

#if LK3_ENABLED
    PURE module function getPaddedMargSL_D1_LK3(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: lpfill
        logical(LKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

#if LK2_ENABLED
    PURE module function getPaddedMargSL_D1_LK2(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: lpfill
        logical(LKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

#if LK1_ENABLED
    PURE module function getPaddedMargSL_D1_LK1(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: lpfill
        logical(LKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getPaddedMargSL_D1_CK5(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: lpfill
        complex(CKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

#if CK4_ENABLED
    PURE module function getPaddedMargSL_D1_CK4(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: lpfill
        complex(CKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

#if CK3_ENABLED
    PURE module function getPaddedMargSL_D1_CK3(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: lpfill
        complex(CKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

#if CK2_ENABLED
    PURE module function getPaddedMargSL_D1_CK2(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: lpfill
        complex(CKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

#if CK1_ENABLED
    PURE module function getPaddedMargSL_D1_CK1(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: lpfill
        complex(CKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPaddedMargSL_D1_RK5(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: lpfill
        real(RKG)                   , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

#if RK4_ENABLED
    PURE module function getPaddedMargSL_D1_RK4(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: lpfill
        real(RKG)                   , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

#if RK3_ENABLED
    PURE module function getPaddedMargSL_D1_RK3(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: lpfill
        real(RKG)                   , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

#if RK2_ENABLED
    PURE module function getPaddedMargSL_D1_RK2(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: lpfill
        real(RKG)                   , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

#if RK1_ENABLED
    PURE module function getPaddedMargSL_D1_RK1(array, lpsize, lpfill, lmsize, lmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSL_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: lpfill
        real(RKG)                   , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+lpsize+lmsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Resize the input `array` by padding it to the left with the
    !>  requested paddings and optionally adding margins to the left of the padded `array`
    !>  optionally filled with the corresponding fills.
    !>
    !>  \param[inout]   array   :   The input/output `allocatable` array of shape `(:)` of either <br>
    !>                              <ul>
    !>                                  <li>    type `character` of kind \SKALL, or <br>
    !>                                  <li>    type `logical` of kind \LKALL, or <br>
    !>                                  <li>    type `integer` of kind \IKALL, or <br>
    !>                                  <li>    type `complex` of kind \CKALL, or <br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ul>
    !>                              or,
    !>                              <ul>
    !>                                  <li>    a **scalar allocatable `character`** of default kind \SK.<br>
    !>                              </ul>
    !>                              If the `lmsize` input argument is missing, then the size of the output `array`
    !>                              is `lenArray + lpsize`, otherwise the length of the output `array` is `lenArray + lmsize + lpsize`.<br>
    !>                              On output, the array will be resized (with the same lower bound as before) and padded and marginalized accordingly.
    !>  \param[in]      lpsize  :   The input scalar of type `integer` of default kind \IK representing the number of `lpfill` to add to the left of the array.
    !>  \param[in]      lpfill  :   The input scalar of the same type and kind as the input `array` containing the value to fill the left padding (if any)
    !>                              of the output `array`. If `array` is of type `character`, then <br>
    !>                                  -# the equality `len(lpfill) == 1` must also hold if `array` is a scalar string.<br>
    !>                                  -# the equality `len(lpfill) == len(array)` must also hold if `array` is an array of strings.<br>
    !>  \param[in]      lmsize  :   The input scalar `integer` of default kind \IK representing the size of the left-margin of the output `array`<br>
    !>                              (**optional**, default = `0`).
    !>  \param[in]      lmfill  :   The input scalar of the same type and kind as the input `array` containing the value to fill the left margin (if any)
    !>                              of newly allocated `array`. If `array` is of type `character`, then <br>
    !>                              <ol>
    !>                                  <li> the equality `len(fill) == 1` must also hold if `array` is a scalar string.<br>
    !>                                  <li> the equality `len(fill) == len(array)` must also hold if `array` is an array of strings.<br>
    !>                              </ol>
    !>                              (**optional**, if missing, the left-margin will remain uninitialized. It can be present **only if** `lmsize` is also present).
    !>  \param[in]      failed  :   The input scalar `logical` of default kind \LK that is `.true.` if the requested array resizing and padding is successful<br>
    !>                              (**optional**, if missing and an error occurs, the processor dictates the program behavior).
    !>
    !>  \interface{setPaddedl}
    !>  \code{.F90}
    !>
    !>      use pm_arrayPad, only: setPaddedl
    !>
    !>      call setPaddedl(array, lpsize, lpfill, failed = failed)
    !>      call setPaddedl(array, lpsize, lpfill, lmsize, failed = failed)
    !>      call setPaddedl(array, lpsize, lpfill, lmsize, lmfill = lmfill, failed = failed)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The left margin elements will not be initialized to any particular
    !>  values unless the corresponding `lmfill` argument is present.
    !>
    !>  \warning
    !>  The input `lpsize` and `lmsize` arguments must be non-negative input arguments.<br>
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
    !>  \example{setPaddedl}
    !>  \include{lineno} example/pm_arrayPad/setPaddedl/main.F90
    !>  \compilef{setPaddedl}
    !>  \output{setPaddedl}
    !>  \include{lineno} example/pm_arrayPad/setPaddedl/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayPad](@ref test_pm_arrayPad)
    !>
    !>  \todo
    !>  \pmed Two new optional input scalar `lbcold` and `ubcold` arguments can be added to procedures to specify
    !>  a subset of the contents of the original array that has to be kept in the newly allocated padded array.
    !>
    !>  \final{setPaddedl}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setPaddedl

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setPaddedAsisSL_D0_SK5(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(:,SKG)            , intent(inout) , allocatable   :: array
        integer(IK)                 , intent(in)                    :: lpsize
        character(1,SKG)            , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setPaddedAsisSL_D0_SK4(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)            , intent(inout) , allocatable   :: array
        integer(IK)                 , intent(in)                    :: lpsize
        character(1,SKG)            , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setPaddedAsisSL_D0_SK3(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)            , intent(inout) , allocatable   :: array
        integer(IK)                 , intent(in)                    :: lpsize
        character(1,SKG)            , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setPaddedAsisSL_D0_SK2(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)            , intent(inout) , allocatable   :: array
        integer(IK)                 , intent(in)                    :: lpsize
        character(1,SKG)            , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setPaddedAsisSL_D0_SK1(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)            , intent(inout) , allocatable   :: array
        integer(IK)                 , intent(in)                    :: lpsize
        character(1,SKG)            , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_SK5(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        character(len(array,IK),SKG), intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_SK4(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        character(len(array,IK),SKG), intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_SK3(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        character(len(array,IK),SKG), intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_SK2(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        character(len(array,IK),SKG), intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_SK1(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        character(len(array,IK),SKG), intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_IK5(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        integer(IKG)                , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_IK4(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        integer(IKG)                , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_IK3(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        integer(IKG)                , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_IK2(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        integer(IKG)                , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_IK1(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        integer(IKG)                , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_LK5(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        logical(LKG)                , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_LK4(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        logical(LKG)                , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_LK3(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        logical(LKG)                , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_LK2(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        logical(LKG)                , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_LK1(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        logical(LKG)                , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_CK5(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        complex(CKG)                , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_CK4(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        complex(CKG)                , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_CK3(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        complex(CKG)                , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_CK2(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        complex(CKG)                , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_CK1(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        complex(CKG)                , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_RK5(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        real(RKG)                   , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_RK4(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        real(RKG)                   , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_RK3(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        real(RKG)                   , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_RK2(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        real(RKG)                   , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPaddedAsisSL_D1_RK1(array, lpsize, lpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSL_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: lpsize
        real(RKG)                   , intent(in)                    :: lpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setPaddedMargSL_D0_SK5(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(:,SKG)            , intent(inout) , allocatable   :: array
        character(1,SKG)            , intent(in)                    :: lpfill
        character(1,SKG)            , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setPaddedMargSL_D0_SK4(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)            , intent(inout) , allocatable   :: array
        character(1,SKG)            , intent(in)                    :: lpfill
        character(1,SKG)            , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setPaddedMargSL_D0_SK3(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)            , intent(inout) , allocatable   :: array
        character(1,SKG)            , intent(in)                    :: lpfill
        character(1,SKG)            , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setPaddedMargSL_D0_SK2(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)            , intent(inout) , allocatable   :: array
        character(1,SKG)            , intent(in)                    :: lpfill
        character(1,SKG)            , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setPaddedMargSL_D0_SK1(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)            , intent(inout) , allocatable   :: array
        character(1,SKG)            , intent(in)                    :: lpfill
        character(1,SKG)            , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setPaddedMargSL_D1_SK5(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: lpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setPaddedMargSL_D1_SK4(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: lpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setPaddedMargSL_D1_SK3(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: lpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setPaddedMargSL_D1_SK2(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: lpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setPaddedMargSL_D1_SK1(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: lpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setPaddedMargSL_D1_IK5(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IKG)                , intent(in)                    :: lpfill
        integer(IKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setPaddedMargSL_D1_IK4(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IKG)                , intent(in)                    :: lpfill
        integer(IKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setPaddedMargSL_D1_IK3(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IKG)                , intent(in)                    :: lpfill
        integer(IKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setPaddedMargSL_D1_IK2(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IKG)                , intent(in)                    :: lpfill
        integer(IKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setPaddedMargSL_D1_IK1(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IKG)                , intent(in)                    :: lpfill
        integer(IKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setPaddedMargSL_D1_LK5(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        logical(LKG)                , intent(in)                    :: lpfill
        logical(LKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setPaddedMargSL_D1_LK4(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        logical(LKG)                , intent(in)                    :: lpfill
        logical(LKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setPaddedMargSL_D1_LK3(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        logical(LKG)                , intent(in)                    :: lpfill
        logical(LKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setPaddedMargSL_D1_LK2(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        logical(LKG)                , intent(in)                    :: lpfill
        logical(LKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setPaddedMargSL_D1_LK1(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        logical(LKG)                , intent(in)                    :: lpfill
        logical(LKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPaddedMargSL_D1_CK5(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        complex(CKG)                , intent(in)                    :: lpfill
        complex(CKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPaddedMargSL_D1_CK4(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        complex(CKG)                , intent(in)                    :: lpfill
        complex(CKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPaddedMargSL_D1_CK3(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        complex(CKG)                , intent(in)                    :: lpfill
        complex(CKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPaddedMargSL_D1_CK2(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        complex(CKG)                , intent(in)                    :: lpfill
        complex(CKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPaddedMargSL_D1_CK1(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        complex(CKG)                , intent(in)                    :: lpfill
        complex(CKG)                , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPaddedMargSL_D1_RK5(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        real(RKG)                   , intent(in)                    :: lpfill
        real(RKG)                   , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPaddedMargSL_D1_RK4(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        real(RKG)                   , intent(in)                    :: lpfill
        real(RKG)                   , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPaddedMargSL_D1_RK3(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        real(RKG)                   , intent(in)                    :: lpfill
        real(RKG)                   , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPaddedMargSL_D1_RK2(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        real(RKG)                   , intent(in)                    :: lpfill
        real(RKG)                   , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPaddedMargSL_D1_RK1(array, lpsize, lpfill, lmsize, lmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSL_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        real(RKG)                   , intent(in)                    :: lpfill
        real(RKG)                   , intent(in)    , optional      :: lmfill
        integer(IK)                 , intent(in)                    :: lpsize, lmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate a resized copy of the input `array` by padding it to the right with the
    !>  requested paddings and optionally adding margins to the right of the padded `array`
    !>  optionally filled with the corresponding fills.
    !>
    !>  \param[in]  array   :   The input `contiguous` array of shape `(:)` of either <br>
    !>                          <ul>
    !>                              <li>    type `character` of kind \SKALL, or <br>
    !>                              <li>    type `integer` of kind \IKALL, or <br>
    !>                              <li>    type `logical` of kind \LKALL, or <br>
    !>                              <li>    type `complex` of kind \CKALL, or <br>
    !>                              <li>    type `real` of kind \RKALL, or <br>
    !>                          </ul>
    !>                          or,
    !>                          <ul>
    !>                              <li>    a **scalar assumed-length `character`** of default kind \SK.<br>
    !>                          </ul>
    !>                          If the `rmsize` input argument is missing, then the size of the output `arrayPadded`
    !>                          is `lenArray + rpsize`, otherwise the length of the output `arrayPadded` is `lenArray + rmsize + rpsize`.
    !>  \param[in]  rpsize  :   The input scalar of type `integer` of default kind \IK representing the number of `rpfill` to add to the right of the array.
    !>  \param[in]  rpfill  :   The input scalar of the same type and kind as the input `array` containing the value to fill the right padding (if any)
    !>                          of the output `arrayPadded`. If `array` is of type `character`, then <br>
    !>                          <ol>
    !>                              <li>    the equality `len(rpfill) == 1` must also hold if `array` is a scalar string.<br>
    !>                              <li>    the equality `len(rpfill) == len(array)` must also hold if `array` is an array of strings.<br>
    !>                          </ol>
    !>  \param[in]  rmsize  :   The input scalar `integer` of default kind \IK representing the size of the right-margin of the output `arrayPadded`<br>
    !>                          (**optional**, default = `0`).
    !>  \param[in]  rmfill  :   The input scalar of the same type and kind as the input `array` containing the value to fill the right margin (if any)
    !>                          of newly allocated `array`. If `array` is of type `character`, then <br>
    !>                          <ol>
    !>                              <li> the equality `len(fill) == 1` must also hold if `array` is a scalar string.<br>
    !>                              <li> the equality `len(fill) == len(array)` must also hold if `array` is an array of strings.<br>
    !>                          </ol>
    !>                          (**optional**, if missing, the right-margin will remain uninitialized. It can be present **only if** `rmsize` is also present).
    !>
    !>  \return
    !>  `arrayPadded`      :   The output object of the same type, kind, and rank as the input `array` with the same lower bound,
    !>                          whose size is `lenArray + rpsize + rmsize`, whose contents are the same as the contents of `array`
    !>                          but padded to the right, optionally with the specified right margin.
    !>
    !>  \interface{getPaddedr}
    !>  \code{.F90}
    !>
    !>      use pm_arrayPad, only: getPaddedr
    !>
    !>      arrayPadded = getPaddedr(array, rpsize, rpfill)
    !>      arrayPadded = getPaddedr(array, rpsize, rpfill, rmsize, rmfill = rmfill)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The right margin elements will not be initialized to any particular values unless the corresponding `rmfill` argument is present.<br>
    !>
    !>  \warning
    !>  The input `rpsize` and `rmsize` arguments must be non-negative input arguments.<br>
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
    !>  \example{getPaddedr}
    !>  \include{lineno} example/pm_arrayPad/getPaddedr/main.F90
    !>  \compilef{getPaddedr}
    !>  \output{getPaddedr}
    !>  \include{lineno} example/pm_arrayPad/getPaddedr/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayPad](@ref test_pm_arrayPad)
    !>
    !>  \todo
    !>  \pmed Two new optional input scalar `lbcold` and `ubcold` arguments can be added to procedures to specify
    !>  a subset of the contents of the original array that has to be kept in the newly allocated padded array.
    !>
    !>  \final{getPaddedr}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getPaddedr

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getPaddedAsisSR_D0_SK5(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: rpsize
        character(1,SKG)            , intent(in)                    :: rpfill
        character(len(array,IK)+rpsize,SKG)                         :: arrayPadded
    end function
#endif

#if SK4_ENABLED
    PURE module function getPaddedAsisSR_D0_SK4(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: rpsize
        character(1,SKG)            , intent(in)                    :: rpfill
        character(len(array,IK)+rpsize,SKG)                         :: arrayPadded
    end function
#endif

#if SK3_ENABLED
    PURE module function getPaddedAsisSR_D0_SK3(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: rpsize
        character(1,SKG)            , intent(in)                    :: rpfill
        character(len(array,IK)+rpsize,SKG)                         :: arrayPadded
    end function
#endif

#if SK2_ENABLED
    PURE module function getPaddedAsisSR_D0_SK2(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: rpsize
        character(1,SKG)            , intent(in)                    :: rpfill
        character(len(array,IK)+rpsize,SKG)                         :: arrayPadded
    end function
#endif

#if SK1_ENABLED
    PURE module function getPaddedAsisSR_D0_SK1(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: rpsize
        character(1,SKG)            , intent(in)                    :: rpfill
        character(len(array,IK)+rpsize,SKG)                         :: arrayPadded
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getPaddedAsisSR_D1_SK5(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        character(len(array,IK),SKG), intent(in)                    :: rpfill
        character(len(array,IK),SKG)                                :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

#if SK4_ENABLED
    PURE module function getPaddedAsisSR_D1_SK4(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        character(len(array,IK),SKG), intent(in)                    :: rpfill
        character(len(array,IK),SKG)                                :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

#if SK3_ENABLED
    PURE module function getPaddedAsisSR_D1_SK3(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        character(len(array,IK),SKG), intent(in)                    :: rpfill
        character(len(array,IK),SKG)                                :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

#if SK2_ENABLED
    PURE module function getPaddedAsisSR_D1_SK2(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        character(len(array,IK),SKG), intent(in)                    :: rpfill
        character(len(array,IK),SKG)                                :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

#if SK1_ENABLED
    PURE module function getPaddedAsisSR_D1_SK1(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        character(len(array,IK),SKG), intent(in)                    :: rpfill
        character(len(array,IK),SKG)                                :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getPaddedAsisSR_D1_IK5(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        integer(IKG)                , intent(in)                    :: rpfill
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

#if IK4_ENABLED
    PURE module function getPaddedAsisSR_D1_IK4(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        integer(IKG)                , intent(in)                    :: rpfill
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

#if IK3_ENABLED
    PURE module function getPaddedAsisSR_D1_IK3(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        integer(IKG)                , intent(in)                    :: rpfill
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

#if IK2_ENABLED
    PURE module function getPaddedAsisSR_D1_IK2(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        integer(IKG)                , intent(in)                    :: rpfill
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

#if IK1_ENABLED
    PURE module function getPaddedAsisSR_D1_IK1(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        integer(IKG)                , intent(in)                    :: rpfill
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getPaddedAsisSR_D1_LK5(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        logical(LKG)                , intent(in)                    :: rpfill
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

#if LK4_ENABLED
    PURE module function getPaddedAsisSR_D1_LK4(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        logical(LKG)                , intent(in)                    :: rpfill
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

#if LK3_ENABLED
    PURE module function getPaddedAsisSR_D1_LK3(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        logical(LKG)                , intent(in)                    :: rpfill
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

#if LK2_ENABLED
    PURE module function getPaddedAsisSR_D1_LK2(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        logical(LKG)                , intent(in)                    :: rpfill
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

#if LK1_ENABLED
    PURE module function getPaddedAsisSR_D1_LK1(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        logical(LKG)                , intent(in)                    :: rpfill
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getPaddedAsisSR_D1_CK5(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        complex(CKG)                , intent(in)                    :: rpfill
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

#if CK4_ENABLED
    PURE module function getPaddedAsisSR_D1_CK4(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        complex(CKG)                , intent(in)                    :: rpfill
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

#if CK3_ENABLED
    PURE module function getPaddedAsisSR_D1_CK3(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        complex(CKG)                , intent(in)                    :: rpfill
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

#if CK2_ENABLED
    PURE module function getPaddedAsisSR_D1_CK2(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        complex(CKG)                , intent(in)                    :: rpfill
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

#if CK1_ENABLED
    PURE module function getPaddedAsisSR_D1_CK1(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        complex(CKG)                , intent(in)                    :: rpfill
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPaddedAsisSR_D1_RK5(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        real(RKG)                   , intent(in)                    :: rpfill
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

#if RK4_ENABLED
    PURE module function getPaddedAsisSR_D1_RK4(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        real(RKG)                   , intent(in)                    :: rpfill
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

#if RK3_ENABLED
    PURE module function getPaddedAsisSR_D1_RK3(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        real(RKG)                   , intent(in)                    :: rpfill
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

#if RK2_ENABLED
    PURE module function getPaddedAsisSR_D1_RK2(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        real(RKG)                   , intent(in)                    :: rpfill
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

#if RK1_ENABLED
    PURE module function getPaddedAsisSR_D1_RK1(array, rpsize, rpfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedAsisSR_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        real(RKG)                   , intent(in)                    :: rpfill
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+rpsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getPaddedMargSR_D0_SK5(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)                    :: rpfill
        character(1,SKG)            , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        character(len(array,IK)+rpsize+rmsize,SKG)                  :: arrayPadded
    end function
#endif

#if SK4_ENABLED
    PURE module function getPaddedMargSR_D0_SK4(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)                    :: rpfill
        character(1,SKG)            , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        character(len(array,IK)+rpsize+rmsize,SKG)                  :: arrayPadded
    end function
#endif

#if SK3_ENABLED
    PURE module function getPaddedMargSR_D0_SK3(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)                    :: rpfill
        character(1,SKG)            , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        character(len(array,IK)+rpsize+rmsize,SKG)                  :: arrayPadded
    end function
#endif

#if SK2_ENABLED
    PURE module function getPaddedMargSR_D0_SK2(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)                    :: rpfill
        character(1,SKG)            , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        character(len(array,IK)+rpsize+rmsize,SKG)                  :: arrayPadded
    end function
#endif

#if SK1_ENABLED
    PURE module function getPaddedMargSR_D0_SK1(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        character(1,SKG)            , intent(in)                    :: rpfill
        character(1,SKG)            , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        character(len(array,IK)+rpsize+rmsize,SKG)                  :: arrayPadded
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getPaddedMargSR_D1_SK5(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: rpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        character(len(array,IK),SKG)                                :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

#if SK4_ENABLED
    PURE module function getPaddedMargSR_D1_SK4(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: rpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        character(len(array,IK),SKG)                                :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

#if SK3_ENABLED
    PURE module function getPaddedMargSR_D1_SK3(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: rpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        character(len(array,IK),SKG)                                :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

#if SK2_ENABLED
    PURE module function getPaddedMargSR_D1_SK2(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: rpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        character(len(array,IK),SKG)                                :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

#if SK1_ENABLED
    PURE module function getPaddedMargSR_D1_SK1(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: rpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        character(len(array,IK),SKG)                                :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getPaddedMargSR_D1_IK5(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: rpfill
        integer(IKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

#if IK4_ENABLED
    PURE module function getPaddedMargSR_D1_IK4(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: rpfill
        integer(IKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

#if IK3_ENABLED
    PURE module function getPaddedMargSR_D1_IK3(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: rpfill
        integer(IKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

#if IK2_ENABLED
    PURE module function getPaddedMargSR_D1_IK2(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: rpfill
        integer(IKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

#if IK1_ENABLED
    PURE module function getPaddedMargSR_D1_IK1(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: rpfill
        integer(IKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        integer(IKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getPaddedMargSR_D1_LK5(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: rpfill
        logical(LKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

#if LK4_ENABLED
    PURE module function getPaddedMargSR_D1_LK4(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: rpfill
        logical(LKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

#if LK3_ENABLED
    PURE module function getPaddedMargSR_D1_LK3(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: rpfill
        logical(LKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

#if LK2_ENABLED
    PURE module function getPaddedMargSR_D1_LK2(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: rpfill
        logical(LKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

#if LK1_ENABLED
    PURE module function getPaddedMargSR_D1_LK1(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: rpfill
        logical(LKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getPaddedMargSR_D1_CK5(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: rpfill
        complex(CKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

#if CK4_ENABLED
    PURE module function getPaddedMargSR_D1_CK4(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: rpfill
        complex(CKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

#if CK3_ENABLED
    PURE module function getPaddedMargSR_D1_CK3(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: rpfill
        complex(CKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

#if CK2_ENABLED
    PURE module function getPaddedMargSR_D1_CK2(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: rpfill
        complex(CKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

#if CK1_ENABLED
    PURE module function getPaddedMargSR_D1_CK1(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: rpfill
        complex(CKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        complex(CKG)                                                :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPaddedMargSR_D1_RK5(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: rpfill
        real(RKG)                   , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

#if RK4_ENABLED
    PURE module function getPaddedMargSR_D1_RK4(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: rpfill
        real(RKG)                   , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

#if RK3_ENABLED
    PURE module function getPaddedMargSR_D1_RK3(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: rpfill
        real(RKG)                   , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

#if RK2_ENABLED
    PURE module function getPaddedMargSR_D1_RK2(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: rpfill
        real(RKG)                   , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

#if RK1_ENABLED
    PURE module function getPaddedMargSR_D1_RK1(array, rpsize, rpfill, rmsize, rmfill) result(arrayPadded)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedMargSR_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: rpfill
        real(RKG)                   , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        real(RKG)                                                   :: arrayPadded(size(array,kind=IK)+rpsize+rmsize)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Resize the input `array` by padding it to the right with the
    !>  requested paddings and optionally adding margins to the right of the padded `array`
    !>  optionally filled with the corresponding fills.
    !>
    !>  \param[inout]   array   :   The input/output `allocatable` array of shape `(:)` of either <br>
    !>                              <ul>
    !>                                  <li>    type `character` of kind \SKALL, or <br>
    !>                                  <li>    type `integer` of kind \IKALL, or <br>
    !>                                  <li>    type `logical` of kind \LKALL, or <br>
    !>                                  <li>    type `complex` of kind \CKALL, or <br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ul>
    !>                              or,
    !>                              <ul>
    !>                                  <li>    a **scalar allocatable `character`** of default kind \SK.<br>
    !>                              </ul>
    !>                              If the `rmsize` input argument is missing, then the size of the output `array`
    !>                              is `lenArray + rpsize`, otherwise the length of the output `array` is `lenArray + rmsize + rpsize`.<br>
    !>                              On output, the array will be resized (with the same lower bound as before) and padded and marginalized accordingly.
    !>  \param[in]      rpsize  :   The input scalar of type `integer` of default kind \IK representing the number of `rpfill` to add to the right of the array.
    !>  \param[in]      rpfill  :   The input scalar of the same type and kind as the input `array` containing the value to fill the right padding (if any)
    !>                              of the output `array`. If `array` is of type `character`, then <br>
    !>                                  -# the equality `len(rpfill) == 1` must also hold if `array` is a scalar string.<br>
    !>                                  -# the equality `len(rpfill) == len(array)` must also hold if `array` is an array of strings.<br>
    !>  \param[in]      rmsize  :   The input scalar `integer` of default kind \IK representing the size of the right-margin of the output `array`<br>
    !>                              (**optional**, default = `0`).
    !>  \param[in]      rmfill  :   The input scalar of the same type and kind as the input `array` containing the value to fill the right margin (if any)
    !>                              of newly allocated `array`. If `array` is of type `character`, then <br>
    !>                              <ol>
    !>                                  <li> the equality `len(fill) == 1` must also hold if `array` is a scalar string.<br>
    !>                                  <li> the equality `len(fill) == len(array)` must also hold if `array` is an array of strings.<br>
    !>                              </ol>
    !>                              (**optional**, if missing, the right-margin will remain uninitialized. It can be present **only if** `rmsize` is also present).
    !>  \param[in]      failed  :   The input scalar `logical` of default kind \LK that is `.true.` if the requested array resizing and padding is successful<br>
    !>                              (**optional**, if missing and an error occurs, the processor dictates the program behavior).
    !>
    !>  \interface{setPaddedr}
    !>  \code{.F90}
    !>
    !>      use pm_arrayPad, only: setPaddedr
    !>
    !>      call setPaddedr(array, rpsize, rpfill)
    !>      call setPaddedr(array, rpsize, rpfill, rmsize, rmfill = rmfill, failed = failed)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The right margin elements will not be initialized to any particular
    !>  values unless the corresponding `rmfill` argument is present.
    !>
    !>  \warning
    !>  The input `rpsize` and `rmsize` arguments must be non-negative input arguments.<br>
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
    !>  \example{setPaddedr}
    !>  \include{lineno} example/pm_arrayPad/setPaddedr/main.F90
    !>  \compilef{setPaddedr}
    !>  \output{setPaddedr}
    !>  \include{lineno} example/pm_arrayPad/setPaddedr/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayPad](@ref test_pm_arrayPad)
    !>
    !>  \todo
    !>  \pmed Two new optional input scalar `lbcold` and `ubcold` arguments can be added to procedures to specify
    !>  a subset of the contents of the original array that has to be kept in the newly allocated padded array.
    !>
    !>  \final{setPaddedr}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setPaddedr

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setPaddedAsisSR_D0_SK5(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(:,SKG)            , intent(inout) , allocatable   :: array
        integer(IK)                 , intent(in)                    :: rpsize
        character(1,SKG)            , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setPaddedAsisSR_D0_SK4(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)            , intent(inout) , allocatable   :: array
        integer(IK)                 , intent(in)                    :: rpsize
        character(1,SKG)            , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setPaddedAsisSR_D0_SK3(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)            , intent(inout) , allocatable   :: array
        integer(IK)                 , intent(in)                    :: rpsize
        character(1,SKG)            , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setPaddedAsisSR_D0_SK2(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)            , intent(inout) , allocatable   :: array
        integer(IK)                 , intent(in)                    :: rpsize
        character(1,SKG)            , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setPaddedAsisSR_D0_SK1(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)            , intent(inout) , allocatable   :: array
        integer(IK)                 , intent(in)                    :: rpsize
        character(1,SKG)            , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_SK5(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        character(len(array,IK),SKG), intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_SK4(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        character(len(array,IK),SKG), intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_SK3(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        character(len(array,IK),SKG), intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_SK2(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        character(len(array,IK),SKG), intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_SK1(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        character(len(array,IK),SKG), intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_IK5(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        integer(IKG)                , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_IK4(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        integer(IKG)                , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_IK3(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        integer(IKG)                , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_IK2(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        integer(IKG)                , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_IK1(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        integer(IKG)                , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_LK5(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        logical(LKG)                , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_LK4(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        logical(LKG)                , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_LK3(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        logical(LKG)                , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_LK2(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        logical(LKG)                , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_LK1(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        logical(LKG)                , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_CK5(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        complex(CKG)                , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_CK4(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        complex(CKG)                , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_CK3(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        complex(CKG)                , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_CK2(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        complex(CKG)                , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_CK1(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        complex(CKG)                , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_RK5(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        real(RKG)                   , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_RK4(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        real(RKG)                   , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_RK3(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        real(RKG)                   , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_RK2(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        real(RKG)                   , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPaddedAsisSR_D1_RK1(array, rpsize, rpfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedAsisSR_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        integer(IK)                 , intent(in)                    :: rpsize
        real(RKG)                   , intent(in)                    :: rpfill
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setPaddedMargSR_D0_SK5(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , allocatable   :: array
        character(1,SKG)            , intent(in)                    :: rpfill
        character(1,SKG)            , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setPaddedMargSR_D0_SK4(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)            , intent(inout) , allocatable   :: array
        character(1,SKG)            , intent(in)                    :: rpfill
        character(1,SKG)            , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setPaddedMargSR_D0_SK3(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)            , intent(inout) , allocatable   :: array
        character(1,SKG)            , intent(in)                    :: rpfill
        character(1,SKG)            , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setPaddedMargSR_D0_SK2(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)            , intent(inout) , allocatable   :: array
        character(1,SKG)            , intent(in)                    :: rpfill
        character(1,SKG)            , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setPaddedMargSR_D0_SK1(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)            , intent(inout) , allocatable   :: array
        character(1,SKG)            , intent(in)                    :: rpfill
        character(1,SKG)            , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setPaddedMargSR_D1_SK5(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: rpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setPaddedMargSR_D1_SK4(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: rpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setPaddedMargSR_D1_SK3(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: rpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setPaddedMargSR_D1_SK2(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: rpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setPaddedMargSR_D1_SK1(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , allocatable   :: array(:)
        character(len(array,IK),SKG), intent(in)                    :: rpfill
        character(len(array,IK),SKG), intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setPaddedMargSR_D1_IK5(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IKG)                , intent(in)                    :: rpfill
        integer(IKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setPaddedMargSR_D1_IK4(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IKG)                , intent(in)                    :: rpfill
        integer(IKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setPaddedMargSR_D1_IK3(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IKG)                , intent(in)                    :: rpfill
        integer(IKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setPaddedMargSR_D1_IK2(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IKG)                , intent(in)                    :: rpfill
        integer(IKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setPaddedMargSR_D1_IK1(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , allocatable   :: array(:)
        integer(IKG)                , intent(in)                    :: rpfill
        integer(IKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setPaddedMargSR_D1_LK5(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        logical(LKG)                , intent(in)                    :: rpfill
        logical(LKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setPaddedMargSR_D1_LK4(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        logical(LKG)                , intent(in)                    :: rpfill
        logical(LKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setPaddedMargSR_D1_LK3(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        logical(LKG)                , intent(in)                    :: rpfill
        logical(LKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setPaddedMargSR_D1_LK2(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        logical(LKG)                , intent(in)                    :: rpfill
        logical(LKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setPaddedMargSR_D1_LK1(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , allocatable   :: array(:)
        logical(LKG)                , intent(in)                    :: rpfill
        logical(LKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPaddedMargSR_D1_CK5(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        complex(CKG)                , intent(in)                    :: rpfill
        complex(CKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPaddedMargSR_D1_CK4(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        complex(CKG)                , intent(in)                    :: rpfill
        complex(CKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPaddedMargSR_D1_CK3(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        complex(CKG)                , intent(in)                    :: rpfill
        complex(CKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPaddedMargSR_D1_CK2(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        complex(CKG)                , intent(in)                    :: rpfill
        complex(CKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPaddedMargSR_D1_CK1(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , allocatable   :: array(:)
        complex(CKG)                , intent(in)                    :: rpfill
        complex(CKG)                , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPaddedMargSR_D1_RK5(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        real(RKG)                   , intent(in)                    :: rpfill
        real(RKG)                   , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPaddedMargSR_D1_RK4(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        real(RKG)                   , intent(in)                    :: rpfill
        real(RKG)                   , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPaddedMargSR_D1_RK3(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        real(RKG)                   , intent(in)                    :: rpfill
        real(RKG)                   , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPaddedMargSR_D1_RK2(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        real(RKG)                   , intent(in)                    :: rpfill
        real(RKG)                   , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPaddedMargSR_D1_RK1(array, rpsize, rpfill, rmsize, rmfill, failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPaddedMargSR_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , allocatable   :: array(:)
        real(RKG)                   , intent(in)                    :: rpfill
        real(RKG)                   , intent(in)    , optional      :: rmfill
        integer(IK)                 , intent(in)                    :: rpsize, rmsize
        logical(LK)                 , intent(out)   , optional      :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayPad ! LCOV_EXCL_LINE