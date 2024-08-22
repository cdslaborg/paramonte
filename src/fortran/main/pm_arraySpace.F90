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
!>  This module contains procedures and generic interfaces for generating arrays with linear or logarithmic spacing.
!>
!>  \see
!>  [pm_arrayRange](@ref pm_arrayRange)<br>
!>
!>  \test
!>  [test_pm_arraySpace](@ref test_pm_arraySpace)
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Tuesday 08:49 PM, August 10, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arraySpace

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_arraySpace"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate `count` evenly spaced points over the interval `[x1, x2]` if `x1 < x2`, or `[x2, x1]` if `x2 < x1`.
    !>
    !>  \param[in]  x1      :   The input scalar of either<br>
    !>                          <ul>
    !>                              <li>    type `complex` of kind \CKALL or
    !>                              <li>    type `real` of kind \RKALL
    !>                          </ul>
    !>                          representing the starting value of the sequence.
    !>  \param[in]  x2      :   The input scalar of the same type and kind as `x1` representing the ending value of the sequence.
    !>  \param[in]  count   :   The input scalar of type `integer` of default kind \IK representing the length of `linSpace` to generate.
    !>  \param[in]  fopen   :   The input scalar of type `logical` of default kind \LK.
    !>                          If `.true.`, the `linSpace` will be <b>f</b>irst-open, meaning that `x1` will NOT be in the output `linSpace` sequence.<br>
    !>                          (**optional**, default = `.false.`)
    !>  \param[in]  lopen   :   The input scalar of type `logical` of default kind \LK.
    !>                          If `.true.`, the `linSpace` will be <b>l</b>ast-open, meaning that `x2` will NOT be in the output `linSpace` sequence.<br>
    !>                          (**optional**, default = `.false.`)
    !>
    !>  \return
    !>  `linSpace`          :   The output array of shape `(1:count)` of the same type and kind as the input `x1`
    !>                          containing the evenly-spaced sequence within the interval specified by `x1` and `x2`.
    !>
    !>  \interface{getLinSpace}
    !>  \code{.F90}
    !>
    !>      use pm_arraySpace, only: getLinSpace
    !>
    !>      linSpace = getLinSpace(x1, x2, count)
    !>      linSpace = getLinSpace(x1, x2, count, fopen)
    !>      linSpace = getLinSpace(x1, x2, count, fopen, lopen)
    !>      linSpace = getLinSpace(x1, x2, count, lopen = lopen)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < count` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  Setting both `fopen = .true.` and `lopen = .true.` will lead to an output `linSpace` whose points are centers of the `count + 1`
    !>  equally-sized bins in the interval `[x1,x2]`.
    !>
    !>  \remark
    !>  If `count == 1` holds while `fopen = .false.` and `lopen = .false.`, then `linSpace = [x1]` on return.
    !>
    !>  \see
    !>  [getLinSpace](@ref pm_arraySpace::getLinSpace)<br>
    !>  [getLogSpace](@ref pm_arraySpace::getLogSpace)<br>
    !>  [setLogSpace](@ref pm_arraySpace::setLogSpace)<br>
    !>
    !>  \example{getLinSpace}
    !>  \include{lineno} example/pm_arraySpace/getLinSpace/main.F90
    !>  \compilef{getLinSpace}
    !>  \output{getLinSpace}
    !>  \include{lineno} example/pm_arraySpace/getLinSpace/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arraySpace](@ref test_pm_arraySpace)
    !>
    !>  \final{getLinSpace}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 08:49 PM, August 10, 2021, Dallas, TX
    interface getLinSpace

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getLinSpace_CK5(x1, x2, count, fopen, lopen) result(linSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLinSpace_CK5
#endif
        use pm_kind, only: CKG => CK5
        implicit none
        complex(CKG), intent(in)            :: x1
        complex(CKG), intent(in)            :: x2
        integer(IK) , intent(in)            :: count
        logical(LK) , intent(in), optional  :: fopen, lopen
        complex(CKG)                        :: linSpace(0 : count - 1)
    end function
#endif

#if CK4_ENABLED
    PURE module function getLinSpace_CK4(x1, x2, count, fopen, lopen) result(linSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLinSpace_CK4
#endif
        use pm_kind, only: CKG => CK4
        implicit none
        complex(CKG), intent(in)            :: x1
        complex(CKG), intent(in)            :: x2
        integer(IK) , intent(in)            :: count
        logical(LK) , intent(in), optional  :: fopen, lopen
        complex(CKG)                        :: linSpace(0 : count - 1)
    end function
#endif

#if CK3_ENABLED
    PURE module function getLinSpace_CK3(x1, x2, count, fopen, lopen) result(linSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLinSpace_CK3
#endif
        use pm_kind, only: CKG => CK3
        implicit none
        complex(CKG), intent(in)            :: x1
        complex(CKG), intent(in)            :: x2
        integer(IK) , intent(in)            :: count
        logical(LK) , intent(in), optional  :: fopen, lopen
        complex(CKG)                        :: linSpace(0 : count - 1)
    end function
#endif

#if CK2_ENABLED
    PURE module function getLinSpace_CK2(x1, x2, count, fopen, lopen) result(linSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLinSpace_CK2
#endif
        use pm_kind, only: CKG => CK2
        implicit none
        complex(CKG), intent(in)            :: x1
        complex(CKG), intent(in)            :: x2
        integer(IK) , intent(in)            :: count
        logical(LK) , intent(in), optional  :: fopen, lopen
        complex(CKG)                        :: linSpace(0 : count - 1)
    end function
#endif

#if CK1_ENABLED
    PURE module function getLinSpace_CK1(x1, x2, count, fopen, lopen) result(linSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLinSpace_CK1
#endif
        use pm_kind, only: CKG => CK1
        implicit none
        complex(CKG), intent(in)            :: x1
        complex(CKG), intent(in)            :: x2
        integer(IK) , intent(in)            :: count
        logical(LK) , intent(in), optional  :: fopen, lopen
        complex(CKG)                        :: linSpace(0 : count - 1)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getLinSpace_RK5(x1, x2, count, fopen, lopen) result(linSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLinSpace_RK5
#endif
        use pm_kind, only: RKG => RK5
        implicit none
        real(RKG)   , intent(in)            :: x1
        real(RKG)   , intent(in)            :: x2
        integer(IK) , intent(in)            :: count
        logical(LK) , intent(in), optional  :: fopen, lopen
        real(RKG)                           :: linSpace(0 : count - 1)
    end function
#endif

#if RK4_ENABLED
    PURE module function getLinSpace_RK4(x1, x2, count, fopen, lopen) result(linSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLinSpace_RK4
#endif
        use pm_kind, only: RKG => RK4
        implicit none
        real(RKG)   , intent(in)            :: x1
        real(RKG)   , intent(in)            :: x2
        integer(IK) , intent(in)            :: count
        logical(LK) , intent(in), optional  :: fopen, lopen
        real(RKG)                           :: linSpace(0 : count - 1)
    end function
#endif

#if RK3_ENABLED
    PURE module function getLinSpace_RK3(x1, x2, count, fopen, lopen) result(linSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLinSpace_RK3
#endif
        use pm_kind, only: RKG => RK3
        implicit none
        real(RKG)   , intent(in)            :: x1
        real(RKG)   , intent(in)            :: x2
        integer(IK) , intent(in)            :: count
        logical(LK) , intent(in), optional  :: fopen, lopen
        real(RKG)                           :: linSpace(0 : count - 1)
    end function
#endif

#if RK2_ENABLED
    PURE module function getLinSpace_RK2(x1, x2, count, fopen, lopen) result(linSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLinSpace_RK2
#endif
        use pm_kind, only: RKG => RK2
        implicit none
        real(RKG)   , intent(in)            :: x1
        real(RKG)   , intent(in)            :: x2
        integer(IK) , intent(in)            :: count
        logical(LK) , intent(in), optional  :: fopen, lopen
        real(RKG)                           :: linSpace(0 : count - 1)
    end function
#endif

#if RK1_ENABLED
    PURE module function getLinSpace_RK1(x1, x2, count, fopen, lopen) result(linSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLinSpace_RK1
#endif
        use pm_kind, only: RKG => RK1
        implicit none
        real(RKG)   , intent(in)            :: x1
        real(RKG)   , intent(in)            :: x2
        integer(IK) , intent(in)            :: count
        logical(LK) , intent(in), optional  :: fopen, lopen
        real(RKG)                           :: linSpace(0 : count - 1)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the `linSpace` output argument with `size(linSpace)` elements of evenly-spaced values over the interval `[x1, x2]` if `x1 < x2`, or `[x2, x1]` if `x2 < x1`.
    !>
    !>  \param[out] linSpace    :   The output `contiguous` array of either<br>
    !>                              <ol>
    !>                                  <li> type `real` of kind \RKALL or
    !>                                  <li> type `complex` of kind \CKALL
    !>                              </ol>
    !>                              containing the evenly spaced numbers.
    !>  \param[in]  x1          :   The input scalar of the same type and kind as `linSpace` representing the starting value of the sequence.
    !>  \param[in]  x2          :   The input scalar of the same type and kind as `x1` representing the ending value of the sequence.
    !>  \param[in]  fopen       :   The input scalar of type `logical` of default kind \LK.
    !>                              If `.true.`, the `linSpace` will be <b>f</b>irst-open, meaning that `x1` will NOT be in the output `linSpace` sequence.<br>
    !>                              (**optional**, default = `.false.`)
    !>  \param[in]  lopen       :   The input scalar of type `logical` of default kind \LK.
    !>                              If `.true.`, the `linSpace` will be <b>l</b>ast-open, meaning that `x2` will NOT be in the output `linSpace` sequence.<br>
    !>                              (**optional**, default = `.false.`)
    !>
    !>  \interface{setLinSpace}
    !>  \code{.F90}
    !>
    !>      use pm_arraySpace, only: setLinSpace
    !>
    !>      call setLinSpace(linSpace, x1, x2, fopen = fopen, lopen = lopen)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \remark
    !>  Setting both `fopen = .true.` and `lopen = .true.` will lead to an output `linSpace` whose points are centers of the `count + 1`
    !>  equally-sized bins in the interval `[x1,x2]`.
    !>
    !>  \remark
    !>  If `linSpace` has a size of `1` while `fopen = .false.` and `lopen = .false.`, then `linSpace = [x1]` on return.
    !>
    !>  \see
    !>  [getLinSpace](@ref pm_arraySpace::getLinSpace)<br>
    !>  [getLogSpace](@ref pm_arraySpace::getLogSpace)<br>
    !>  [setLogSpace](@ref pm_arraySpace::setLogSpace)<br>
    !>
    !>  \example{setLinSpace}
    !>  \include{lineno} example/pm_arraySpace/setLinSpace/main.F90
    !>  \compilef{setLinSpace}
    !>  \output{setLinSpace}
    !>  \include{lineno} example/pm_arraySpace/setLinSpace/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arraySpace](@ref test_pm_arraySpace)
    !>
    !>  \final{setLinSpace}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 08:49 PM, August 10, 2021, Dallas, TX
    interface setLinSpace

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module subroutine setLinSpace_CK5(linSpace, x1, x2, fopen, lopen)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLinSpace_CK5
#endif
        use pm_kind, only: CKG => CK5
        implicit none
        complex(CKG), intent(out)   , contiguous    :: linSpace(0:)
        complex(CKG), intent(in)                    :: x1
        complex(CKG), intent(in)                    :: x2
        logical(LK) , intent(in)    , optional      :: fopen, lopen
    end subroutine
#endif

#if CK4_ENABLED
    pure module subroutine setLinSpace_CK4(linSpace, x1, x2, fopen, lopen)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLinSpace_CK4
#endif
        use pm_kind, only: CKG => CK4
        implicit none
        complex(CKG), intent(out)   , contiguous    :: linSpace(0:)
        complex(CKG), intent(in)                    :: x1
        complex(CKG), intent(in)                    :: x2
        logical(LK) , intent(in)    , optional      :: fopen, lopen
    end subroutine
#endif

#if CK3_ENABLED
    pure module subroutine setLinSpace_CK3(linSpace, x1, x2, fopen, lopen)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLinSpace_CK3
#endif
        use pm_kind, only: CKG => CK3
        implicit none
        complex(CKG), intent(out)   , contiguous    :: linSpace(0:)
        complex(CKG), intent(in)                    :: x1
        complex(CKG), intent(in)                    :: x2
        logical(LK) , intent(in)    , optional      :: fopen, lopen
    end subroutine
#endif

#if CK2_ENABLED
    pure module subroutine setLinSpace_CK2(linSpace, x1, x2, fopen, lopen)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLinSpace_CK2
#endif
        use pm_kind, only: CKG => CK2
        implicit none
        complex(CKG), intent(out)   , contiguous    :: linSpace(0:)
        complex(CKG), intent(in)                    :: x1
        complex(CKG), intent(in)                    :: x2
        logical(LK) , intent(in)    , optional      :: fopen, lopen
    end subroutine
#endif

#if CK1_ENABLED
    pure module subroutine setLinSpace_CK1(linSpace, x1, x2, fopen, lopen)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLinSpace_CK1
#endif
        use pm_kind, only: CKG => CK1
        implicit none
        complex(CKG), intent(out)   , contiguous    :: linSpace(0:)
        complex(CKG), intent(in)                    :: x1
        complex(CKG), intent(in)                    :: x2
        logical(LK) , intent(in)    , optional      :: fopen, lopen
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module subroutine setLinSpace_RK5(linSpace, x1, x2, fopen, lopen)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLinSpace_RK5
#endif
        use pm_kind, only: RKG => RK5
        implicit none
        real(RKG)   , intent(out)   , contiguous    :: linSpace(0:)
        real(RKG)   , intent(in)                    :: x1
        real(RKG)   , intent(in)                    :: x2
        logical(LK) , intent(in)    , optional      :: fopen, lopen
    end subroutine
#endif

#if RK4_ENABLED
    pure module subroutine setLinSpace_RK4(linSpace, x1, x2, fopen, lopen)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLinSpace_RK4
#endif
        use pm_kind, only: RKG => RK4
        implicit none
        real(RKG)   , intent(out)   , contiguous    :: linSpace(0:)
        real(RKG)   , intent(in)                    :: x1
        real(RKG)   , intent(in)                    :: x2
        logical(LK) , intent(in)    , optional      :: fopen, lopen
    end subroutine
#endif

#if RK3_ENABLED
    pure module subroutine setLinSpace_RK3(linSpace, x1, x2, fopen, lopen)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLinSpace_RK3
#endif
        use pm_kind, only: RKG => RK3
        implicit none
        real(RKG)   , intent(out)   , contiguous    :: linSpace(0:)
        real(RKG)   , intent(in)                    :: x1
        real(RKG)   , intent(in)                    :: x2
        logical(LK) , intent(in)    , optional      :: fopen, lopen
    end subroutine
#endif

#if RK2_ENABLED
    pure module subroutine setLinSpace_RK2(linSpace, x1, x2, fopen, lopen)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLinSpace_RK2
#endif
        use pm_kind, only: RKG => RK2
        implicit none
        real(RKG)   , intent(out)   , contiguous    :: linSpace(0:)
        real(RKG)   , intent(in)                    :: x1
        real(RKG)   , intent(in)                    :: x2
        logical(LK) , intent(in)    , optional      :: fopen, lopen
    end subroutine
#endif

#if RK1_ENABLED
    pure module subroutine setLinSpace_RK1(linSpace, x1, x2, fopen, lopen)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLinSpace_RK1
#endif
        use pm_kind, only: RKG => RK1
        implicit none
        real(RKG)   , intent(out)   , contiguous    :: linSpace(0:)
        real(RKG)   , intent(in)                    :: x1
        real(RKG)   , intent(in)                    :: x2
        logical(LK) , intent(in)    , optional      :: fopen, lopen
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate `count` evenly-logarithmically-spaced points over the interval `[base**logx1, base**logx2]` if `logx1 < logx2`,
    !>  or over the interval `[logx2, logx1]` if `logx2 < logx1`.
    !>
    !>  \param[in]  logx1   :   The input scalar of type `real` of kind \RKALL representing the
    !>                          starting exponent such that `base**logx1` is the starting value of the sequence.
    !>  \param[in]  logx2   :   The input scalar of the same type and kind as `logx1` representing the
    !>                          ending exponent such that `base**logx1` is the last value of the sequence.
    !>  \param[in]  count   :   The input scalar of type `integer` of default kind \IK representing the length of `logSpace` to generate.
    !>  \param[in]  fopen   :   The input scalar of type `logical` of default kind \LK.
    !>                          If `.true.`, the `logSpace` will be <b>f</b>irst-open, meaning that `logx1` will NOT be in the output `logSpace` sequence<br>
    !>                          (**optional**, default = `.false.`)
    !>  \param[in]  lopen   :   The input scalar of type `logical` of default kind \LK.
    !>                          If `.true.`, the `logSpace` will be <b>l</b>ast-open, meaning that `logx2` will NOT be in the output `logSpace` sequence<br>
    !>                          (**optional**, default = `.false.`)
    !>  \param[in]  base    :   The input scalar of the same type and kind as `logx1` representing the base of the logarithm<br>
    !>                          (**optional**, default = `exp(1)`, that is, Euler`s number).
    !>
    !>  \return
    !>  `logSpace`          :   The output array of shape `(1:count)` of the same type and kind as the input `logx1`
    !>                          containing the log-evenly-spaced sequence within the interval specified by `logx1` and `logx2`.
    !>
    !>  \interface{getLogSpace}
    !>  \code{.F90}
    !>
    !>      use pm_arraySpace, only: getLogSpace
    !>
    !>      logSpace = getLogSpace(logx1, logx2, count, fopen = fopen, lopen = lopen, base = base)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < count` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  Setting both `fopen = .true.` and `lopen = .true.` will lead to an output `logSpace` whose points are centers of the `count`
    !>  equally-sized bins in the interval `[logx1,logx2]`.
    !>
    !>  \see
    !>  [setLogSpace](@ref pm_arraySpace::setLogSpace)<br>
    !>  [getLinSpace](@ref pm_arraySpace::getLinSpace)<br>
    !>  [setLinSpace](@ref pm_arraySpace::setLinSpace)<br>
    !>
    !>  \example{getLogSpace}
    !>  \include{lineno} example/pm_arraySpace/getLogSpace/main.F90
    !>  \compilef{getLogSpace}
    !>  \output{getLogSpace}
    !>  \include{lineno} example/pm_arraySpace/getLogSpace/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arraySpace](@ref test_pm_arraySpace)
    !>
    !>  \todo
    !>  \pmed This generic interface can be expanded to include procedures with real-valued input argument `space` instead of `count`.
    !>
    !>  \final{getLogSpace}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 11:34 PM, August 10, 2021, Dallas, TX
    interface getLogSpace

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getLogSpace_CK5(logx1, logx2, count, fopen, lopen, base) result(logSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSpace_CK5
#endif
        use pm_kind, only: CKG => CK5
        implicit none
        complex(CKG), intent(in)            :: logx1
        complex(CKG), intent(in)            :: logx2
        integer(IK) , intent(in)            :: count
        logical(LK) , intent(in), optional  :: fopen,lopen
        real(CKG)   , intent(in), optional  :: base
        complex(CKG)                        :: logSpace(count)
    end function
#endif

#if CK4_ENABLED
    PURE module function getLogSpace_CK4(logx1, logx2, count, fopen, lopen, base) result(logSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSpace_CK4
#endif
        use pm_kind, only: CKG => CK4
        implicit none
        complex(CKG), intent(in)            :: logx1
        complex(CKG), intent(in)            :: logx2
        integer(IK) , intent(in)            :: count
        logical(LK) , intent(in), optional  :: fopen,lopen
        real(CKG)   , intent(in), optional  :: base
        complex(CKG)                        :: logSpace(count)
    end function
#endif

#if CK3_ENABLED
    PURE module function getLogSpace_CK3(logx1, logx2, count, fopen, lopen, base) result(logSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSpace_CK3
#endif
        use pm_kind, only: CKG => CK3
        implicit none
        complex(CKG), intent(in)            :: logx1
        complex(CKG), intent(in)            :: logx2
        integer(IK) , intent(in)            :: count
        logical(LK) , intent(in), optional  :: fopen,lopen
        real(CKG)   , intent(in), optional  :: base
        complex(CKG)                        :: logSpace(count)
    end function
#endif

#if CK2_ENABLED
    PURE module function getLogSpace_CK2(logx1, logx2, count, fopen, lopen, base) result(logSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSpace_CK2
#endif
        use pm_kind, only: CKG => CK2
        implicit none
        complex(CKG), intent(in)            :: logx1
        complex(CKG), intent(in)            :: logx2
        integer(IK) , intent(in)            :: count
        logical(LK) , intent(in), optional  :: fopen,lopen
        real(CKG)   , intent(in), optional  :: base
        complex(CKG)                        :: logSpace(count)
    end function
#endif

#if CK1_ENABLED
    PURE module function getLogSpace_CK1(logx1, logx2, count, fopen, lopen, base) result(logSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSpace_CK1
#endif
        use pm_kind, only: CKG => CK1
        implicit none
        complex(CKG), intent(in)            :: logx1
        complex(CKG), intent(in)            :: logx2
        integer(IK) , intent(in)            :: count
        logical(LK) , intent(in), optional  :: fopen,lopen
        real(CKG)   , intent(in), optional  :: base
        complex(CKG)                        :: logSpace(count)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getLogSpace_RK5(logx1, logx2, count, fopen, lopen, base) result(logSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSpace_RK5
#endif
        use pm_kind, only: RKG => RK5
        implicit none
        real(RKG)   , intent(in)            :: logx1
        real(RKG)   , intent(in)            :: logx2
        integer(IK) , intent(in)            :: count
        logical(LK) , intent(in), optional  :: fopen,lopen
        real(RKG)   , intent(in), optional  :: base
        real(RKG)                           :: logSpace(count)
    end function
#endif

#if RK4_ENABLED
    PURE module function getLogSpace_RK4(logx1, logx2, count, fopen, lopen, base) result(logSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSpace_RK4
#endif
        use pm_kind, only: RKG => RK4
        implicit none
        real(RKG)   , intent(in)            :: logx1
        real(RKG)   , intent(in)            :: logx2
        integer(IK) , intent(in)            :: count
        logical(LK) , intent(in), optional  :: fopen,lopen
        real(RKG)   , intent(in), optional  :: base
        real(RKG)                           :: logSpace(count)
    end function
#endif

#if RK3_ENABLED
    PURE module function getLogSpace_RK3(logx1, logx2, count, fopen, lopen, base) result(logSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSpace_RK3
#endif
        use pm_kind, only: RKG => RK3
        implicit none
        real(RKG)   , intent(in)            :: logx1
        real(RKG)   , intent(in)            :: logx2
        integer(IK) , intent(in)            :: count
        logical(LK) , intent(in), optional  :: fopen,lopen
        real(RKG)   , intent(in), optional  :: base
        real(RKG)                           :: logSpace(count)
    end function
#endif

#if RK2_ENABLED
    PURE module function getLogSpace_RK2(logx1, logx2, count, fopen, lopen, base) result(logSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSpace_RK2
#endif
        use pm_kind, only: RKG => RK2
        implicit none
        real(RKG)   , intent(in)            :: logx1
        real(RKG)   , intent(in)            :: logx2
        integer(IK) , intent(in)            :: count
        logical(LK) , intent(in), optional  :: fopen,lopen
        real(RKG)   , intent(in), optional  :: base
        real(RKG)                           :: logSpace(count)
    end function
#endif

#if RK1_ENABLED
    PURE module function getLogSpace_RK1(logx1, logx2, count, fopen, lopen, base) result(logSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSpace_RK1
#endif
        use pm_kind, only: RKG => RK1
        implicit none
        real(RKG)   , intent(in)            :: logx1
        real(RKG)   , intent(in)            :: logx2
        integer(IK) , intent(in)            :: count
        logical(LK) , intent(in), optional  :: fopen,lopen
        real(RKG)   , intent(in), optional  :: base
        real(RKG)                           :: logSpace(count)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the `logSpace` output argument with `size(logSpace)` elements of logarithmically-evenly-spaced values over
    !>  the interval `[base**logx1, base**logx2]` if `logx1 < logx2`, or over the interval `[logx2, logx1]` if `logx2 < logx1`.
    !>
    !>  \param[out] logSpace    :   The output `contiguous` array of shape `(:)` of type either `real` of kind \RKALL or `complex` of kind \CKALL
    !>                              containing the log-evenly-spaced sequence within the interval specified by `logx1` and `logx2`.
    !>  \param[in]  logx1       :   The input scalar of the same type and kind as `logSpace` representing the
    !>                              starting exponent such that `base**logx1` is the starting value of the sequence.
    !>  \param[in]  logx2       :   The input scalar of the same type and kind as `logSpace` representing the
    !>                              ending exponent such that `base**logx1` is the last value of the sequence.
    !>  \param[in]  fopen       :   The input scalar of type `logical` of default kind \LK.
    !>                              If `.true.`, the `logSpace` will be <b>f</b>irst-open, meaning that `logx1` will NOT be in the output `logSpace` sequence<br>
    !>                              (**optional**, default = `.false.`)
    !>  \param[in]  lopen       :   The input scalar of type `logical` of default kind \LK.
    !>                              If `.true.`, the `logSpace` will be <b>l</b>ast-open, meaning that `logx2` will NOT be in the output `logSpace` sequence<br>
    !>                              (**optional**, default = `.false.`)
    !>  \param[in]  base        :   The input scalar of type `real` but the same kind as `logSpace` representing the base of the logarithm<br>
    !>                              (**optional**, default = `exp(1)`, that is, Euler`s number).
    !>
    !>  \interface{setLogSpace}
    !>  \code{.F90}
    !>
    !>      use pm_arraySpace, only: setLogSpace
    !>
    !>      call setLogSpace(logSpace, logx1, logx2, fopen = fopen, lopen = lopen, base = base)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \remark
    !>  Setting both `fopen = .true.` and `lopen = .true.` will lead to an output `logSpace` whose points are centers of the `count`
    !>  equally-sized bins in the interval `[logx1,logx2]`.
    !>
    !>  \see
    !>  [getLogSpace](@ref pm_arraySpace::getLogSpace)<br>
    !>  [getLinSpace](@ref pm_arraySpace::getLinSpace)<br>
    !>  [setLinSpace](@ref pm_arraySpace::setLinSpace)<br>
    !>
    !>  \example{setLogSpace}
    !>  \include{lineno} example/pm_arraySpace/setLogSpace/main.F90
    !>  \compilef{setLogSpace}
    !>  \output{setLogSpace}
    !>  \include{lineno} example/pm_arraySpace/setLogSpace/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arraySpace](@ref test_pm_arraySpace)
    !>
    !>  \todo
    !>  \plow This generic interface can be expanded to include procedures with real-valued input argument `space` instead of `count`.
    !>
    !>  \final{setLogSpace}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 11:34 PM, August 10, 2021, Dallas, TX
    interface setLogSpace

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module subroutine setLogSpace_CK5(logSpace, logx1, logx2, fopen, lopen, base)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogSpace_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(out)   , contiguous    :: logSpace(:)
        complex(CKG), intent(in)                    :: logx1
        complex(CKG), intent(in)                    :: logx2
        logical(LK) , intent(in)    , optional      :: fopen,lopen
        real(CKG)   , intent(in)    , optional      :: base
    end subroutine
#endif

#if CK4_ENABLED
    pure module subroutine setLogSpace_CK4(logSpace, logx1, logx2, fopen, lopen, base)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogSpace_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(out)   , contiguous    :: logSpace(:)
        complex(CKG), intent(in)                    :: logx1
        complex(CKG), intent(in)                    :: logx2
        logical(LK) , intent(in)    , optional      :: fopen,lopen
        real(CKG)   , intent(in)    , optional      :: base
    end subroutine
#endif

#if CK3_ENABLED
    pure module subroutine setLogSpace_CK3(logSpace, logx1, logx2, fopen, lopen, base)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogSpace_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(out)   , contiguous    :: logSpace(:)
        complex(CKG), intent(in)                    :: logx1
        complex(CKG), intent(in)                    :: logx2
        logical(LK) , intent(in)    , optional      :: fopen,lopen
        real(CKG)   , intent(in)    , optional      :: base
    end subroutine
#endif

#if CK2_ENABLED
    pure module subroutine setLogSpace_CK2(logSpace, logx1, logx2, fopen, lopen, base)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogSpace_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(out)   , contiguous    :: logSpace(:)
        complex(CKG), intent(in)                    :: logx1
        complex(CKG), intent(in)                    :: logx2
        logical(LK) , intent(in)    , optional      :: fopen,lopen
        real(CKG)   , intent(in)    , optional      :: base
    end subroutine
#endif

#if CK1_ENABLED
    pure module subroutine setLogSpace_CK1(logSpace, logx1, logx2, fopen, lopen, base)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogSpace_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(out)   , contiguous    :: logSpace(:)
        complex(CKG), intent(in)                    :: logx1
        complex(CKG), intent(in)                    :: logx2
        logical(LK) , intent(in)    , optional      :: fopen,lopen
        real(CKG)   , intent(in)    , optional      :: base
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module subroutine setLogSpace_RK5(logSpace, logx1, logx2, fopen, lopen, base)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogSpace_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)   , contiguous    :: logSpace(:)
        real(RKG)   , intent(in)                    :: logx1
        real(RKG)   , intent(in)                    :: logx2
        logical(LK) , intent(in)    , optional      :: fopen,lopen
        real(RKG)   , intent(in)    , optional      :: base
    end subroutine
#endif

#if RK4_ENABLED
    pure module subroutine setLogSpace_RK4(logSpace, logx1, logx2, fopen, lopen, base)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogSpace_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)   , contiguous    :: logSpace(:)
        real(RKG)   , intent(in)                    :: logx1
        real(RKG)   , intent(in)                    :: logx2
        logical(LK) , intent(in)    , optional      :: fopen,lopen
        real(RKG)   , intent(in)    , optional      :: base
    end subroutine
#endif

#if RK3_ENABLED
    pure module subroutine setLogSpace_RK3(logSpace, logx1, logx2, fopen, lopen, base)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogSpace_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)   , contiguous    :: logSpace(:)
        real(RKG)   , intent(in)                    :: logx1
        real(RKG)   , intent(in)                    :: logx2
        logical(LK) , intent(in)    , optional      :: fopen,lopen
        real(RKG)   , intent(in)    , optional      :: base
    end subroutine
#endif

#if RK2_ENABLED
    pure module subroutine setLogSpace_RK2(logSpace, logx1, logx2, fopen, lopen, base)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogSpace_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)   , contiguous    :: logSpace(:)
        real(RKG)   , intent(in)                    :: logx1
        real(RKG)   , intent(in)                    :: logx2
        logical(LK) , intent(in)    , optional      :: fopen,lopen
        real(RKG)   , intent(in)    , optional      :: base
    end subroutine
#endif

#if RK1_ENABLED
    pure module subroutine setLogSpace_RK1(logSpace, logx1, logx2, fopen, lopen, base)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogSpace_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)   , contiguous    :: logSpace(:)
        real(RKG)   , intent(in)                    :: logx1
        real(RKG)   , intent(in)                    :: logx2
        logical(LK) , intent(in)    , optional      :: fopen,lopen
        real(RKG)   , intent(in)    , optional      :: base
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arraySpace ! LCOV_EXCL_LINE