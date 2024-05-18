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
!>  This module contains procedures and generic interfaces for flattening
!>  (duplicating the elements of) an array according to a user-specified `weight`.
!>
!>  \see
!>  [pm_arrayCompact](@ref pm_arrayCompact)<br>
!>
!>  \test
!>  [test_pm_arrayVerbose](@ref test_pm_arrayVerbose)
!>
!>  \todo
!>  \pmed
!>  A generic subroutine interface corresponding to the generic function interface [getVerbose](@ref pm_arrayVerbose::getVerbose) might be added in future.<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Saturday 1:30 AM, August 20, 2016, Institute for Computational Engineering and Sciences, UT Austin, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayVerbose

!>  \cond excluded
!   \bug
!   The following bypasses the bug reported below that creates a conflict between Intel and gfortran.
#if     __INTEL_COMPILER
#define LEN_STR :
#else
#define LEN_STR len(array,IK)
#endif
#if 0
#define CONTIGUOUS__, contiguous
#else
#define CONTIGUOUS__
#endif
!>  \endcond excluded

    use pm_kind, only: SK, IK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_arrayVerbose"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate an equally-weighted (verbose or flattened) array of the input weighted `array` of rank 1 or 2.
    !>
    !>  \details
    !>  The output object contains elements of the input `array` duplicated as many times as indicated by the corresponding element of `weight`.
    !>
    !>  \param[in]  array       :   The input `contiguous` array of shape `(1:nsam)` or `(1:ndim, 1:nsam)` of either
    !>                              <ul>
    !>                                  <li>    type `character` of kind \SKALL or <br>
    !>                                  <li>    type `integer` of kind \IKALL or <br>
    !>                                  <li>    type `logical` of kind \LKALL or <br>
    !>                                  <li>    type `complex` of kind \CKALL or <br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ul>
    !>                              or,
    !>                              <ul>
    !>                                  <li>    a scalar `character` of kind \SKALL,
    !>                              </ul>
    !>                              containing the weighted elements to be flattened.<br>
    !>  \param[in]  weight      :   The input `contiguous` array of shape `(1:nsam)` of type `integer` of default kind \IK of the same size as `size(array, dim)`, i.e., along the axis `array` must be flattened.<br>
    !>                              It contains the weights of elements of `array`.<br>
    !>                              The values of the elements of `weight` are allowed to be negative or zero, in which case,
    !>                              the corresponding elements will be excluded from the output `verbose`.
    !>  \param[in]  weisum      :   The input scalar of type `integer` of default kind \IK representing the sum of all **non-negative** elements of `weight`.<br>
    !>                              The value of `weisum` can be readily obtained via `sum(weight, mask = weight > 0)`.
    !>  \param[in]  dim         :   The input scalar of type `integer` of default kind \IK representing the axis of `array(:,:)` along which `array` must be flattened.<br>
    !>                              (**optional**, it must be present <b>if and only if</b> `array` is of shape `(:,:)`. If present, `size(weight) == size(array, dim)` must hold.)
    !>
    !>  \return
    !>  `verbose`               :   The output array of the same type, kind, and shape as the input `array` containing the flattened version of the input `array`.<br>
    !>                              <ul>
    !>                                  <li>    If `array` is of rank `2`, then the condition `size(verbose, dim) == weisum .and. size(verbose, 3 - dim) == size(array, 3 - dim)` holds.<br>
    !>                                  <li>    If `array` is of rank `1`, then the condition `size(verbose) == weisum` holds.<br>
    !>                              </ul>
    !>
    !>  \interface{getVerbose}
    !>  \code{.F90}
    !>
    !>      use pm_arrayVerbose, only: getVerbose
    !>
    !>      verbose = getVerbose(array, weight(1:len(array)), weisum) ! scalar character objects.
    !>      verbose(:) = getVerbose(array(1:size(array)), weight(1:size(array)), weisum) ! all intrinsic array objects.
    !>      verbose(1 : merge(weisum, size(array, 3 - dim), dim == 1), 1 : merge(weisum, size(array, 3 - dim), dim == 2)) = getVerbose(array(:,:), weight(1:size(array,dim)), weisum, dim) ! all intrinsic array objects.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(weight) == size(array, dim)` must hold for the input arguments `array` and `weight`.<br>
    !>  The condition `weisum == sum(weight, mask = weight > 0)` must hold for the input arguments `weight` and `weisum`.<br>
    !>  The condition `dim == 1 .or. dim == 2` must hold.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getCompact](@ref pm_arrayCompact::getCompact)<br>
    !>  [setCompact](@ref pm_arrayCompact::setCompact)<br>
    !>
    !>  \example{getVerbose}
    !>  \include{lineno} example/pm_arrayVerbose/getVerbose/main.F90
    !>  \compilef{getVerbose}
    !>  \output{getVerbose}
    !>  \include{lineno} example/pm_arrayVerbose/getVerbose/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayVerbose](@ref test_pm_arrayVerbose)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.2.0}, \gfortran{10-11}
    !>  \desc
    !>  See the bug described in [getUnique](@ref pm_arrayUnique::getUnique).<br>
    !>  \remedy
    !>  See the bug described in [getUnique](@ref pm_arrayUnique::getUnique).<br>
    !>
    !>  \final{getVerbose}
    !>
    !>  \author
    !>  \AmirShahmoradi, Saturday 1:30 AM, August 20, 2016, Institute for Computational Engineering and Sciences, UT Austin, TX
    interface getVerbose

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getVerbose_D0_SK5(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)    , intent(in)                    :: array
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        character(weisum,SKC)                               :: verbose
    end function
#endif

#if SK4_ENABLED
    PURE module function getVerbose_D0_SK4(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)    , intent(in)                    :: array
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        character(weisum,SKC)                               :: verbose
    end function
#endif

#if SK3_ENABLED
    PURE module function getVerbose_D0_SK3(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)    , intent(in)                    :: array
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        character(weisum,SKC)                               :: verbose
    end function
#endif

#if SK2_ENABLED
    PURE module function getVerbose_D0_SK2(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)    , intent(in)                    :: array
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        character(weisum,SKC)                               :: verbose
    end function
#endif

#if SK1_ENABLED
    PURE module function getVerbose_D0_SK1(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)    , intent(in)                    :: array
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        character(weisum,SKC)                               :: verbose
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getVerbose_D1_SK5(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)    , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        character(len(array,IK),SKC)                        :: verbose(weisum)
    end function
#endif

#if SK4_ENABLED
    PURE module function getVerbose_D1_SK4(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)    , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        character(len(array,IK),SKC)                        :: verbose(weisum)
    end function
#endif

#if SK3_ENABLED
    PURE module function getVerbose_D1_SK3(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)    , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        character(len(array,IK),SKC)                        :: verbose(weisum)
    end function
#endif

#if SK2_ENABLED
    PURE module function getVerbose_D1_SK2(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)    , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        character(len(array,IK),SKC)                        :: verbose(weisum)
    end function
#endif

#if SK1_ENABLED
    PURE module function getVerbose_D1_SK1(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)    , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        character(len(array,IK),SKC)                        :: verbose(weisum)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getVerbose_D1_IK5(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)        , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        integer(IKC)                                        :: verbose(weisum)
    end function
#endif

#if IK4_ENABLED
    PURE module function getVerbose_D1_IK4(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)        , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        integer(IKC)                                        :: verbose(weisum)
    end function
#endif

#if IK3_ENABLED
    PURE module function getVerbose_D1_IK3(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)        , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        integer(IKC)                                        :: verbose(weisum)
    end function
#endif

#if IK2_ENABLED
    PURE module function getVerbose_D1_IK2(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)        , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        integer(IKC)                                        :: verbose(weisum)
    end function
#endif

#if IK1_ENABLED
    PURE module function getVerbose_D1_IK1(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)        , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        integer(IKC)                                        :: verbose(weisum)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getVerbose_D1_LK5(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)        , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        logical(LKC)                                        :: verbose(weisum)
    end function
#endif

#if LK4_ENABLED
    PURE module function getVerbose_D1_LK4(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)        , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        logical(LKC)                                        :: verbose(weisum)
    end function
#endif

#if LK3_ENABLED
    PURE module function getVerbose_D1_LK3(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)        , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        logical(LKC)                                        :: verbose(weisum)
    end function
#endif

#if LK2_ENABLED
    PURE module function getVerbose_D1_LK2(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)        , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        logical(LKC)                                        :: verbose(weisum)
    end function
#endif

#if LK1_ENABLED
    PURE module function getVerbose_D1_LK1(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)        , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        logical(LKC)                                        :: verbose(weisum)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getVerbose_D1_CK5(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        complex(CKC)                                        :: verbose(weisum)
    end function
#endif

#if CK4_ENABLED
    PURE module function getVerbose_D1_CK4(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        complex(CKC)                                        :: verbose(weisum)
    end function
#endif

#if CK3_ENABLED
    PURE module function getVerbose_D1_CK3(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        complex(CKC)                                        :: verbose(weisum)
    end function
#endif

#if CK2_ENABLED
    PURE module function getVerbose_D1_CK2(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        complex(CKC)                                        :: verbose(weisum)
    end function
#endif

#if CK1_ENABLED
    PURE module function getVerbose_D1_CK1(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        complex(CKC)                                        :: verbose(weisum)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getVerbose_D1_RK5(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        real(RKC)                                           :: verbose(weisum)
    end function
#endif

#if RK4_ENABLED
    PURE module function getVerbose_D1_RK4(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        real(RKC)                                           :: verbose(weisum)
    end function
#endif

#if RK3_ENABLED
    PURE module function getVerbose_D1_RK3(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        real(RKC)                                           :: verbose(weisum)
    end function
#endif

#if RK2_ENABLED
    PURE module function getVerbose_D1_RK2(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        real(RKC)                                           :: verbose(weisum)
    end function
#endif

#if RK1_ENABLED
    PURE module function getVerbose_D1_RK1(array, weight, weisum) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in)    CONTIGUOUS__    :: array(:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum
        real(RKC)                                           :: verbose(weisum)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getVerbose_D2_SK5(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)    , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        character(len(array,IK),SKC)                        :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

#if SK4_ENABLED
    PURE module function getVerbose_D2_SK4(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)    , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        character(len(array,IK),SKC)                        :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

#if SK3_ENABLED
    PURE module function getVerbose_D2_SK3(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)    , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        character(len(array,IK),SKC)                        :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

#if SK2_ENABLED
    PURE module function getVerbose_D2_SK2(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)    , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        character(len(array,IK),SKC)                        :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

#if SK1_ENABLED
    PURE module function getVerbose_D2_SK1(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)    , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        character(len(array,IK),SKC)                        :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getVerbose_D2_IK5(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)        , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        integer(IKC)                                        :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

#if IK4_ENABLED
    PURE module function getVerbose_D2_IK4(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)        , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        integer(IKC)                                        :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

#if IK3_ENABLED
    PURE module function getVerbose_D2_IK3(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)        , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        integer(IKC)                                        :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

#if IK2_ENABLED
    PURE module function getVerbose_D2_IK2(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)        , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        integer(IKC)                                        :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

#if IK1_ENABLED
    PURE module function getVerbose_D2_IK1(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)        , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        integer(IKC)                                        :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getVerbose_D2_LK5(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)        , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        logical(LKC)                                        :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

#if LK4_ENABLED
    PURE module function getVerbose_D2_LK4(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)        , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        logical(LKC)                                        :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

#if LK3_ENABLED
    PURE module function getVerbose_D2_LK3(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)        , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        logical(LKC)                                        :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

#if LK2_ENABLED
    PURE module function getVerbose_D2_LK2(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)        , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        logical(LKC)                                        :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

#if LK1_ENABLED
    PURE module function getVerbose_D2_LK1(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)        , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        logical(LKC)                                        :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getVerbose_D2_CK5(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        complex(CKC)                                        :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getVerbose_D2_CK4(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        complex(CKC)                                        :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getVerbose_D2_CK3(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        complex(CKC)                                        :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getVerbose_D2_CK2(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        complex(CKC)                                        :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getVerbose_D2_CK1(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        complex(CKC)                                        :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getVerbose_D2_RK5(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        real(RKC)                                           :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getVerbose_D2_RK4(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        real(RKC)                                           :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getVerbose_D2_RK3(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        real(RKC)                                           :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getVerbose_D2_RK2(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        real(RKC)                                           :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getVerbose_D2_RK1(array, weight, weisum, dim) result(verbose)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVerbose_D2_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in)    CONTIGUOUS__    :: array(:,:)
        integer(IK)         , intent(in)    CONTIGUOUS__    :: weight(:)
        integer(IK)         , intent(in)                    :: weisum, dim
        real(RKC)                                           :: verbose(merge(weisum, size(array, 3 - dim, IK), dim == 1_IK), merge(weisum, size(array, 3 - dim, IK), dim == 2_IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getVerbose

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayVerbose ! LCOV_EXCL_LINE