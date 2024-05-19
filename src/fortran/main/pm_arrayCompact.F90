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

!>  \cond excluded
!   \bug
!   The following bypasses the bug reported below that creates a conflict between Intel and gfortran.
#if     __INTEL_COMPILER
#define LEN_STR :
#else
#define LEN_STR len(array,IK)
#endif
!>  \endcond excluded

!>  \brief
!>  This module contains procedures and generic interfaces for condensing
!>  (removing duplicate sequential the elements of) an array of arbitrary intrinsic type.<br>
!>
!>  \see
!>  [pm_arrayVerbose](@ref pm_arrayVerbose)<br>
!>
!>  \test
!>  [test_pm_arrayCompact](@ref test_pm_arrayCompact)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Saturday 1:48 AM, August 20, 2016, Institute for Computational Engineering and Sciences, UT Austin, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayCompact

    use pm_kind, only: SK, IK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_arrayCompact"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate a compact version of the input `array` where all **sequentially** duplicate entries along
    !>  the specified dimension of `array` are condensed to a single entry in the output `compact`.
    !>
    !>  \param[in]  array       :   The input `contiguous` array of shape `(:)` or `(:,:)` of either
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL or
    !>                                  <li>    type `integer` of kind \IKALL or
    !>                                  <li>    type `logical` of kind \LKALL or
    !>                                  <li>    type `complex` of kind \CKALL or
    !>                                  <li>    type `real` of kind \RKALL
    !>                              </ol>
    !>                              or,
    !>                              <ol>
    !>                                  <li>    scalar `character` of kind \SKALL,
    !>                              </ol>
    !>                              containing the array that has to be compactified.
    !>  \param[in]  dim         :   The input scalar of type `integer` of default kind \IK representing the axis of `array(:,:)` along which `array` must be compactified.<br>
    !>                              (**optional**, it must be present <b>if and only if</b> `array` is of shape `(:,:)`.)
    !>
    !>  \return
    !>  `compact`               :   The output array of the same type, kind, and shape as the input `array` containing the compactified version of the input `array`.<br>
    !>
    !>  \interface{getCompact}
    !>  \code{.F90}
    !>
    !>      use pm_arrayCompact, only: getCompact
    !>
    !>      compact = getCompact(array) ! scalar character objects.
    !>      compact(:) = getCompact(array(:)) ! all intrinsic array objects.
    !>      compact(:,:) = getCompact(array(:,:), dim) ! all intrinsic array objects.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `dim == 1 .or. dim == 2` must hold.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [setCompact](@ref pm_arrayCompact::setCompact)<br>
    !>  [getVerbose](@ref pm_arrayVerbose::getVerbose)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_arrayCompact/getCompact/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_arrayCompact/getCompact/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayCompact](@ref test_pm_arrayCompact)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{10.3-12}, \ifort{2021-2022}
    !>  \desc
    !>  \ifort and \gfortran share a common bug with opposing behavior.<br>
    !>  \ifort cannot handle assumed-length `allocatable` output arguments of type `character`.<br>
    !>  \gfortran cannot handle deferred-length `allocatable` output arguments of type `character`.<br>
    !>  \remedy
    !>  For now, a preprocessor macro defines two separate interfaces for the two compilers so that both compilers can compile this file.<br>
    !>  This minor interface difference should not impact the usage of this module with different compilers.<br>
    !>
    !>  \final{getCompact}
    !>
    !>  \author
    !>  \AmirShahmoradi, Saturday 1:48 AM, August 20, 2016, Institute for Computational Engineering and Sciences, UT Austin, TX
    interface getCompact

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getCompact_D0_SK5(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        character(:,SKG)            , allocatable                   :: compact
    end function
#endif

#if SK4_ENABLED
    PURE module function getCompact_D0_SK4(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        character(:,SKG)            , allocatable                   :: compact
    end function
#endif

#if SK3_ENABLED
    PURE module function getCompact_D0_SK3(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        character(:,SKG)            , allocatable                   :: compact
    end function
#endif

#if SK2_ENABLED
    PURE module function getCompact_D0_SK2(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        character(:,SKG)            , allocatable                   :: compact
    end function
#endif

#if SK1_ENABLED
    PURE module function getCompact_D0_SK1(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        character(:,SKG)            , allocatable                   :: compact
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getCompact_D1_SK5(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(LEN_STR,SKG)                      , allocatable   :: compact(:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getCompact_D1_SK4(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(LEN_STR,SKG)                      , allocatable   :: compact(:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getCompact_D1_SK3(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(LEN_STR,SKG)                      , allocatable   :: compact(:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getCompact_D1_SK2(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(LEN_STR,SKG)                      , allocatable   :: compact(:)
    end function
#endif

#if SK1_ENABLED
    PURE module function getCompact_D1_SK1(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(LEN_STR,SKG)                      , allocatable   :: compact(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getCompact_D1_IK5(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , allocatable                   :: compact(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getCompact_D1_IK4(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , allocatable                   :: compact(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getCompact_D1_IK3(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , allocatable                   :: compact(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getCompact_D1_IK2(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , allocatable                   :: compact(:)
    end function
#endif

#if IK1_ENABLED
    PURE module function getCompact_D1_IK1(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , allocatable                   :: compact(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getCompact_D1_LK5(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , allocatable                   :: compact(:)
    end function
#endif

#if LK4_ENABLED
    PURE module function getCompact_D1_LK4(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , allocatable                   :: compact(:)
    end function
#endif

#if LK3_ENABLED
    PURE module function getCompact_D1_LK3(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , allocatable                   :: compact(:)
    end function
#endif

#if LK2_ENABLED
    PURE module function getCompact_D1_LK2(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , allocatable                   :: compact(:)
    end function
#endif

#if LK1_ENABLED
    PURE module function getCompact_D1_LK1(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , allocatable                   :: compact(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCompact_D1_CK5(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , allocatable                   :: compact(:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getCompact_D1_CK4(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , allocatable                   :: compact(:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getCompact_D1_CK3(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , allocatable                   :: compact(:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getCompact_D1_CK2(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , allocatable                   :: compact(:)
    end function
#endif

#if CK1_ENABLED
    PURE module function getCompact_D1_CK1(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , allocatable                   :: compact(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCompact_D1_RK5(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , allocatable                   :: compact(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getCompact_D1_RK4(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , allocatable                   :: compact(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getCompact_D1_RK3(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , allocatable                   :: compact(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getCompact_D1_RK2(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , allocatable                   :: compact(:)
    end function
#endif

#if RK1_ENABLED
    PURE module function getCompact_D1_RK1(array) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , allocatable                   :: compact(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getCompact_D2_SK5(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        character(LEN_STR,SKG)                      , allocatable   :: compact(:,:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getCompact_D2_SK4(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        character(LEN_STR,SKG)                      , allocatable   :: compact(:,:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getCompact_D2_SK3(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        character(LEN_STR,SKG)                      , allocatable   :: compact(:,:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getCompact_D2_SK2(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        character(LEN_STR,SKG)                      , allocatable   :: compact(:,:)
    end function
#endif

#if SK1_ENABLED
    PURE module function getCompact_D2_SK1(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        character(LEN_STR,SKG)                      , allocatable   :: compact(:,:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getCompact_D2_IK5(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        integer(IKG)                                , allocatable   :: compact(:,:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getCompact_D2_IK4(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        integer(IKG)                                , allocatable   :: compact(:,:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getCompact_D2_IK3(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        integer(IKG)                                , allocatable   :: compact(:,:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getCompact_D2_IK2(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        integer(IKG)                                , allocatable   :: compact(:,:)
    end function
#endif

#if IK1_ENABLED
    PURE module function getCompact_D2_IK1(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        integer(IKG)                                , allocatable   :: compact(:,:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getCompact_D2_LK5(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        logical(LKG)                                , allocatable   :: compact(:,:)
    end function
#endif

#if LK4_ENABLED
    PURE module function getCompact_D2_LK4(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        logical(LKG)                                , allocatable   :: compact(:,:)
    end function
#endif

#if LK3_ENABLED
    PURE module function getCompact_D2_LK3(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        logical(LKG)                                , allocatable   :: compact(:,:)
    end function
#endif

#if LK2_ENABLED
    PURE module function getCompact_D2_LK2(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        logical(LKG)                                , allocatable   :: compact(:,:)
    end function
#endif

#if LK1_ENABLED
    PURE module function getCompact_D2_LK1(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        logical(LKG)                                , allocatable   :: compact(:,:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCompact_D2_CK5(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        complex(CKG)                                , allocatable   :: compact(:,:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getCompact_D2_CK4(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        complex(CKG)                                , allocatable   :: compact(:,:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getCompact_D2_CK3(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        complex(CKG)                                , allocatable   :: compact(:,:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getCompact_D2_CK2(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        complex(CKG)                                , allocatable   :: compact(:,:)
    end function
#endif

#if CK1_ENABLED
    PURE module function getCompact_D2_CK1(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        complex(CKG)                                , allocatable   :: compact(:,:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCompact_D2_RK5(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        real(RKG)                                   , allocatable   :: compact(:,:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getCompact_D2_RK4(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        real(RKG)                                   , allocatable   :: compact(:,:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getCompact_D2_RK3(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        real(RKG)                                   , allocatable   :: compact(:,:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getCompact_D2_RK2(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        real(RKG)                                   , allocatable   :: compact(:,:)
    end function
#endif

#if RK1_ENABLED
    PURE module function getCompact_D2_RK1(array, dim) result(compact)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompact_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)                    :: dim
        real(RKG)                                   , allocatable   :: compact(:,:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getCompact

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate a compacted version of the input `array` where all **sequentially** duplicate entries along
    !>  the specified dimension of `array` are condensed to a single entry in the output `compact`.
    !>
    !>  \param[inout]   array       :   The input `contiguous` array of shape `(:)` or `(:,:)` of either
    !>                                  <ul>
    !>                                      <li>    type `character` of kind \SKALL or
    !>                                      <li>    type `integer` of kind \IKALL or
    !>                                      <li>    type `logical` of kind \LKALL or
    !>                                      <li>    type `complex` of kind \CKALL or
    !>                                      <li>    type `real` of kind \RKALL
    !>                                  </ul>
    !>                                  or scalar `character` of kind \SKALL.<br>
    !>                                  On input, it contains the elements to be condensed.<br>
    !>                                  On output, the first `csize` entries of array contain the condensed elements of the original input `array`,
    !>                                  such that **no** two adjacent entries in `array` from entry `1` to entry `csize` are duplicates.
    !>  \param[out] weight          :   The output `contiguous` array of shape `(:)` of type `integer` of default kind \IK of the same size as `array`.<br>
    !>                                  On output, the first `csize` elements of `weight` contain the weights of the elements of the output condensed `array` such that
    !>                                  [getVerbose(array(1:csize), weight(1:csize), sum(weight))](@ref pm_arrayVerbose::getVerbose) returns the original `array`.<br>
    !>  \param[out] csize           :   The output scalar of type `integer` of default kind \IK such that,<br>
    !>                                  <ul>
    !>                                      <li>    `weight(1:csize)` contains the weight of the condensed elements of the input `array`, and<br>
    !>                                      <li>    `array(1:csize)` or `array(:,1:csize)` or `array(1:csize,:)` contains the condensed elements of the input `array`.<br>
    !>                                  </ul>
    !>  \param[in]  dim             :   The input scalar of type `integer` of default kind \IK representing the axis of `array(:,:)` along which `array` must be compacted.<br>
    !>                                  (**optional**, it must be present <b>if and only if</b> `array` is of shape `(:,:)`.)
    !>
    !>  \interface{setCompact}
    !>  \code{.F90}
    !>
    !>      use pm_arrayCompact, only: setCompact
    !>
    !>      call setCompact(array, weight(:), csize) ! scalar character objects.
    !>      call setCompact(array(:), weight(:), csize) ! all intrinsic array objects.
    !>      call setCompact(array(:,:), weight(:), csize, dim) ! all intrinsic array objects.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `dim == 1 .or. dim == 2` must hold.<br>
    !>  The condition `size(weight) == size(array)` must hold when `array` is rank `1`.<br>
    !>  The condition `size(weight) == size(array, dim)` must hold when `array` is rank `2`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getCompact](@ref pm_arrayCompact::getCompact)<br>
    !>  [getVerbose](@ref pm_arrayVerbose::getVerbose)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_arrayCompact/setCompact/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_arrayCompact/setCompact/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayCompact](@ref test_pm_arrayCompact)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Saturday 1:48 AM, August 20, 2016, Institute for Computational Engineering and Sciences, UT Austin, TX
    interface setCompact

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setCompact_D0_SK5(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setCompact_D0_SK4(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setCompact_D0_SK3(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setCompact_D0_SK2(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setCompact_D0_SK1(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setCompact_D1_SK5(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setCompact_D1_SK4(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setCompact_D1_SK3(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setCompact_D1_SK2(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setCompact_D1_SK1(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setCompact_D1_IK5(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setCompact_D1_IK4(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setCompact_D1_IK3(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setCompact_D1_IK2(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setCompact_D1_IK1(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setCompact_D1_LK5(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setCompact_D1_LK4(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setCompact_D1_LK3(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setCompact_D1_LK2(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setCompact_D1_LK1(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCompact_D1_CK5(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCompact_D1_CK4(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCompact_D1_CK3(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCompact_D1_CK2(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCompact_D1_CK1(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCompact_D1_RK5(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCompact_D1_RK4(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCompact_D1_RK3(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCompact_D1_RK2(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCompact_D1_RK1(array, weight, csize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setCompact_D2_SK5(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setCompact_D2_SK4(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setCompact_D2_SK3(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setCompact_D2_SK2(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setCompact_D2_SK1(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setCompact_D2_IK5(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setCompact_D2_IK4(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setCompact_D2_IK3(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setCompact_D2_IK2(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setCompact_D2_IK1(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setCompact_D2_LK5(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setCompact_D2_LK4(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setCompact_D2_LK3(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setCompact_D2_LK2(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setCompact_D2_LK1(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCompact_D2_CK5(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCompact_D2_CK4(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCompact_D2_CK3(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCompact_D2_CK2(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCompact_D2_CK1(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCompact_D2_RK5(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCompact_D2_RK4(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCompact_D2_RK3(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCompact_D2_RK2(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCompact_D2_RK1(array, weight, csize, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCompact_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(out)   , contiguous    :: weight(:)
        integer(IK)                 , intent(out)                   :: csize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#if SK5_ENABLED
!    PURE module subroutine setCompactNew_D0_SK5(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D0_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        character(*,SKG)            , intent(in)                    :: array
!        integer(IK)                 , intent(in)    , allocatable   :: weight(:)
!        character(:,SKG)            , intent(out)   , allocatable   :: compact
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setCompactNew_D0_SK4(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D0_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        character(*,SKG)            , intent(in)                    :: array
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        character(:,SKG)            , intent(out)   , allocatable   :: compact
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setCompactNew_D0_SK3(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D0_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        character(*,SKG)            , intent(in)                    :: array
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        character(:,SKG)            , intent(out)   , allocatable   :: compact
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setCompactNew_D0_SK2(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D0_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        character(*,SKG)            , intent(in)                    :: array
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        character(:,SKG)            , intent(out)   , allocatable   :: compact
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setCompactNew_D0_SK1(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D0_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        character(*,SKG)            , intent(in)                    :: array
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        character(:,SKG)            , intent(out)   , allocatable   :: compact
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module subroutine setCompactNew_D1_SK5(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        character(len(array),SKG)   , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setCompactNew_D1_SK4(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        character(len(array),SKG)   , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setCompactNew_D1_SK3(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        character(len(array),SKG)   , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setCompactNew_D1_SK2(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        character(len(array),SKG)   , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setCompactNew_D1_SK1(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        character(len(array),SKG)   , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module subroutine setCompactNew_D1_IK5(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_IK5
!#endif
!        use pm_kind, only: IKG => IK5
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IKG)                , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    PURE module subroutine setCompactNew_D1_IK4(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_IK4
!#endif
!        use pm_kind, only: IKG => IK4
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IKG)                , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    PURE module subroutine setCompactNew_D1_IK3(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_IK3
!#endif
!        use pm_kind, only: IKG => IK3
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IKG)                , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    PURE module subroutine setCompactNew_D1_IK2(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_IK2
!#endif
!        use pm_kind, only: IKG => IK2
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IKG)                , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    PURE module subroutine setCompactNew_D1_IK1(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_IK1
!#endif
!        use pm_kind, only: IKG => IK1
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IKG)                , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if LK5_ENABLED
!    PURE module subroutine setCompactNew_D1_LK5(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_LK5
!#endif
!        use pm_kind, only: LKG => LK5
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        logical(LKG)                , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!#if LK4_ENABLED
!    PURE module subroutine setCompactNew_D1_LK4(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_LK4
!#endif
!        use pm_kind, only: LKG => LK4
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        logical(LKG)                , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!#if LK3_ENABLED
!    PURE module subroutine setCompactNew_D1_LK3(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_LK3
!#endif
!        use pm_kind, only: LKG => LK3
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        logical(LKG)                , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!#if LK2_ENABLED
!    PURE module subroutine setCompactNew_D1_LK2(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_LK2
!#endif
!        use pm_kind, only: LKG => LK2
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        logical(LKG)                , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!#if LK1_ENABLED
!    PURE module subroutine setCompactNew_D1_LK1(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_LK1
!#endif
!        use pm_kind, only: LKG => LK1
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        logical(LKG)                , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if CK5_ENABLED
!    PURE module subroutine setCompactNew_D1_CK5(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_CK5
!#endif
!        use pm_kind, only: CKG => CK5
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        complex(CKG)                , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!#if CK4_ENABLED
!    PURE module subroutine setCompactNew_D1_CK4(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_CK4
!#endif
!        use pm_kind, only: CKG => CK4
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        complex(CKG)                , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!#if CK3_ENABLED
!    PURE module subroutine setCompactNew_D1_CK3(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_CK3
!#endif
!        use pm_kind, only: CKG => CK3
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        complex(CKG)                , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!#if CK2_ENABLED
!    PURE module subroutine setCompactNew_D1_CK2(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_CK2
!#endif
!        use pm_kind, only: CKG => CK2
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        complex(CKG)                , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!#if CK1_ENABLED
!    PURE module subroutine setCompactNew_D1_CK1(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_CK1
!#endif
!        use pm_kind, only: CKG => CK1
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        complex(CKG)                , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module subroutine setCompactNew_D1_RK5(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        real(RKG)                   , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE module subroutine setCompactNew_D1_RK4(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        real(RKG)                   , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE module subroutine setCompactNew_D1_RK3(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        real(RKG)                   , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE module subroutine setCompactNew_D1_RK2(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        real(RKG)                   , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE module subroutine setCompactNew_D1_RK1(array, weight, compact)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D1_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        real(RKG)                   , intent(out)   , allocatable   :: compact(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module subroutine setCompactNew_D2_SK5(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        character(len(array),SKG)   , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setCompactNew_D2_SK4(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        character(len(array),SKG)   , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setCompactNew_D2_SK3(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        character(len(array),SKG)   , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setCompactNew_D2_SK2(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        character(len(array),SKG)   , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setCompactNew_D2_SK1(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        character(len(array),SKG)   , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module subroutine setCompactNew_D2_IK5(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_IK5
!#endif
!        use pm_kind, only: IKG => IK5
!        integer(IKG)                , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        integer(IKG)                , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    PURE module subroutine setCompactNew_D2_IK4(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_IK4
!#endif
!        use pm_kind, only: IKG => IK4
!        integer(IKG)                , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        integer(IKG)                , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    PURE module subroutine setCompactNew_D2_IK3(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_IK3
!#endif
!        use pm_kind, only: IKG => IK3
!        integer(IKG)                , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        integer(IKG)                , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    PURE module subroutine setCompactNew_D2_IK2(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_IK2
!#endif
!        use pm_kind, only: IKG => IK2
!        integer(IKG)                , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        integer(IKG)                , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    PURE module subroutine setCompactNew_D2_IK1(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_IK1
!#endif
!        use pm_kind, only: IKG => IK1
!        integer(IKG)                , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        integer(IKG)                , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if LK5_ENABLED
!    PURE module subroutine setCompactNew_D2_LK5(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_LK5
!#endif
!        use pm_kind, only: LKG => LK5
!        logical(LKG)                , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        logical(LKG)                , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!#if LK4_ENABLED
!    PURE module subroutine setCompactNew_D2_LK4(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_LK4
!#endif
!        use pm_kind, only: LKG => LK4
!        logical(LKG)                , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        logical(LKG)                , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!#if LK3_ENABLED
!    PURE module subroutine setCompactNew_D2_LK3(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_LK3
!#endif
!        use pm_kind, only: LKG => LK3
!        logical(LKG)                , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        logical(LKG)                , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!#if LK2_ENABLED
!    PURE module subroutine setCompactNew_D2_LK2(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_LK2
!#endif
!        use pm_kind, only: LKG => LK2
!        logical(LKG)                , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        logical(LKG)                , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!#if LK1_ENABLED
!    PURE module subroutine setCompactNew_D2_LK1(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_LK1
!#endif
!        use pm_kind, only: LKG => LK1
!        logical(LKG)                , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        logical(LKG)                , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if CK5_ENABLED
!    PURE module subroutine setCompactNew_D2_CK5(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_CK5
!#endif
!        use pm_kind, only: CKG => CK5
!        complex(CKG)                , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        complex(CKG)                , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!#if CK4_ENABLED
!    PURE module subroutine setCompactNew_D2_CK4(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_CK4
!#endif
!        use pm_kind, only: CKG => CK4
!        complex(CKG)                , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        complex(CKG)                , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!#if CK3_ENABLED
!    PURE module subroutine setCompactNew_D2_CK3(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_CK3
!#endif
!        use pm_kind, only: CKG => CK3
!        complex(CKG)                , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        complex(CKG)                , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!#if CK2_ENABLED
!    PURE module subroutine setCompactNew_D2_CK2(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_CK2
!#endif
!        use pm_kind, only: CKG => CK2
!        complex(CKG)                , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        complex(CKG)                , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!#if CK1_ENABLED
!    PURE module subroutine setCompactNew_D2_CK1(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_CK1
!#endif
!        use pm_kind, only: CKG => CK1
!        complex(CKG)                , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        complex(CKG)                , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module subroutine setCompactNew_D2_RK5(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)                   , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        real(RKG)                   , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE module subroutine setCompactNew_D2_RK4(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)                   , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        real(RKG)                   , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE module subroutine setCompactNew_D2_RK3(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)                   , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        real(RKG)                   , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE module subroutine setCompactNew_D2_RK2(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)                   , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        real(RKG)                   , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE module subroutine setCompactNew_D2_RK1(array, weight, compact, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setCompactNew_D2_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)                   , intent(in)    , contiguous    :: array(:,:)
!        integer(IK)                 , intent(out)   , allocatable   :: weight(:)
!        integer(IK)                 , intent(in)                    :: dim
!        real(RKG)                   , intent(out)   , allocatable   :: compact(:,:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCompact

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayCompact ! LCOV_EXCL_LINE