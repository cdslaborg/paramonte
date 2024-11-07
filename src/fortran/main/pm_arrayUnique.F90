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
!>  This module contains procedures and generic interfaces for finding unique values of an input array of various types.
!>
!>  \test
!>  [test_pm_arrayUnique](@ref test_pm_arrayUnique)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayUnique

    use pm_container, only: cvi_type
    use pm_kind, only: SK, IK, LK
    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_arrayUnique"

!>  \cond excluded
! \bug
! The following bypasses the bug reported below that creates a conflict between Intel and gfortran.
#if __INTEL_COMPILER
#define LEN_ARRAY :
#else
#define LEN_ARRAY len(array)
#endif
!>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` for each element of the input sequence whose value is unique among all sequence element values, otherwise return `.false.`.<br>
    !>
    !>  \details
    !>  The uniqueness of the values can be optionally determined by the user by specifying
    !>  the external input function that checks the equivalence of a pair of elements of the input `array`.<br>
    !>
    !>  \param[in]  array       :   The input `contiguous` array of shape `(:)` of either <br>
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
    !>                              whose elements uniqueness are to be tested.
    !>  \param      iseq        :   The `external` user-specified function that takes two input **scalar** arguments of the same type and kind as the input `array`.<br>
    !>                              It returns a scalar `logical` of default kind \LK that is `.true.` if the two input arguments are equivalent (e.g., equal)
    !>                              according to the user-defined criterion, otherwise, it is `.false.`.<br>
    !>                              The following illustrates the generic interface of `iseq`,
    !>                              \code{.F90}
    !>                                  function iseq(element1, element2) result(equivalent)
    !>                                      use pm_kind, only: LK
    !>                                      TYPE(KIND)  , intent(in)    :: element1, element2
    !>                                      logical(LK)                 :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `setA`, which can be one of the following,
    !>                              \code{.F90}
    !>                                  use pm_kind, only: SK, IK, CK, RK
    !>                                  character(*, SK), intent(in)    :: element1, element2 ! when `array` is a string vector.
    !>                                  character(1, SK), intent(in)    :: element1, element2 ! when `array` is a string scalar.
    !>                                  integer(IK)     , intent(in)    :: element1, element2
    !>                                  logical(LK)     , intent(in)    :: element1, element2
    !>                                  complex(CK)     , intent(in)    :: element1, element2
    !>                                  real(RK)        , intent(in)    :: element1, element2
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              This user-defined equivalence check is extremely useful where a user-defined equivalence test other than exact equality
    !>                              or identity is needed, for example, when the array elements should match only within a given threshold or,
    !>                              when the case-sensitivity in character comparisons do not matter.<br>
    !>                              In such cases, the user can define a custom equivalence criterion within the user-defined external function `iseq` to achieve the goal.<br>
    !>                              (**optional**, the default equivalence operator is `.eqv.` if the input `array` is `logical`, otherwise `==`.)
    !>
    !>  \return
    !>  `unique`                :   The output `logical` vector of default kind \LK of the same size as the length of the input sequence,
    !>                              each element of which is `.true.` if and only if the corresponding element of the input sequence has a unique value.<br>
    !>
    !>  \interface{isUnique}
    !>  \code{.F90}
    !>
    !>      use pm_arrayUnique, only: isUnique
    !>
    !>      unique(:) = isUnique(array)
    !>      unique(:) = isUnique(array, iseq)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \warning
    !>  The procedures under this generic interface are `impure` when the user-specified `external` procedure `iseq` is specified as input argument.
    !>
    !>  \warning
    !>  Note that in Fortran, trailing blanks are ignored in character comparison, that is, `"Fortran" == "Fortran "` yields `.true.`.
    !>
    !>  \see
    !>  [getUnique](@ref pm_arrayUnique::getUnique)<br>
    !>  [setUnique](@ref pm_arrayUnique::setUnique)<br>
    !>  [isUnique](@ref pm_arrayUnique::isUnique)<br>
    !>  [isUniqueAny](@ref pm_arrayUnique::isUniqueAny)<br>
    !>
    !>  \example{isUnique}
    !>  \include{lineno} example/pm_arrayUnique/isUnique/main.F90
    !>  \compilef{isUnique}
    !>  \output{isUnique}
    !>  \include{lineno} example/pm_arrayUnique/isUnique/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayUnique](@ref test_pm_arrayUnique)
    !>
    !>  \final{isUnique}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isUnique

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function isUniqueDefCom_D0_SK5(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: array
        logical(LK)                                             :: unique(len(array, IK))
    end function
#endif

#if SK4_ENABLED
    pure module function isUniqueDefCom_D0_SK4(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: array
        logical(LK)                                             :: unique(len(array, IK))
    end function
#endif

#if SK3_ENABLED
    pure module function isUniqueDefCom_D0_SK3(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: array
        logical(LK)                                             :: unique(len(array, IK))
    end function
#endif

#if SK2_ENABLED
    pure module function isUniqueDefCom_D0_SK2(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: array
        logical(LK)                                             :: unique(len(array, IK))
    end function
#endif

#if SK1_ENABLED
    pure module function isUniqueDefCom_D0_SK1(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: array
        logical(LK)                                             :: unique(len(array, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function isUniqueDefCom_D1_SK5(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if SK4_ENABLED
    pure module function isUniqueDefCom_D1_SK4(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if SK3_ENABLED
    pure module function isUniqueDefCom_D1_SK3(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if SK2_ENABLED
    pure module function isUniqueDefCom_D1_SK2(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if SK1_ENABLED
    pure module function isUniqueDefCom_D1_SK1(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module function isUniqueDefCom_D1_IK5(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if IK4_ENABLED
    pure module function isUniqueDefCom_D1_IK4(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if IK3_ENABLED
    pure module function isUniqueDefCom_D1_IK3(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if IK2_ENABLED
    pure module function isUniqueDefCom_D1_IK2(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if IK1_ENABLED
    pure module function isUniqueDefCom_D1_IK1(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module function isUniqueDefCom_D1_LK5(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if LK4_ENABLED
    pure module function isUniqueDefCom_D1_LK4(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if LK3_ENABLED
    pure module function isUniqueDefCom_D1_LK3(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if LK2_ENABLED
    pure module function isUniqueDefCom_D1_LK2(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if LK1_ENABLED
    pure module function isUniqueDefCom_D1_LK1(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function isUniqueDefCom_D1_CK5(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if CK4_ENABLED
    pure module function isUniqueDefCom_D1_CK4(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if CK3_ENABLED
    pure module function isUniqueDefCom_D1_CK3(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if CK2_ENABLED
    pure module function isUniqueDefCom_D1_CK2(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if CK1_ENABLED
    pure module function isUniqueDefCom_D1_CK1(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module function isUniqueDefCom_D1_RK5(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if RK4_ENABLED
    pure module function isUniqueDefCom_D1_RK4(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if RK3_ENABLED
    pure module function isUniqueDefCom_D1_RK3(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if RK2_ENABLED
    pure module function isUniqueDefCom_D1_RK2(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if RK1_ENABLED
    pure module function isUniqueDefCom_D1_RK1(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueDefCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function isUniqueCusCom_D0_SK5(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: array
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(len(array, IK))
    end function
#endif

#if SK4_ENABLED
    module function isUniqueCusCom_D0_SK4(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: array
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(len(array, IK))
    end function
#endif

#if SK3_ENABLED
    module function isUniqueCusCom_D0_SK3(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: array
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(len(array, IK))
    end function
#endif

#if SK2_ENABLED
    module function isUniqueCusCom_D0_SK2(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: array
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(len(array, IK))
    end function
#endif

#if SK1_ENABLED
    module function isUniqueCusCom_D0_SK1(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: array
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(len(array, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function isUniqueCusCom_D1_SK5(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if SK4_ENABLED
    module function isUniqueCusCom_D1_SK4(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if SK3_ENABLED
    module function isUniqueCusCom_D1_SK3(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if SK2_ENABLED
    module function isUniqueCusCom_D1_SK2(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if SK1_ENABLED
    module function isUniqueCusCom_D1_SK1(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function isUniqueCusCom_D1_IK5(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if IK4_ENABLED
    module function isUniqueCusCom_D1_IK4(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if IK3_ENABLED
    module function isUniqueCusCom_D1_IK3(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if IK2_ENABLED
    module function isUniqueCusCom_D1_IK2(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if IK1_ENABLED
    module function isUniqueCusCom_D1_IK1(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function isUniqueCusCom_D1_LK5(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if LK4_ENABLED
    module function isUniqueCusCom_D1_LK4(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if LK3_ENABLED
    module function isUniqueCusCom_D1_LK3(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if LK2_ENABLED
    module function isUniqueCusCom_D1_LK2(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if LK1_ENABLED
    module function isUniqueCusCom_D1_LK1(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function isUniqueCusCom_D1_CK5(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if CK4_ENABLED
    module function isUniqueCusCom_D1_CK4(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if CK3_ENABLED
    module function isUniqueCusCom_D1_CK3(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if CK2_ENABLED
    module function isUniqueCusCom_D1_CK2(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if CK1_ENABLED
    module function isUniqueCusCom_D1_CK1(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function isUniqueCusCom_D1_RK5(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if RK4_ENABLED
    module function isUniqueCusCom_D1_RK4(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if RK3_ENABLED
    module function isUniqueCusCom_D1_RK3(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if RK2_ENABLED
    module function isUniqueCusCom_D1_RK2(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

#if RK1_ENABLED
    module function isUniqueCusCom_D1_RK1(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueCusCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: unique(size(array, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface isUnique

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if and only if all elements of the input sequence values are unique.<br>
    !>
    !>  \details
    !>  The uniqueness of the values can be optionally determined by the user by specifying
    !>  the external input function that checks the equivalence of a pair of elements of the input `array`.<br>
    !>
    !>  \param[in]  array       :   The input `contiguous` array of shape `(:)` of either <br>
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
    !>                              whose elements uniqueness are to be tested.
    !>  \param      iseq        :   The `external` user-specified function that takes two input **scalar** arguments of the same type and kind as the input `array`.<br>
    !>                              It returns a scalar `logical` of default kind \LK that is `.true.` if the two input arguments are equivalent (e.g., equal)
    !>                              according to the user-defined criterion, otherwise, it is `.false.`.<br>
    !>                              The following illustrates the generic interface of `iseq`,
    !>                              \code{.F90}
    !>                                  function iseq(element1, element2) result(equivalent)
    !>                                      use pm_kind, only: LK
    !>                                      TYPE(KIND)  , intent(in)    :: element1, element2
    !>                                      logical(LK)                 :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `setA`, which can be one of the following,
    !>                              \code{.F90}
    !>                                  use pm_kind, only: SK, IK, CK, RK
    !>                                  character(*, SK), intent(in)    :: element1, element2 ! when `array` is a string vector.
    !>                                  character(1, SK), intent(in)    :: element1, element2 ! when `array` is a string scalar.
    !>                                  integer(IK)     , intent(in)    :: element1, element2
    !>                                  logical(LK)     , intent(in)    :: element1, element2
    !>                                  complex(CK)     , intent(in)    :: element1, element2
    !>                                  real(RK)        , intent(in)    :: element1, element2
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              This user-defined equivalence check is extremely useful where a user-defined equivalence test other than exact equality
    !>                              or identity is needed, for example, when the array elements should match only within a given threshold or,
    !>                              when the case-sensitivity in character comparisons do not matter.<br>
    !>                              In such cases, the user can define a custom equivalence criterion within the user-defined external function `iseq` to achieve the goal.<br>
    !>                              (**optional**, the default equivalence operator is `.eqv.` if the input `array` is `logical`, otherwise `==`.)
    !>
    !>  \return
    !>  `uniqueAll`             :   The output scalar `logical` of default kind \LK that is `.true.` if and only if all elements of the input sequence are unique.<br>
    !>
    !>  \interface{isUniqueAll}
    !>  \code{.F90}
    !>
    !>      use pm_arrayUnique, only: isUniqueAll
    !>
    !>      uniqueAll = isUniqueAll(array)
    !>      uniqueAll = isUniqueAll(array, iseq)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \warning
    !>  The procedures under this generic interface are `impure` when the user-specified `external` procedure `iseq` is specified as input argument.
    !>
    !>  \warning
    !>  Note that in Fortran, trailing blanks are ignored in character comparison, that is, `"Fortran" == "Fortran "` yields `.true.`.
    !>
    !>  \see
    !>  [getUnique](@ref pm_arrayUnique::getUnique)<br>
    !>  [setUnique](@ref pm_arrayUnique::setUnique)<br>
    !>  [isUniqueAll](@ref pm_arrayUnique::isUniqueAll)<br>
    !>  [isUniqueAny](@ref pm_arrayUnique::isUniqueAny)<br>
    !>
    !>  \example{isUniqueAll}
    !>  \include{lineno} example/pm_arrayUnique/isUniqueAll/main.F90
    !>  \compilef{isUniqueAll}
    !>  \output{isUniqueAll}
    !>  \include{lineno} example/pm_arrayUnique/isUniqueAll/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayUnique](@ref test_pm_arrayUnique)
    !>
    !>  \final{isUniqueAll}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isUniqueAll

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function isUniqueAllDefCom_D0_SK5(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: array
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if SK4_ENABLED
    PURE module function isUniqueAllDefCom_D0_SK4(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: array
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if SK3_ENABLED
    PURE module function isUniqueAllDefCom_D0_SK3(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: array
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if SK2_ENABLED
    PURE module function isUniqueAllDefCom_D0_SK2(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: array
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if SK1_ENABLED
    PURE module function isUniqueAllDefCom_D0_SK1(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: array
        logical(LK)                                             :: uniqueAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function isUniqueAllDefCom_D1_SK5(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if SK4_ENABLED
    PURE module function isUniqueAllDefCom_D1_SK4(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if SK3_ENABLED
    PURE module function isUniqueAllDefCom_D1_SK3(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if SK2_ENABLED
    PURE module function isUniqueAllDefCom_D1_SK2(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if SK1_ENABLED
    PURE module function isUniqueAllDefCom_D1_SK1(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function isUniqueAllDefCom_D1_IK5(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if IK4_ENABLED
    PURE module function isUniqueAllDefCom_D1_IK4(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if IK3_ENABLED
    PURE module function isUniqueAllDefCom_D1_IK3(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if IK2_ENABLED
    PURE module function isUniqueAllDefCom_D1_IK2(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if IK1_ENABLED
    PURE module function isUniqueAllDefCom_D1_IK1(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function isUniqueAllDefCom_D1_LK5(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if LK4_ENABLED
    PURE module function isUniqueAllDefCom_D1_LK4(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if LK3_ENABLED
    PURE module function isUniqueAllDefCom_D1_LK3(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if LK2_ENABLED
    PURE module function isUniqueAllDefCom_D1_LK2(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if LK1_ENABLED
    PURE module function isUniqueAllDefCom_D1_LK1(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function isUniqueAllDefCom_D1_CK5(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if CK4_ENABLED
    PURE module function isUniqueAllDefCom_D1_CK4(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if CK3_ENABLED
    PURE module function isUniqueAllDefCom_D1_CK3(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if CK2_ENABLED
    PURE module function isUniqueAllDefCom_D1_CK2(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if CK1_ENABLED
    PURE module function isUniqueAllDefCom_D1_CK1(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function isUniqueAllDefCom_D1_RK5(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if RK4_ENABLED
    PURE module function isUniqueAllDefCom_D1_RK4(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if RK3_ENABLED
    PURE module function isUniqueAllDefCom_D1_RK3(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if RK2_ENABLED
    PURE module function isUniqueAllDefCom_D1_RK2(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if RK1_ENABLED
    PURE module function isUniqueAllDefCom_D1_RK1(array) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllDefCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function isUniqueAllCusCom_D0_SK5(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: array
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if SK4_ENABLED
    module function isUniqueAllCusCom_D0_SK4(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: array
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if SK3_ENABLED
    module function isUniqueAllCusCom_D0_SK3(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: array
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if SK2_ENABLED
    module function isUniqueAllCusCom_D0_SK2(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: array
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if SK1_ENABLED
    module function isUniqueAllCusCom_D0_SK1(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: array
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function isUniqueAllCusCom_D1_SK5(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if SK4_ENABLED
    module function isUniqueAllCusCom_D1_SK4(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if SK3_ENABLED
    module function isUniqueAllCusCom_D1_SK3(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if SK2_ENABLED
    module function isUniqueAllCusCom_D1_SK2(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if SK1_ENABLED
    module function isUniqueAllCusCom_D1_SK1(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function isUniqueAllCusCom_D1_IK5(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if IK4_ENABLED
    module function isUniqueAllCusCom_D1_IK4(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if IK3_ENABLED
    module function isUniqueAllCusCom_D1_IK3(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if IK2_ENABLED
    module function isUniqueAllCusCom_D1_IK2(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if IK1_ENABLED
    module function isUniqueAllCusCom_D1_IK1(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function isUniqueAllCusCom_D1_LK5(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if LK4_ENABLED
    module function isUniqueAllCusCom_D1_LK4(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if LK3_ENABLED
    module function isUniqueAllCusCom_D1_LK3(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if LK2_ENABLED
    module function isUniqueAllCusCom_D1_LK2(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if LK1_ENABLED
    module function isUniqueAllCusCom_D1_LK1(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function isUniqueAllCusCom_D1_CK5(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if CK4_ENABLED
    module function isUniqueAllCusCom_D1_CK4(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if CK3_ENABLED
    module function isUniqueAllCusCom_D1_CK3(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if CK2_ENABLED
    module function isUniqueAllCusCom_D1_CK2(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if CK1_ENABLED
    module function isUniqueAllCusCom_D1_CK1(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function isUniqueAllCusCom_D1_RK5(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if RK4_ENABLED
    module function isUniqueAllCusCom_D1_RK4(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if RK3_ENABLED
    module function isUniqueAllCusCom_D1_RK3(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if RK2_ENABLED
    module function isUniqueAllCusCom_D1_RK2(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

#if RK1_ENABLED
    module function isUniqueAllCusCom_D1_RK1(array, iseq) result(uniqueAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAllCusCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface isUniqueAll

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if and only if at least one element of the input sequence is unique among others.<br>
    !>
    !>  \details
    !>  The uniqueness of the values can be optionally determined by the user by specifying
    !>  the external input function that checks the equivalence of a pair of elements of the input `array`.<br>
    !>
    !>  \param[in]  array       :   The input `contiguous` array of shape `(:)` of either <br>
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
    !>                              whose elements uniqueness are to be tested.
    !>  \param      iseq        :   The `external` user-specified function that takes two input **scalar** arguments of the same type and kind as the input `array`.<br>
    !>                              It returns a scalar `logical` of default kind \LK that is `.true.` if the two input arguments are equivalent (e.g., equal)
    !>                              according to the user-defined criterion, otherwise, it is `.false.`.<br>
    !>                              The following illustrates the generic interface of `iseq`,
    !>                              \code{.F90}
    !>                                  function iseq(element1, element2) result(equivalent)
    !>                                      use pm_kind, only: LK
    !>                                      TYPE(KIND)  , intent(in)    :: element1, element2
    !>                                      logical(LK)                 :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `setA`, which can be one of the following,
    !>                              \code{.F90}
    !>                                  use pm_kind, only: SK, IK, CK, RK
    !>                                  character(*, SK), intent(in)    :: element1, element2 ! when `array` is a string vector.
    !>                                  character(1, SK), intent(in)    :: element1, element2 ! when `array` is a string scalar.
    !>                                  integer(IK)     , intent(in)    :: element1, element2
    !>                                  logical(LK)     , intent(in)    :: element1, element2
    !>                                  complex(CK)     , intent(in)    :: element1, element2
    !>                                  real(RK)        , intent(in)    :: element1, element2
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              This user-defined equivalence check is extremely useful where a user-defined equivalence test other than exact equality
    !>                              or identity is needed, for example, when the array elements should match only within a given threshold or,
    !>                              when the case-sensitivity in character comparisons do not matter.<br>
    !>                              In such cases, the user can define a custom equivalence criterion within the user-defined external function `iseq` to achieve the goal.<br>
    !>                              (**optional**, the default equivalence operator is `.eqv.` if the input `array` is `logical`, otherwise `==`.)
    !>
    !>  \return
    !>  `uniqueAny`             :   The output scalar `logical` of default kind \LK that is `.true.` if and only if any elements of the input sequence are unique.<br>
    !>
    !>  \interface{isUniqueAny}
    !>  \code{.F90}
    !>
    !>      use pm_arrayUnique, only: isUniqueAny
    !>
    !>      uniqueAny = isUniqueAny(array)
    !>      uniqueAny = isUniqueAny(array, iseq)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \warning
    !>  The procedures under this generic interface are `impure` when the user-specified `external` procedure `iseq` is specified as input argument.
    !>
    !>  \warning
    !>  Note that in Fortran, trailing blanks are ignored in character comparison, that is, `"Fortran" == "Fortran "` yields `.true.`.
    !>
    !>  \see
    !>  [getUnique](@ref pm_arrayUnique::getUnique)<br>
    !>  [setUnique](@ref pm_arrayUnique::setUnique)<br>
    !>  [isUniqueAny](@ref pm_arrayUnique::isUniqueAny)<br>
    !>
    !>  \example{isUniqueAny}
    !>  \include{lineno} example/pm_arrayUnique/isUniqueAny/main.F90
    !>  \compilef{isUniqueAny}
    !>  \output{isUniqueAny}
    !>  \include{lineno} example/pm_arrayUnique/isUniqueAny/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayUnique](@ref test_pm_arrayUnique)
    !>
    !>  \final{isUniqueAny}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isUniqueAny

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function isUniqueAnyDefCom_D0_SK5(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: array
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if SK4_ENABLED
    PURE module function isUniqueAnyDefCom_D0_SK4(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: array
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if SK3_ENABLED
    PURE module function isUniqueAnyDefCom_D0_SK3(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: array
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if SK2_ENABLED
    PURE module function isUniqueAnyDefCom_D0_SK2(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: array
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if SK1_ENABLED
    PURE module function isUniqueAnyDefCom_D0_SK1(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: array
        logical(LK)                                             :: uniqueAny
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function isUniqueAnyDefCom_D1_SK5(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if SK4_ENABLED
    PURE module function isUniqueAnyDefCom_D1_SK4(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if SK3_ENABLED
    PURE module function isUniqueAnyDefCom_D1_SK3(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if SK2_ENABLED
    PURE module function isUniqueAnyDefCom_D1_SK2(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if SK1_ENABLED
    PURE module function isUniqueAnyDefCom_D1_SK1(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function isUniqueAnyDefCom_D1_IK5(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if IK4_ENABLED
    PURE module function isUniqueAnyDefCom_D1_IK4(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if IK3_ENABLED
    PURE module function isUniqueAnyDefCom_D1_IK3(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if IK2_ENABLED
    PURE module function isUniqueAnyDefCom_D1_IK2(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if IK1_ENABLED
    PURE module function isUniqueAnyDefCom_D1_IK1(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function isUniqueAnyDefCom_D1_LK5(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if LK4_ENABLED
    PURE module function isUniqueAnyDefCom_D1_LK4(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if LK3_ENABLED
    PURE module function isUniqueAnyDefCom_D1_LK3(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if LK2_ENABLED
    PURE module function isUniqueAnyDefCom_D1_LK2(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if LK1_ENABLED
    PURE module function isUniqueAnyDefCom_D1_LK1(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function isUniqueAnyDefCom_D1_CK5(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if CK4_ENABLED
    PURE module function isUniqueAnyDefCom_D1_CK4(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if CK3_ENABLED
    PURE module function isUniqueAnyDefCom_D1_CK3(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if CK2_ENABLED
    PURE module function isUniqueAnyDefCom_D1_CK2(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if CK1_ENABLED
    PURE module function isUniqueAnyDefCom_D1_CK1(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function isUniqueAnyDefCom_D1_RK5(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if RK4_ENABLED
    PURE module function isUniqueAnyDefCom_D1_RK4(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if RK3_ENABLED
    PURE module function isUniqueAnyDefCom_D1_RK3(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if RK2_ENABLED
    PURE module function isUniqueAnyDefCom_D1_RK2(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if RK1_ENABLED
    PURE module function isUniqueAnyDefCom_D1_RK1(array) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyDefCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in), contiguous    :: array(:)
        logical(LK)                                             :: uniqueAny
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function isUniqueAnyCusCom_D0_SK5(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: array
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if SK4_ENABLED
    module function isUniqueAnyCusCom_D0_SK4(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: array
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if SK3_ENABLED
    module function isUniqueAnyCusCom_D0_SK3(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: array
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if SK2_ENABLED
    module function isUniqueAnyCusCom_D0_SK2(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: array
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if SK1_ENABLED
    module function isUniqueAnyCusCom_D0_SK1(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: array
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function isUniqueAnyCusCom_D1_SK5(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if SK4_ENABLED
    module function isUniqueAnyCusCom_D1_SK4(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if SK3_ENABLED
    module function isUniqueAnyCusCom_D1_SK3(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if SK2_ENABLED
    module function isUniqueAnyCusCom_D1_SK2(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if SK1_ENABLED
    module function isUniqueAnyCusCom_D1_SK1(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function isUniqueAnyCusCom_D1_IK5(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if IK4_ENABLED
    module function isUniqueAnyCusCom_D1_IK4(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if IK3_ENABLED
    module function isUniqueAnyCusCom_D1_IK3(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if IK2_ENABLED
    module function isUniqueAnyCusCom_D1_IK2(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if IK1_ENABLED
    module function isUniqueAnyCusCom_D1_IK1(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function isUniqueAnyCusCom_D1_LK5(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if LK4_ENABLED
    module function isUniqueAnyCusCom_D1_LK4(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if LK3_ENABLED
    module function isUniqueAnyCusCom_D1_LK3(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if LK2_ENABLED
    module function isUniqueAnyCusCom_D1_LK2(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if LK1_ENABLED
    module function isUniqueAnyCusCom_D1_LK1(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function isUniqueAnyCusCom_D1_CK5(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if CK4_ENABLED
    module function isUniqueAnyCusCom_D1_CK4(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if CK3_ENABLED
    module function isUniqueAnyCusCom_D1_CK3(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if CK2_ENABLED
    module function isUniqueAnyCusCom_D1_CK2(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if CK1_ENABLED
    module function isUniqueAnyCusCom_D1_CK1(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function isUniqueAnyCusCom_D1_RK5(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if RK4_ENABLED
    module function isUniqueAnyCusCom_D1_RK4(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if RK3_ENABLED
    module function isUniqueAnyCusCom_D1_RK3(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if RK2_ENABLED
    module function isUniqueAnyCusCom_D1_RK2(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

#if RK1_ENABLED
    module function isUniqueAnyCusCom_D1_RK1(array, iseq) result(uniqueAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isUniqueAnyCusCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LK)                                             :: uniqueAny
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface isUniqueAny

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a vector of unique values in the input array.<br>
    !>
    !>  \details
    !>  The uniqueness of the values can be optionally determined by the user by specifying
    !>  the external input function that checks the equivalence of a pair of elements of the input `array`.
    !>
    !>  \param[in]  array       :   The input `contiguous` array of shape `(:)` of either <br>
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
    !>                              whose unique elements are to be returned.
    !>  \param      iseq        :   The `external` user-specified function that takes two input **scalar** arguments of the same type and kind as the input `array`.<br>
    !>                              It returns a scalar `logical` of default kind \LK that is `.true.` if the two input arguments are equivalent (e.g., equal)
    !>                              according to the user-defined criterion, otherwise, it is `.false.`.<br>
    !>                              The following illustrates the generic interface of `iseq`,
    !>                              \code{.F90}
    !>                                  function iseq(element1, element2) result(equivalent)
    !>                                      use pm_kind, only: LK
    !>                                      TYPE(KIND)  , intent(in)    :: element1, element2
    !>                                      logical(LK)                 :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `setA`, which can be one of the following,
    !>                              \code{.F90}
    !>                                  use pm_kind, only: SK, IK, CK, RK
    !>                                  character(*, SK), intent(in)    :: element1, element2 ! when `array` is a string vector.
    !>                                  character(1, SK), intent(in)    :: element1, element2 ! when `array` is a string scalar.
    !>                                  integer(IK)     , intent(in)    :: element1, element2
    !>                                  logical(LK)     , intent(in)    :: element1, element2
    !>                                  complex(CK)     , intent(in)    :: element1, element2
    !>                                  real(RK)        , intent(in)    :: element1, element2
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              This user-defined equivalence check is extremely useful where a user-defined equivalence test other than exact equality
    !>                              or identity is needed, for example, when the array elements should match only within a given threshold or,
    !>                              when the case-sensitivity in character comparisons do not matter.<br>
    !>                              In such cases, the user can define a custom equivalence criterion within the user-defined external function `iseq` to achieve the goal.<br>
    !>                              (**optional**, the default equivalence operator is `.eqv.` if the input `array` is `logical`, otherwise `==`.)
    !>
    !>  \return
    !>  `unique`                :   The output `allocatable` array of the same type and kind as the input `array` containing all of its unique elements.
    !>
    !>  \interface{getUnique}
    !>  \code{.F90}
    !>
    !>      use pm_arrayUnique, only: getUnique
    !>
    !>      unique = getUnique(array)
    !>      unique = getUnique(array, iseq)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \warning
    !>  The procedures under this generic interface are `impure` when the user-specified `external` procedure `iseq` is specified as input argument.
    !>
    !>  \warning
    !>  Note that in Fortran, trailing blanks are ignored in character comparison, that is, `"Fortran" == "Fortran "` yields `.true.`.
    !>
    !>  \remark
    !>  The functions under this generic interface are slightly slower than the [setUnique](@ref pm_arrayUnique::setUnique) subroutine implementations.
    !>
    !>  \note
    !>  If needed, use [setSorted](@ref pm_arraySort::setSorted) to sort the output unique values (in ascending order).
    !>
    !>  \see
    !>  [getUnique](@ref pm_arrayUnique::getUnique)<br>
    !>  [setUnique](@ref pm_arrayUnique::setUnique)<br>
    !>  [isUniqueAll](@ref pm_arrayUnique::isUniqueAll)<br>
    !>  [isUniqueAny](@ref pm_arrayUnique::isUniqueAny)<br>
    !>
    !>  \example{getUnique}
    !>  \include{lineno} example/pm_arrayUnique/getUnique/main.F90
    !>  \compilef{getUnique}
    !>  \output{getUnique}
    !>  \include{lineno} example/pm_arrayUnique/getUnique/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayUnique](@ref test_pm_arrayUnique)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.2.0}, \gfortran{10-11}
    !>  \desc
    !>  The \ifort{2021.2.0} has a bug for the following interface definition,<br>
    !>  \code{.F90}
    !>      character(len(array),SK)   , allocatable   :: arrayNew(:)
    !>  \endcode
    !>  leading to an **internal compiler error**.<br>
    !>  \remedy
    !>  For now, the remedy seems to be to redefine the interface as,
    !>  \code{.F90}
    !>      character(:, SK), allocatable   :: arrayNew(:)
    !>  \endcode
    !>  and changing the allocation method accordingly in the implementation to,
    !>  \code{.F90}
    !>      allocate(character(len(array, kind = IK)) :: arrayNew(lenArrayNew))
    !>  \endcode
    !>  However, this introduces `internal compiler error: Segmentation fault` with \gfortran.
    !>  Here is a code snippet to regenerate the bug in \ifort (uncomment the commented line to reproduce the gfortran bug),
    !>  \code{.F90}
    !>
    !>      module pm_explicitLenResult
    !>          implicit none
    !>          interface
    !>              pure module function bug(array) result(arrayNew)
    !>                  character(*, SK), intent(in), contiguous    :: array(:)
    !>                  character(len(array),SK)    , allocatable   :: arrayNew(:) ! catastrophic internal error with ifort 2021.2. Fine with gfortran 10.3
    !>                 !character(:, SK)            , allocatable   :: arrayNew(:) ! catastrophic internal error with gfortran 10.3. Fine with ifort 2021.2
    !>              end function
    !>          end interface
    !>      end module pm_explicitLenResult
    !>
    !>      submodule (pm_explicitLenResult) routines
    !>          implicit none
    !>      contains
    !>          module procedure bug
    !>             allocate(arrayNew, source = array)
    !>          end procedure
    !>      end submodule routines
    !>
    !>      program main
    !>          use pm_explicitLenResult, only: bug
    !>          character(2) :: array(3) = ["AA", "BB", "CC"]
    !>          character(2), allocatable :: arrayNew(:)
    !>          arrayNew = bug(array)
    !>      end program main
    !>
    !>  \endcode
    !>  \remedy
    !>  It turns out that both \gfortran and \ifort do not tolerate the separation of interface from implementation in the above code snippet.<br>
    !>  If one duplicates the interface in the implementation submodule, then both compilers compile and run the code with no errors.<br>
    !>  This is the remedy that is currently used in this [getRemoved](@ref pm_arrayRemove::getRemoved) generic interface
    !>  (interface duplication where the bug exists).<br>
    !>  Here is a workaround example for the bug in the above code snippet,
    !>  \code{.F90}
    !>
    !>      module pm_explicitLenResult
    !>          implicit none
    !>          interface
    !>              pure module function bug(array) result(arrayNew)
    !>                  character(*, SK), intent(in), contiguous    :: array(:)
    !>                  character(len(array),SK)   , allocatable   :: arrayNew(:) ! catastrophic internal error with ifort 2021.2. Fine with gfortran 10.3
    !>              end function
    !>          end interface
    !>      end module pm_explicitLenResult
    !>
    !>      submodule (pm_explicitLenResult) routines
    !>          implicit none
    !>      contains
    !>          module procedure bug
    !>             allocate(arrayNew, source = array)
    !>          end procedure
    !>      end submodule routines
    !>
    !>      program main
    !>          use pm_explicitLenResult, only: bug
    !>          character(2) :: array(3) = ["AA", "BB", "CC"]
    !>          character(2), allocatable :: arrayNew(:)
    !>          arrayNew = bug(array)
    !>      end program main
    !>
    !>  \endcode
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to 2D input objects.
    !>
    !>  \todo
    !>  \pvhigh The internal compiler error with `ifort` and `gfortran` has to be fixed in the future versions.<br>
    !>
    !>  \final{getUnique}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getUnique

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getUniArbDefCom_D0_SK5(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: array
        character(:,SKG)            , allocatable               :: unique
    end function
#endif

#if SK4_ENABLED
    PURE module function getUniArbDefCom_D0_SK4(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: array
        character(:,SKG)            , allocatable               :: unique
    end function
#endif

#if SK3_ENABLED
    PURE module function getUniArbDefCom_D0_SK3(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: array
        character(:,SKG)            , allocatable               :: unique
    end function
#endif

#if SK2_ENABLED
    PURE module function getUniArbDefCom_D0_SK2(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: array
        character(:,SKG)            , allocatable               :: unique
    end function
#endif

#if SK1_ENABLED
    PURE module function getUniArbDefCom_D0_SK1(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: array
        character(:,SKG)            , allocatable               :: unique
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getUniArbDefCom_D1_SK5(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        character(LEN_ARRAY,SKG)    , allocatable               :: unique(:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getUniArbDefCom_D1_SK4(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        character(LEN_ARRAY,SKG)    , allocatable               :: unique(:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getUniArbDefCom_D1_SK3(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        character(LEN_ARRAY,SKG)    , allocatable               :: unique(:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getUniArbDefCom_D1_SK2(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        character(LEN_ARRAY,SKG)    , allocatable               :: unique(:)
    end function
#endif

#if SK1_ENABLED
    PURE module function getUniArbDefCom_D1_SK1(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        character(LEN_ARRAY,SKG)    , allocatable               :: unique(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getUniArbDefCom_D1_IK5(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in), contiguous    :: array(:)
        integer(IKG)                , allocatable               :: unique(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getUniArbDefCom_D1_IK4(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in), contiguous    :: array(:)
        integer(IKG)                , allocatable               :: unique(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getUniArbDefCom_D1_IK3(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in), contiguous    :: array(:)
        integer(IKG)                , allocatable               :: unique(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getUniArbDefCom_D1_IK2(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in), contiguous    :: array(:)
        integer(IKG)                , allocatable               :: unique(:)
    end function
#endif

#if IK1_ENABLED
    PURE module function getUniArbDefCom_D1_IK1(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in), contiguous    :: array(:)
        integer(IKG)                , allocatable               :: unique(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getUniArbDefCom_D1_LK5(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in), contiguous    :: array(:)
        logical(LKG)                , allocatable               :: unique(:)
    end function
#endif

#if LK4_ENABLED
    PURE module function getUniArbDefCom_D1_LK4(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in), contiguous    :: array(:)
        logical(LKG)                , allocatable               :: unique(:)
    end function
#endif

#if LK3_ENABLED
    PURE module function getUniArbDefCom_D1_LK3(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in), contiguous    :: array(:)
        logical(LKG)                , allocatable               :: unique(:)
    end function
#endif

#if LK2_ENABLED
    PURE module function getUniArbDefCom_D1_LK2(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in), contiguous    :: array(:)
        logical(LKG)                , allocatable               :: unique(:)
    end function
#endif

#if LK1_ENABLED
    PURE module function getUniArbDefCom_D1_LK1(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in), contiguous    :: array(:)
        logical(LKG)                , allocatable               :: unique(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getUniArbDefCom_D1_CK5(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in), contiguous    :: array(:)
        complex(CKG)                , allocatable               :: unique(:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getUniArbDefCom_D1_CK4(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in), contiguous    :: array(:)
        complex(CKG)                , allocatable               :: unique(:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getUniArbDefCom_D1_CK3(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in), contiguous    :: array(:)
        complex(CKG)                , allocatable               :: unique(:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getUniArbDefCom_D1_CK2(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in), contiguous    :: array(:)
        complex(CKG)                , allocatable               :: unique(:)
    end function
#endif

#if CK1_ENABLED
    PURE module function getUniArbDefCom_D1_CK1(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in), contiguous    :: array(:)
        complex(CKG)                , allocatable               :: unique(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getUniArbDefCom_D1_RK5(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in), contiguous    :: array(:)
        real(RKG)                   , allocatable               :: unique(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getUniArbDefCom_D1_RK4(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in), contiguous    :: array(:)
        real(RKG)                   , allocatable               :: unique(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getUniArbDefCom_D1_RK3(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in), contiguous    :: array(:)
        real(RKG)                   , allocatable               :: unique(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getUniArbDefCom_D1_RK2(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in), contiguous    :: array(:)
        real(RKG)                   , allocatable               :: unique(:)
    end function
#endif

#if RK1_ENABLED
    PURE module function getUniArbDefCom_D1_RK1(array) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbDefCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in), contiguous    :: array(:)
        real(RKG)                   , allocatable               :: unique(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getUniArbCusCom_D0_SK5(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: array
        procedure(logical(LK))                                  :: iseq
        character(:,SKG)            , allocatable               :: unique
    end function
#endif

#if SK4_ENABLED
    module function getUniArbCusCom_D0_SK4(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: array
        procedure(logical(LK))                                  :: iseq
        character(:,SKG)            , allocatable               :: unique
    end function
#endif

#if SK3_ENABLED
    module function getUniArbCusCom_D0_SK3(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: array
        procedure(logical(LK))                                  :: iseq
        character(:,SKG)            , allocatable               :: unique
    end function
#endif

#if SK2_ENABLED
    module function getUniArbCusCom_D0_SK2(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: array
        procedure(logical(LK))                                  :: iseq
        character(:,SKG)            , allocatable               :: unique
    end function
#endif

#if SK1_ENABLED
    module function getUniArbCusCom_D0_SK1(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: array
        procedure(logical(LK))                                  :: iseq
        character(:,SKG)            , allocatable               :: unique
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getUniArbCusCom_D1_SK5(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)    , allocatable               :: unique(:)
    end function
#endif

#if SK4_ENABLED
    module function getUniArbCusCom_D1_SK4(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)    , allocatable               :: unique(:)
    end function
#endif

#if SK3_ENABLED
    module function getUniArbCusCom_D1_SK3(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)    , allocatable               :: unique(:)
    end function
#endif

#if SK2_ENABLED
    module function getUniArbCusCom_D1_SK2(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)    , allocatable               :: unique(:)
    end function
#endif

#if SK1_ENABLED
    module function getUniArbCusCom_D1_SK1(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)    , allocatable               :: unique(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getUniArbCusCom_D1_IK5(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                , allocatable               :: unique(:)
    end function
#endif

#if IK4_ENABLED
    module function getUniArbCusCom_D1_IK4(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                , allocatable               :: unique(:)
    end function
#endif

#if IK3_ENABLED
    module function getUniArbCusCom_D1_IK3(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                , allocatable               :: unique(:)
    end function
#endif

#if IK2_ENABLED
    module function getUniArbCusCom_D1_IK2(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                , allocatable               :: unique(:)
    end function
#endif

#if IK1_ENABLED
    module function getUniArbCusCom_D1_IK1(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                , allocatable               :: unique(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getUniArbCusCom_D1_LK5(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                , allocatable               :: unique(:)
    end function
#endif

#if LK4_ENABLED
    module function getUniArbCusCom_D1_LK4(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                , allocatable               :: unique(:)
    end function
#endif

#if LK3_ENABLED
    module function getUniArbCusCom_D1_LK3(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                , allocatable               :: unique(:)
    end function
#endif

#if LK2_ENABLED
    module function getUniArbCusCom_D1_LK2(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                , allocatable               :: unique(:)
    end function
#endif

#if LK1_ENABLED
    module function getUniArbCusCom_D1_LK1(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                , allocatable               :: unique(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getUniArbCusCom_D1_CK5(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                , allocatable               :: unique(:)
    end function
#endif

#if CK4_ENABLED
    module function getUniArbCusCom_D1_CK4(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                , allocatable               :: unique(:)
    end function
#endif

#if CK3_ENABLED
    module function getUniArbCusCom_D1_CK3(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                , allocatable               :: unique(:)
    end function
#endif

#if CK2_ENABLED
    module function getUniArbCusCom_D1_CK2(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                , allocatable               :: unique(:)
    end function
#endif

#if CK1_ENABLED
    module function getUniArbCusCom_D1_CK1(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                , allocatable               :: unique(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getUniArbCusCom_D1_RK5(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        real(RKG)                   , allocatable               :: unique(:)
    end function
#endif

#if RK4_ENABLED
    module function getUniArbCusCom_D1_RK4(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        real(RKG)                   , allocatable               :: unique(:)
    end function
#endif

#if RK3_ENABLED
    module function getUniArbCusCom_D1_RK3(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        real(RKG)                   , allocatable               :: unique(:)
    end function
#endif

#if RK2_ENABLED
    module function getUniArbCusCom_D1_RK2(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        real(RKG)                   , allocatable               :: unique(:)
    end function
#endif

#if RK1_ENABLED
    module function getUniArbCusCom_D1_RK1(array, iseq) result(unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniArbCusCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in), contiguous    :: array(:)
        procedure(logical(LK))                                  :: iseq
        real(RKG)                   , allocatable               :: unique(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getUnique

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a vector of unique values in the input array **in place of** the array itself.<br>
    !>
    !>  \details
    !>  The uniqueness of the values can be optionally determined by the user by specifying
    !>  the external input function that checks the equivalence of a pair of elements of the input `array`.<br>
    !>  Furthermore, optionally return the count of each of the unique values in the input `array`,
    !>  as well as the indices of occurrences of each of the unique values.<br>
    !>
    !>  \param[in]  array       :   The input `contiguous` array of shape `(:)` of either <br>
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
    !>                              whose unique elements are to be returned.
    !>  \param[out] unique      :   The output array of the same type and kind as the input `array` containing all of its unique elements.<br>
    !>                              <ul>
    !>                                  <li>    If the output argument `lenUnique` is present, then `unique` must be a `contiguous` array (or scalar string) whose length equals the length of the input `array`.<br>
    !>                                  <li>    If the output argument `lenUnique` is missing, then `unique` must be an `allocatable` array (or scalar string).<br>
    !>                              </ul>
    !>  \param[out] lenUnique   :   The output scalar of type `integer` of default kind  of the same type and kind as the input `array` containing all of its unique elements.<br>
    !>                              (**optional**, if present, the specified arrays `unique`, `count`, `index` must be preallocated to the length of `array`.)<br>
    !>  \param[out] count       :   The output array of type `integer` of default kind \IK containing the corresponding counts of each of
    !>                              the uniques elements of `array` in the output `unique` array.<br>
    !>                              <ul>
    !>                                  <li>    If the output argument `lenUnique` is present, then `count` must be a `contiguous` array (or scalar string) whose length equals the length of the input `array`.<br>
    !>                                  <li>    If the output argument `lenUnique` is missing, then `count` must be an `allocatable` array (or scalar string).<br>
    !>                              </ul>
    !>  \param[out] index       :   The output jagged-array of type [cvi_type](@ref pm_container::cvi_type) of default kind \IK each
    !>                              element of which contains the vector of indices of the positions of occurrences of the uniques values in `array`.<br>
    !>                              <ul>
    !>                                  <li>    If the output argument `lenUnique` is present, then `index` must be a `contiguous` array (or scalar string) whose length equals the length of the input `array`.<br>
    !>                                  <li>    If the output argument `lenUnique` is missing, then `index` must be an `allocatable` array (or scalar string).<br>
    !>                              </ul>
    !>                              (**optional**, it can be missing.)
    !>  \param[in]  order       :   The input `integer` of default kind \IK.<br>
    !>                              If  `0`, the output array `count` will be not be sorted.<br>
    !>                              If `-1`, the output array `count` (and along with it, `unique` and `index`) will be sorted in **descending** order.<br>
    !>                              If `+1`, the output array `count` (and along with it, `unique` and `index`) will be sorted in  **ascending** order.<br>
    !>                              (**optional**, default = `0_IK`)
    !>  \param      iseq        :   The `external` user-specified function that takes two input **scalar** arguments of the same type and kind as the input `array`.<br>
    !>                              It returns a scalar `logical` of default kind \LK that is `.true.` if the two input arguments are equivalent (e.g., equal)
    !>                              according to the user-defined criterion, otherwise, it is `.false.`.<br>
    !>                              The following illustrates the generic interface of `iseq`,
    !>                              \code{.F90}
    !>                                  function iseq(element1, element2) result(equivalent)
    !>                                      use pm_kind, only: LK
    !>                                      TYPE(KIND)  , intent(in)    :: element1, element2
    !>                                      logical(LK)                 :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `setA`, which can be one of the following,
    !>                              \code{.F90}
    !>                                  use pm_kind, only: SK, IK, CK, RK
    !>                                  character(*, SK), intent(in)    :: element1, element2 ! when `array` is a string vector.
    !>                                  character(1, SK), intent(in)    :: element1, element2 ! when `array` is a string scalar.
    !>                                  integer(IK)     , intent(in)    :: element1, element2
    !>                                  logical(LK)     , intent(in)    :: element1, element2
    !>                                  complex(CK)     , intent(in)    :: element1, element2
    !>                                  real(RK)        , intent(in)    :: element1, element2
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              This user-defined equivalence check is extremely useful where a user-defined equivalence test other than exact equality
    !>                              or identity is needed, for example, when the array elements should match only within a given threshold or,
    !>                              when the case-sensitivity in character comparisons do not matter.<br>
    !>                              In such cases, the user can define a custom equivalence criterion within the user-defined external function `iseq` to achieve the goal.<br>
    !>                              (**optional**, the default equivalence operator is `.eqv.` if the input `array` is `logical`, otherwise `==`.)
    !>
    !>  \interface{setUnique}
    !>  \code{.F90}
    !>
    !>      use pm_arrayUnique, only: setUnique
    !>
    !>      call setUnique(array, unique(:), count(:), index = index(:), order = order) ! Requires allocatable `unique`, `count`, and `index` arguments.
    !>      call setUnique(array, unique(:), count(:), iseq, index = index(:), order = order) ! Requires allocatable `unique`, `count`, and `index` arguments.
    !>      call setUnique(array, unique(:), lenUnique, count(:), index = index(:), order = order) ! Avoids argument allocations for performance-critical tasks.
    !>      call setUnique(array, unique(:), lenUnique, count(:), iseq, index = index(:), order = order) ! Avoids argument allocations for performance-critical tasks.
    !>      !
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \warning
    !>  The procedures under this generic interface are `impure` when the user-specified `external` procedure `iseq` is specified as input argument.<br>
    !>  The condition `lenUnique <= len/size(unique)` must hold for the corresponding input/output arguments.<br>
    !>  The condition `lenUnique <= len/size(index)` must hold for the corresponding input/output arguments.<br>
    !>  The condition `lenUnique <= len/size(count)` must hold for the corresponding input/output arguments.<br>
    !>  \vericons
    !>
    !>  \warning
    !>  Note that in Fortran, trailing blanks are ignored in character comparison, that is, `"Fortran" == "Fortran "` yields `.true.`.
    !>
    !>  \remark
    !>  The functions under this generic interface are slightly faster than the [getUnique](@ref pm_arrayUnique::getUnique) functional implementations.
    !>
    !>  \see
    !>  [getUnique](@ref pm_arrayUnique::getUnique)<br>
    !>  [setUnique](@ref pm_arrayUnique::setUnique)<br>
    !>  [isUniqueAll](@ref pm_arrayUnique::isUniqueAll)<br>
    !>  [isUniqueAny](@ref pm_arrayUnique::isUniqueAny)<br>
    !>
    !>  \example{setUnique}
    !>  \include{lineno} example/pm_arrayUnique/setUnique/main.F90
    !>  \compilef{setUnique}
    !>  \output{setUnique}
    !>  \include{lineno} example/pm_arrayUnique/setUnique/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayUnique](@ref test_pm_arrayUnique)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.7}
    !>  \desc
    !>  The \ifort{2021.7} on Microsoft WSL (and possibly other platforms)
    !>  cannot properly call `getRemapped(unique, ...)` when `unique` is a contiguous `character` array.<br>
    !>  \remedy
    !>  For now, the remappings are done separately using the Fortran syntax to bypass the problem.<br>
    !>  This bug cost half a day of human life to identify.<br>
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to 2D input objects.<br>
    !>
    !>  \todo
    !>  \phigh To avoid extra data copy and improve performance, an extra optional `lenUnique` output argument could be added to the procedures such
    !>  that when it is present, the output arrays `unique` and `count` will not be resized from `size(array)` to the correct length `lenUnique`.<br>
    !>  In such a case, the onus would be on the user to ensure only the elements `(1:lenUnique)` of the output arrays are used for any subsequent
    !>  work as only these elements are meaningful. This would, however, come with the benefit of extra efficiency.<br>
    !>
    !>  \final{setUnique}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setUnique

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setUniArbDefCom_D0_SK5(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                            :: array
        character(:,SKG)        , intent(out)   , allocatable           :: unique
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setUniArbDefCom_D0_SK4(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                            :: array
        character(:,SKG)        , intent(out)   , allocatable           :: unique
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setUniArbDefCom_D0_SK3(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                            :: array
        character(:,SKG)        , intent(out)   , allocatable           :: unique
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setUniArbDefCom_D0_SK2(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                            :: array
        character(:,SKG)        , intent(out)   , allocatable           :: unique
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setUniArbDefCom_D0_SK1(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                            :: array
        character(:,SKG)        , intent(out)   , allocatable           :: unique
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setUniArbDefCom_D1_SK5(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous            :: array(:)
        character(*,SKG)        , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setUniArbDefCom_D1_SK4(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous            :: array(:)
        character(*,SKG)        , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setUniArbDefCom_D1_SK3(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous            :: array(:)
        character(*,SKG)        , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setUniArbDefCom_D1_SK2(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous            :: array(:)
        character(*,SKG)        , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setUniArbDefCom_D1_SK1(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous            :: array(:)
        character(*,SKG)        , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setUniArbDefCom_D1_IK5(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous            :: array(:)
        integer(IKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setUniArbDefCom_D1_IK4(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous            :: array(:)
        integer(IKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setUniArbDefCom_D1_IK3(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous            :: array(:)
        integer(IKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine

#endif

#if IK2_ENABLED
    PURE module subroutine setUniArbDefCom_D1_IK2(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous            :: array(:)
        integer(IKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setUniArbDefCom_D1_IK1(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous            :: array(:)
        integer(IKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setUniArbDefCom_D1_LK5(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous            :: array(:)
        logical(LKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setUniArbDefCom_D1_LK4(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous            :: array(:)
        logical(LKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setUniArbDefCom_D1_LK3(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous            :: array(:)
        logical(LKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setUniArbDefCom_D1_LK2(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous            :: array(:)
        logical(LKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setUniArbDefCom_D1_LK1(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous            :: array(:)
        logical(LKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUniArbDefCom_D1_CK5(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous            :: array(:)
        complex(CKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUniArbDefCom_D1_CK4(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous            :: array(:)
        complex(CKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUniArbDefCom_D1_CK3(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous            :: array(:)
        complex(CKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUniArbDefCom_D1_CK2(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous            :: array(:)
        complex(CKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUniArbDefCom_D1_CK1(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous            :: array(:)
        complex(CKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUniArbDefCom_D1_RK5(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous            :: array(:)
        real(RKG)               , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUniArbDefCom_D1_RK4(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous            :: array(:)
        real(RKG)               , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUniArbDefCom_D1_RK3(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous            :: array(:)
        real(RKG)               , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUniArbDefCom_D1_RK2(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous            :: array(:)
        real(RKG)               , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUniArbDefCom_D1_RK1(array, unique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbDefCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous            :: array(:)
        real(RKG)               , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setUniArbCusCom_D0_SK5(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                            :: array
        character(:,SKG)        , intent(out)   , allocatable           :: unique
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setUniArbCusCom_D0_SK4(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                            :: array
        character(:,SKG)        , intent(out)   , allocatable           :: unique
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setUniArbCusCom_D0_SK3(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                            :: array
        character(:,SKG)        , intent(out)   , allocatable           :: unique
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setUniArbCusCom_D0_SK2(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                            :: array
        character(:,SKG)        , intent(out)   , allocatable           :: unique
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setUniArbCusCom_D0_SK1(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                            :: array
        character(:,SKG)        , intent(out)   , allocatable           :: unique
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setUniArbCusCom_D1_SK5(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous            :: array(:)
        character(*,SKG)        , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setUniArbCusCom_D1_SK4(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous            :: array(:)
        character(*,SKG)        , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setUniArbCusCom_D1_SK3(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous            :: array(:)
        character(*,SKG)        , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setUniArbCusCom_D1_SK2(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous            :: array(:)
        character(*,SKG)        , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setUniArbCusCom_D1_SK1(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous            :: array(:)
        character(*,SKG)        , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setUniArbCusCom_D1_IK5(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous            :: array(:)
        integer(IKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setUniArbCusCom_D1_IK4(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous            :: array(:)
        integer(IKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setUniArbCusCom_D1_IK3(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous            :: array(:)
        integer(IKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine

#endif

#if IK2_ENABLED
    module subroutine setUniArbCusCom_D1_IK2(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous            :: array(:)
        integer(IKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setUniArbCusCom_D1_IK1(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous            :: array(:)
        integer(IKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setUniArbCusCom_D1_LK5(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous            :: array(:)
        logical(LKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setUniArbCusCom_D1_LK4(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous            :: array(:)
        logical(LKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setUniArbCusCom_D1_LK3(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous            :: array(:)
        logical(LKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setUniArbCusCom_D1_LK2(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous            :: array(:)
        logical(LKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setUniArbCusCom_D1_LK1(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous            :: array(:)
        logical(LKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setUniArbCusCom_D1_CK5(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous            :: array(:)
        complex(CKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setUniArbCusCom_D1_CK4(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous            :: array(:)
        complex(CKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setUniArbCusCom_D1_CK3(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous            :: array(:)
        complex(CKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setUniArbCusCom_D1_CK2(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous            :: array(:)
        complex(CKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setUniArbCusCom_D1_CK1(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous            :: array(:)
        complex(CKG)            , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setUniArbCusCom_D1_RK5(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous            :: array(:)
        real(RKG)               , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setUniArbCusCom_D1_RK4(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous            :: array(:)
        real(RKG)               , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setUniArbCusCom_D1_RK3(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous            :: array(:)
        real(RKG)               , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setUniArbCusCom_D1_RK2(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous            :: array(:)
        real(RKG)               , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setUniArbCusCom_D1_RK1(array, unique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniArbCusCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous            :: array(:)
        real(RKG)               , intent(out)   , allocatable           :: unique(:)
        integer(IK)             , intent(out)   , allocatable           :: count(:)
        type(cvi_type)          , intent(out)   , allocatable, optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
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
    PURE module subroutine setUniFixDefCom_D0_SK5(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                            :: array
        character(*,SKG)        , intent(out)                           :: unique
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setUniFixDefCom_D0_SK4(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                            :: array
        character(*,SKG)        , intent(out)                           :: unique
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setUniFixDefCom_D0_SK3(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                            :: array
        character(*,SKG)        , intent(out)                           :: unique
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setUniFixDefCom_D0_SK2(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                            :: array
        character(*,SKG)        , intent(out)                           :: unique
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setUniFixDefCom_D0_SK1(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                            :: array
        character(*,SKG)        , intent(out)                           :: unique
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setUniFixDefCom_D1_SK5(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous            :: array(:)
        character(*,SKG)        , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setUniFixDefCom_D1_SK4(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous            :: array(:)
        character(*,SKG)        , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setUniFixDefCom_D1_SK3(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous            :: array(:)
        character(*,SKG)        , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setUniFixDefCom_D1_SK2(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous            :: array(:)
        character(*,SKG)        , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setUniFixDefCom_D1_SK1(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous            :: array(:)
        character(*,SKG)        , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setUniFixDefCom_D1_IK5(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous            :: array(:)
        integer(IKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setUniFixDefCom_D1_IK4(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous            :: array(:)
        integer(IKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setUniFixDefCom_D1_IK3(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous            :: array(:)
        integer(IKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine

#endif

#if IK2_ENABLED
    PURE module subroutine setUniFixDefCom_D1_IK2(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous            :: array(:)
        integer(IKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setUniFixDefCom_D1_IK1(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous            :: array(:)
        integer(IKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setUniFixDefCom_D1_LK5(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous            :: array(:)
        logical(LKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setUniFixDefCom_D1_LK4(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous            :: array(:)
        logical(LKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setUniFixDefCom_D1_LK3(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous            :: array(:)
        logical(LKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setUniFixDefCom_D1_LK2(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous            :: array(:)
        logical(LKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setUniFixDefCom_D1_LK1(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous            :: array(:)
        logical(LKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setUniFixDefCom_D1_CK5(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous            :: array(:)
        complex(CKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setUniFixDefCom_D1_CK4(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous            :: array(:)
        complex(CKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setUniFixDefCom_D1_CK3(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous            :: array(:)
        complex(CKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setUniFixDefCom_D1_CK2(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous            :: array(:)
        complex(CKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setUniFixDefCom_D1_CK1(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous            :: array(:)
        complex(CKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUniFixDefCom_D1_RK5(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous            :: array(:)
        real(RKG)               , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUniFixDefCom_D1_RK4(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous            :: array(:)
        real(RKG)               , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUniFixDefCom_D1_RK3(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous            :: array(:)
        real(RKG)               , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUniFixDefCom_D1_RK2(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous            :: array(:)
        real(RKG)               , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUniFixDefCom_D1_RK1(array, unique, lenUnique, count, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixDefCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous            :: array(:)
        real(RKG)               , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setUniFixCusCom_D0_SK5(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                            :: array
        character(*,SKG)        , intent(out)                           :: unique
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setUniFixCusCom_D0_SK4(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                            :: array
        character(*,SKG)        , intent(out)                           :: unique
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setUniFixCusCom_D0_SK3(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                            :: array
        character(*,SKG)        , intent(out)                           :: unique
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setUniFixCusCom_D0_SK2(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                            :: array
        character(*,SKG)        , intent(out)                           :: unique
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setUniFixCusCom_D0_SK1(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                            :: array
        character(*,SKG)        , intent(out)                           :: unique
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setUniFixCusCom_D1_SK5(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous            :: array(:)
        character(*,SKG)        , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setUniFixCusCom_D1_SK4(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous            :: array(:)
        character(*,SKG)        , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setUniFixCusCom_D1_SK3(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous            :: array(:)
        character(*,SKG)        , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setUniFixCusCom_D1_SK2(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous            :: array(:)
        character(*,SKG)        , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setUniFixCusCom_D1_SK1(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous            :: array(:)
        character(*,SKG)        , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setUniFixCusCom_D1_IK5(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous            :: array(:)
        integer(IKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setUniFixCusCom_D1_IK4(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous            :: array(:)
        integer(IKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setUniFixCusCom_D1_IK3(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous            :: array(:)
        integer(IKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine

#endif

#if IK2_ENABLED
    module subroutine setUniFixCusCom_D1_IK2(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous            :: array(:)
        integer(IKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setUniFixCusCom_D1_IK1(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous            :: array(:)
        integer(IKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setUniFixCusCom_D1_LK5(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous            :: array(:)
        logical(LKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setUniFixCusCom_D1_LK4(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous            :: array(:)
        logical(LKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setUniFixCusCom_D1_LK3(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous            :: array(:)
        logical(LKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setUniFixCusCom_D1_LK2(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous            :: array(:)
        logical(LKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setUniFixCusCom_D1_LK1(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous            :: array(:)
        logical(LKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setUniFixCusCom_D1_CK5(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous            :: array(:)
        complex(CKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setUniFixCusCom_D1_CK4(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous            :: array(:)
        complex(CKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setUniFixCusCom_D1_CK3(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous            :: array(:)
        complex(CKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setUniFixCusCom_D1_CK2(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous            :: array(:)
        complex(CKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setUniFixCusCom_D1_CK1(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous            :: array(:)
        complex(CKG)            , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setUniFixCusCom_D1_RK5(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous            :: array(:)
        real(RKG)               , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setUniFixCusCom_D1_RK4(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous            :: array(:)
        real(RKG)               , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setUniFixCusCom_D1_RK3(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous            :: array(:)
        real(RKG)               , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setUniFixCusCom_D1_RK2(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous            :: array(:)
        real(RKG)               , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setUniFixCusCom_D1_RK1(array, unique, lenUnique, count, iseq, index, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUniFixCusCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous            :: array(:)
        real(RKG)               , intent(out)   , contiguous            :: unique(:)
        integer(IK)             , intent(out)                           :: lenUnique
        integer(IK)             , intent(out)   , contiguous            :: count(:)
        type(cvi_type)          , intent(out)   , contiguous , optional :: index(:)
        integer(IK)             , intent(in)                 , optional :: order
        procedure(logical(LK))                                          :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setUnique

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    !>  \brief
!    !>  Find the unique values in the input integer vector.
!    !>
!    !>  \param[in]       lenVector       :   The size of the input square matrix - `nd` by `nd`.
!    !>  \param[in]       Vector          :   The input integer vector.
!    !>  \param[out]      lenUnique       :   The length of `UniqueValue`, that is, the total number of unique values.
!    !>  \param[out]      UniqueValue     :   The vector of unique values identified in the input vector.
!    !>  \param[out]      uniqueCount     :   The counts of each unique value in the input vector.
!    !>  \param[out]      UniqueIndex     :   A jaggedArray of type [cvi_type](@ref pm_container::cvi_type) of length `lenUnique`,
!    !>                                       the ith element of which is a vector of length `uniqueCount(i)` that contains the
!    !>                                       indices of `Vector` where `UniqueValue(i)` occur (**optional**).
!    !>  \param[in]       sorting         :   The input `integer` of default kind \IK.
!    !>                                       If  `0`, the output `uniqueCount` will be not be sorted.
!    !>                                       If `-1`, the output `uniqueCount` will be sorted in **descending** order and along with it `UniqueValue` and `UniqueIndex`.
!    !>                                       If `+1`, the output `uniqueCount` will be sorted in  **ascending** order and along with it `UniqueValue` and `UniqueIndex`.
!    !>                                       (**optional**, default = `0_IK`)
!    !>
!    !>  \warnpure
!    !>
!    !>  \warning
!    !>  To avoid extra data copy and improve performance, the output arrays `uniqueCount` and `UniqueValue` will not be
!    !>  resized from `(1:lenVector)` to `(1:lenUnique)`. The onus is on the user to ensure only the elements `(1:lenUnique)`
!    !>  are used for any subsequent work as only these elements are meaningful. **However**, if `Err` argument is present,
!    !>  the two aforementioned arrays will be ordered and then automatically resized to `(1:lenUnique)`.
!    !>
!    PURE subroutine findUniqueValueCount_IK(Vector, lenUnique, UniqueValue, uniqueCount, UniqueIndex, sorting)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: findUniqueValueCount_IK
!#endif
!        use pm_container, only: IV => cvi_type
!        use pm_arraySort, only: setSorted
!        use pm_kind, only: IK
!
!        implicit none
!
!        integer(IK) , intent(in)    , contiguous                :: Vector(:)
!        integer(IK) , intent(out)                               :: lenUnique
!        integer(IK) , intent(out)   , allocatable               :: uniqueValue(:)
!        integer(IK) , intent(out)   , allocatable               :: uniqueCount(:)
!        type(IV)    , intent(out)   , allocatable   , optional  :: uniqueIndex(:)
!        integer(IK) , intent(in)                    , optional  :: sorting
!
!
!        character(*, SK), parameter     :: PROCEDURE_NAME = MODULE_NAME//SK_"@findUniqueValueCount_IK()"
!
!        integer(IK)                     :: lenVector
!        integer(IK)     , allocatable   :: indx(:)
!        integer(IK)                     :: ivec, iuniq, counter
!        integer(IK)                     :: sorting_def
!        logical(LK)                         :: isUnique
!
!        lenVector = size(Vector, kind = IK)
!
!        if (present(sorting)) then
!            sorting_def = sorting
!        else
!            sorting_def = 0_IK
!        end if
!
!        allocate(UniqueValue(lenVector))
!        allocate(uniqueCount(lenVector), source = 0_IK)
!
!        lenUnique = 0_IK
!
!        do ivec = 1, lenVector
!            isUnique = .true._LK
!            loopSearchUnique: do iuniq = 1, lenUnique
!                if (UniqueValue(iuniq)==Vector(ivec)) then
!                    uniqueCount(iuniq) = uniqueCount(iuniq) + 1
!                    isUnique = .false._LK
!                    exit loopSearchUnique
!                end if
!            end do loopSearchUnique
!            if (isUnique) then
!                lenUnique = lenUnique + 1
!                UniqueValue(lenUnique) = Vector(ivec)
!                uniqueCount(lenUnique) = uniqueCount(lenUnique) + 1
!            end if
!        end do
!
!        if (sorting_def /= 0_IK) then
!            allocate(indx(lenUnique))
!            call setSorted(uniqueCount(1:lenUnique), indx)
!            if (sorting_def < 0_IK) then
!                uniqueCount = uniqueCount(indx(lenUnique:1:-1))
!                UniqueValue = UniqueValue(indx(lenUnique:1:-1))
!            else
!                uniqueCount = uniqueCount(indx)
!                UniqueValue = UniqueValue(indx)
!            end if
!            deallocate(indx)
!        end if
!
!        if (present(UniqueIndex)) then
!            allocate(UniqueIndex(lenUnique))
!            loopUniqueIndex: do iuniq = 1, lenUnique
!                allocate(UniqueIndex(iuniq)%val(uniqueCount(iuniq)))
!                counter = 1_IK
!                loopOverVector: do ivec = 1, lenVector
!                    if (UniqueValue(iuniq) == Vector(ivec)) then
!                        UniqueIndex(iuniq)%val(counter) = ivec
!                        counter = counter + 1_IK
!                    end if
!                    if (counter > uniqueCount(iuniq)) exit loopOverVector
!                end do loopOverVector
!            end do loopUniqueIndex
!        end if
!
!    end subroutine findUniqueValueCount_IK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayUnique