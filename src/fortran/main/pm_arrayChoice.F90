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
!>  This module contains procedures and generic interfaces for selecting uniformly-distributed or
!>  arbitrarily-distributed random choices from a given list of intrinsic type of arbitrary kind.<br>
!>
!>  \see
!>  [pm_arrayRemap](@ref pm_arrayRemap)<br>
!>  [pm_arrayChange](@ref pm_arrayChange)<br>
!>  [pm_arrayChoice](@ref pm_arrayChoice)<br>
!>  [pm_arrayShuffle](@ref pm_arrayShuffle)<br>
!>  [pm_distUnif](@ref pm_distUnif)<br>
!>
!>  \test
!>  [test_pm_arrayChoice](@ref test_pm_arrayChoice)
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayChoice

    use pm_kind, only: SK, IK, LK
    use pm_distUnif, only: rngf_type, xoshiro256ssw_type

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_arrayChoice"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Select a single (or multiple) element(s) from the input array of intrinsic type of arbitrary kind randomly
    !>  uniformly or optionally according to the specified Probability Mass Function (PMF) of the input `array`.<br>
    !>
    !>  \param[in]  array   :   The input `contiguous` array of non-zero size of shape `(:)` of either <br>
    !>                          <ol>
    !>                              <li>    type `character` of kind \SKALL, or <br>
    !>                              <li>    type `logical` of kind \LKALL, or <br>
    !>                              <li>    type `integer` of kind \IKALL, or <br>
    !>                              <li>    type `complex` of kind \CKALL, or <br>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ol>
    !>                          or a scalar of,
    !>                          <ol>
    !>                              <li>    type `character` of kind \SKALL of non-zero length, <br>
    !>                              <li>    type `integer` of kind \IKALL of non-zero value, <br>
    !>                          </ol>
    !>                          whose element(s) will be selected randomly and returned in the output `choice`.<br>
    !>                          The input `array` must be of non-zero length or size (unless it is a scalar non-negative integer of arbitrary kind).<br>
    !>  \param[in]  s1      :   The input non-negative scalar `integer` of default kind \IK, representing the number of output choices
    !>                          (i.e., elements selected from the input `array`).<br>
    !>                          (**optional**. If missing, the output is a scalar (of rank `0`).)
    !>  \param[in]  unique  :   The input scalar `logical` of default kind \LK.<br>
    !>                          If `.true.`, the elements of the output `choice` will be **uniquely** selected from the input `array`.<br>
    !>                          (**optional**, default = `.false._LK`. It can be present **if and only if** the input argument `s1` is also present.)
    !>
    !>  \return
    !>  `choice`            :   The output object of the same type and kind as the input `array`
    !>                          whose value is randomly selected from the elements of `array`.<br>
    !>                          <ol>
    !>                              <li>    It is a scalar (of rank `0`, or a scalar `character` of length `1`) if the input argument `s1` is missing.
    !>                              <li>    It is an array of size `s1` (or a scalar `character` of length `s1`) if the input argument `s1` is present.
    !>                          </ol>
    !>
    !>  \interface{getChoice}
    !>  \code{.F90}
    !>
    !>      use pm_arrayChoice, only: getChoice
    !>
    !>      choice = getChoice(array)
    !>      choice(1:s1) = getChoice(array, s1, unique = unique)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < size(array)` for non-`character` input `array` or `0 < len(array)` for `character` input `array` must hold.<br>
    !>  \vericons
    !>
    !>  \remark
    !>  The procedures under this generic interface provide only a convenient functional interface to random choices.<br>
    !>  For high performance applications, use the subroutine interface [setChoice](@ref pm_arrayChoice::setChoice).<br>
    !>
    !>  \devnote
    !>  The choice of the input argument name `s1` is deliberate to allow future
    !>  expansion of the interface to output `choice` of higher rank than `1`.<br>
    !>
    !>  \warning
    !>  The condition `0 < size(array)` for non-`character` input `array` or `0 < len(array)` for `character` input `array` must hold.<br>
    !>  The condition `0 < size(choice)` for non-`character` input `choice` or `0 < len(choice)` for `character` input `choice` must hold.<br>
    !>  The condition `len(choice) == len(array)` for an input `choice` vector of `character` must hold.<br>
    !>  \vericons
    !>
    !>  \see
    !>  [isHead](@ref pm_distBern::isHead)<br>
    !>  [getChange](@ref pm_arrayChange::getChange)<br>
    !>  [setChange](@ref pm_arrayChange::setChange)<br>
    !>  [getChoice](@ref pm_arrayChoice::getChoice)<br>
    !>  [setChoice](@ref pm_arrayChoice::setChoice)<br>
    !>  [getUnifRand](@ref pm_distUnif::getUnifRand)<br>
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand)<br>
    !>  [getShuffled](@ref pm_arrayShuffle::getShuffled)<br>
    !>  [setShuffled](@ref pm_arrayShuffle::setShuffled)<br>
    !>  [getRemapped](@ref pm_arrayRemap::getRemapped)<br>
    !>  [setRemapped](@ref pm_arrayRemap::setRemapped)<br>
    !>
    !>  \example{getChoice}
    !>  \include{lineno} example/pm_arrayChoice/getChoice/main.F90
    !>  \compilef{getChoice}
    !>  \output{getChoice}
    !>  \include{lineno} example/pm_arrayChoice/getChoice/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayChoice](@ref test_pm_arrayChoice)
    !>
    !>  \finmain{getChoice}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    interface getChoice

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getChoiceRNGD_D0_D0_SK5(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D0_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(in)                    :: array
        character(1,SKC)                                            :: choice
    end function
#endif

#if SK4_ENABLED
    module function getChoiceRNGD_D0_D0_SK4(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D0_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(in)                    :: array
        character(1,SKC)                                            :: choice
    end function
#endif

#if SK3_ENABLED
    module function getChoiceRNGD_D0_D0_SK3(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D0_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(in)                    :: array
        character(1,SKC)                                            :: choice
    end function
#endif

#if SK2_ENABLED
    module function getChoiceRNGD_D0_D0_SK2(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D0_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(in)                    :: array
        character(1,SKC)                                            :: choice
    end function
#endif

#if SK1_ENABLED
    module function getChoiceRNGD_D0_D0_SK1(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D0_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(in)                    :: array
        character(1,SKC)                                            :: choice
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getChoiceRNGD_D1_D0_SK5(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC)                                :: choice
    end function
#endif

#if SK4_ENABLED
    module function getChoiceRNGD_D1_D0_SK4(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC)                                :: choice
    end function
#endif

#if SK3_ENABLED
    module function getChoiceRNGD_D1_D0_SK3(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC)                                :: choice
    end function
#endif

#if SK2_ENABLED
    module function getChoiceRNGD_D1_D0_SK2(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC)                                :: choice
    end function
#endif

#if SK1_ENABLED
    module function getChoiceRNGD_D1_D0_SK1(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC)                                :: choice
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getChoiceRNGD_D1_D0_IK5(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                                                :: choice
    end function
#endif

#if IK4_ENABLED
    module function getChoiceRNGD_D1_D0_IK4(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                                                :: choice
    end function
#endif

#if IK3_ENABLED
    module function getChoiceRNGD_D1_D0_IK3(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                                                :: choice
    end function
#endif

#if IK2_ENABLED
    module function getChoiceRNGD_D1_D0_IK2(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                                                :: choice
    end function
#endif

#if IK1_ENABLED
    module function getChoiceRNGD_D1_D0_IK1(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                                                :: choice
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getChoiceRNGD_D1_D0_LK5(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                                                :: choice
    end function
#endif

#if LK4_ENABLED
    module function getChoiceRNGD_D1_D0_LK4(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                                                :: choice
    end function
#endif

#if LK3_ENABLED
    module function getChoiceRNGD_D1_D0_LK3(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                                                :: choice
    end function
#endif

#if LK2_ENABLED
    module function getChoiceRNGD_D1_D0_LK2(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                                                :: choice
    end function
#endif

#if LK1_ENABLED
    module function getChoiceRNGD_D1_D0_LK1(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                                                :: choice
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getChoiceRNGD_D1_D0_CK5(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                                                :: choice
    end function
#endif

#if CK4_ENABLED
    module function getChoiceRNGD_D1_D0_CK4(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                                                :: choice
    end function
#endif

#if CK3_ENABLED
    module function getChoiceRNGD_D1_D0_CK3(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                                                :: choice
    end function
#endif

#if CK2_ENABLED
    module function getChoiceRNGD_D1_D0_CK2(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                                                :: choice
    end function
#endif

#if CK1_ENABLED
    module function getChoiceRNGD_D1_D0_CK1(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                                                :: choice
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getChoiceRNGD_D1_D0_RK5(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                                                   :: choice
    end function
#endif

#if RK4_ENABLED
    module function getChoiceRNGD_D1_D0_RK4(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                                                   :: choice
    end function
#endif

#if RK3_ENABLED
    module function getChoiceRNGD_D1_D0_RK3(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                                                   :: choice
    end function
#endif

#if RK2_ENABLED
    module function getChoiceRNGD_D1_D0_RK2(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                                                   :: choice
    end function
#endif

#if RK1_ENABLED
    module function getChoiceRNGD_D1_D0_RK1(array) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                                                   :: choice
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getChoiceRNGD_D0_S1_SK5(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D0_S1_SK5
#endif
        use pm_kind, only: SKC => SK5
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)                    :: array
        character(s1,SKC)                                           :: choice
    end function
#endif

#if SK4_ENABLED
    module function getChoiceRNGD_D0_S1_SK4(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D0_S1_SK4
#endif
        use pm_kind, only: SKC => SK4
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)                    :: array
        character(s1,SKC)                                           :: choice
    end function
#endif

#if SK3_ENABLED
    module function getChoiceRNGD_D0_S1_SK3(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D0_S1_SK3
#endif
        use pm_kind, only: SKC => SK3
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)                    :: array
        character(s1,SKC)                                           :: choice
    end function
#endif

#if SK2_ENABLED
    module function getChoiceRNGD_D0_S1_SK2(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D0_S1_SK2
#endif
        use pm_kind, only: SKC => SK2
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)                    :: array
        character(s1,SKC)                                           :: choice
    end function
#endif

#if SK1_ENABLED
    module function getChoiceRNGD_D0_S1_SK1(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D0_S1_SK1
#endif
        use pm_kind, only: SKC => SK1
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)                    :: array
        character(s1,SKC)                                           :: choice
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getChoiceRNGD_D1_D1_SK5(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC)                                :: choice(s1)
    end function
#endif

#if SK4_ENABLED
    module function getChoiceRNGD_D1_D1_SK4(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC)                                :: choice(s1)
    end function
#endif

#if SK3_ENABLED
    module function getChoiceRNGD_D1_D1_SK3(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC)                                :: choice(s1)
    end function
#endif

#if SK2_ENABLED
    module function getChoiceRNGD_D1_D1_SK2(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC)                                :: choice(s1)
    end function
#endif

#if SK1_ENABLED
    module function getChoiceRNGD_D1_D1_SK1(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC)                                :: choice(s1)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getChoiceRNGD_D1_D1_IK5(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                                                :: choice(s1)
    end function
#endif

#if IK4_ENABLED
    module function getChoiceRNGD_D1_D1_IK4(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                                                :: choice(s1)
    end function
#endif

#if IK3_ENABLED
    module function getChoiceRNGD_D1_D1_IK3(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                                                :: choice(s1)
    end function
#endif

#if IK2_ENABLED
    module function getChoiceRNGD_D1_D1_IK2(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                                                :: choice(s1)
    end function
#endif

#if IK1_ENABLED
    module function getChoiceRNGD_D1_D1_IK1(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                                                :: choice(s1)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getChoiceRNGD_D1_D1_LK5(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                                                :: choice(s1)
    end function
#endif

#if LK4_ENABLED
    module function getChoiceRNGD_D1_D1_LK4(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                                                :: choice(s1)
    end function
#endif

#if LK3_ENABLED
    module function getChoiceRNGD_D1_D1_LK3(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                                                :: choice(s1)
    end function
#endif

#if LK2_ENABLED
    module function getChoiceRNGD_D1_D1_LK2(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                                                :: choice(s1)
    end function
#endif

#if LK1_ENABLED
    module function getChoiceRNGD_D1_D1_LK1(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                                                :: choice(s1)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getChoiceRNGD_D1_D1_CK5(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                                                :: choice(s1)
    end function
#endif

#if CK4_ENABLED
    module function getChoiceRNGD_D1_D1_CK4(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                                                :: choice(s1)
    end function
#endif

#if CK3_ENABLED
    module function getChoiceRNGD_D1_D1_CK3(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                                                :: choice(s1)
    end function
#endif

#if CK2_ENABLED
    module function getChoiceRNGD_D1_D1_CK2(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                                                :: choice(s1)
    end function
#endif

#if CK1_ENABLED
    module function getChoiceRNGD_D1_D1_CK1(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                                                :: choice(s1)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getChoiceRNGD_D1_D1_RK5(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                                                   :: choice(s1)
    end function
#endif

#if RK4_ENABLED
    module function getChoiceRNGD_D1_D1_RK4(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                                                   :: choice(s1)
    end function
#endif

#if RK3_ENABLED
    module function getChoiceRNGD_D1_D1_RK3(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                                                   :: choice(s1)
    end function
#endif

#if RK2_ENABLED
    module function getChoiceRNGD_D1_D1_RK2(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                                                   :: choice(s1)
    end function
#endif

#if RK1_ENABLED
    module function getChoiceRNGD_D1_D1_RK1(array, s1, unique) result(choice)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChoiceRNGD_D1_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)                 , intent(in)                    :: s1
        logical(LK)                 , intent(in)    , optional      :: unique
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                                                   :: choice(s1)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Select a single (or multiple) element(s) from the input array of intrinsic type of arbitrary kind randomly
    !>  uniformly or optionally according to the specified Probability Mass Function (PMF) of the input `array`.<br>
    !>
    !>  \param[inout]   rng     :   The input/output scalar that can be an object of,
    !>                              <ol>
    !>                                  <li>    type [rngf_type](@ref pm_distUnif::rngf_type),
    !>                                          implying the use of intrinsic Fortran uniform RNG.<br>
    !>                                  <li>    type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type),
    !>                                          implying the use of [xoshiro256**](https://prng.di.unimi.it/) uniform RNG.<br>
    !>                              </ol>
    !>  \param[out]     choice  :   The output scalar or vector of the same type and kind as the input `array`
    !>                              whose value is randomly selected from the elements of `array`.<br>
    !>  \param[in]      array   :   The input `contiguous` array of non-zero size of shape `(:)` of either <br>
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL, or <br>
    !>                                  <li>    type `logical` of kind \LKALL, or <br>
    !>                                  <li>    type `integer` of kind \IKALL, or <br>
    !>                                  <li>    type `complex` of kind \CKALL, or <br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ol>
    !>                              or a scalar of,
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL of non-zero length, <br>
    !>                                  <li>    type `integer` of kind \IKALL of non-zero value, <br>
    !>                              </ol>
    !>                              whose element(s) will be selected randomly and returned in the output `choice`.<br>
    !>                              The input `array` must be of non-zero length or size (unless it is a scalar non-negative integer of arbitrary kind).<br>
    !>  \param[in]      unique  :   The input scalar `logical` of default kind \LK.<br>
    !>                              If `.true.`, the elements of the output `choice` will be **uniquely** selected from the input `array`.<br>
    !>                              If `.true.`, the size/length of the output `choice` must be larger than or equal to the size/length of the input `array`,
    !>                              otherwise, unique selection impossible.<br>
    !>                              (**optional**, default = `.false._LK`. It can be present **if and only if** the input argument is array-like (i.e., either scalar `character` or vector).)
    !>
    !>  \interface{setChoice}
    !>  \code{.F90}
    !>
    !>      use pm_arrayChoice, only: setChoice
    !>
    !>      call setChoice(rng, choice, array) ! `choice` is a scalar.
    !>      call setChoice(rng, choice, array, unique = unique) ! `choice` is array-like (sliceable, e.g., vector, or a scalar `character`) and longer than or equal-length to `array`.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < size(array)` for non-`character` input `array` or `0 < len(array)` for `character` input `array` must hold.<br>
    !>  The condition `len(choice) == len(array)` for an input `choice` vector of `character` must hold.<br>
    !>  \vericons
    !>
    !>  \see
    !>  [isHead](@ref pm_distBern::isHead)<br>
    !>  [getChange](@ref pm_arrayChange::getChange)<br>
    !>  [setChange](@ref pm_arrayChange::setChange)<br>
    !>  [getChoice](@ref pm_arrayChoice::getChoice)<br>
    !>  [setChoice](@ref pm_arrayChoice::setChoice)<br>
    !>  [getUnifRand](@ref pm_distUnif::getUnifRand)<br>
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand)<br>
    !>  [getShuffled](@ref pm_arrayShuffle::getShuffled)<br>
    !>  [setShuffled](@ref pm_arrayShuffle::setShuffled)<br>
    !>  [getRemapped](@ref pm_arrayRemap::getRemapped)<br>
    !>  [setRemapped](@ref pm_arrayRemap::setRemapped)<br>
    !>
    !>  \example{setChoice}
    !>  \include{lineno} example/pm_arrayChoice/setChoice/main.F90
    !>  \compilef{setChoice}
    !>  \output{setChoice}
    !>  \include{lineno} example/pm_arrayChoice/setChoice/main.out.F90
    !>
    !>  \finmain{setChoice}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    interface setChoice

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setChoiceRNGF_D0_D0_SK5(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D0_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)                    :: array
        character(*,SKC)            , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setChoiceRNGF_D0_D0_SK4(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D0_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)                    :: array
        character(*,SKC)            , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setChoiceRNGF_D0_D0_SK3(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D0_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)                    :: array
        character(*,SKC)            , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setChoiceRNGF_D0_D0_SK2(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D0_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)                    :: array
        character(*,SKC)            , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setChoiceRNGF_D0_D0_SK1(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D0_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)                    :: array
        character(*,SKC)            , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setChoiceRNGF_D1_D0_SK5(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC), intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setChoiceRNGF_D1_D0_SK4(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC), intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setChoiceRNGF_D1_D0_SK3(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC), intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setChoiceRNGF_D1_D0_SK2(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC), intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setChoiceRNGF_D1_D0_SK1(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC), intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setChoiceRNGF_D1_D0_IK5(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setChoiceRNGF_D1_D0_IK4(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setChoiceRNGF_D1_D0_IK3(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setChoiceRNGF_D1_D0_IK2(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setChoiceRNGF_D1_D0_IK1(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setChoiceRNGF_D1_D0_LK5(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setChoiceRNGF_D1_D0_LK4(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setChoiceRNGF_D1_D0_LK3(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setChoiceRNGF_D1_D0_LK2(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setChoiceRNGF_D1_D0_LK1(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setChoiceRNGF_D1_D0_CK5(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setChoiceRNGF_D1_D0_CK4(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setChoiceRNGF_D1_D0_CK3(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setChoiceRNGF_D1_D0_CK2(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setChoiceRNGF_D1_D0_CK1(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setChoiceRNGF_D1_D0_RK5(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                   , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setChoiceRNGF_D1_D0_RK4(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                   , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setChoiceRNGF_D1_D0_RK3(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                   , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setChoiceRNGF_D1_D0_RK2(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                   , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setChoiceRNGF_D1_D0_RK1(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                   , intent(out)                   :: choice
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setChoiceRNGF_D1_D1_SK5(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC), intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setChoiceRNGF_D1_D1_SK4(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC), intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setChoiceRNGF_D1_D1_SK3(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC), intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setChoiceRNGF_D1_D1_SK2(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC), intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setChoiceRNGF_D1_D1_SK1(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC), intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setChoiceRNGF_D1_D1_IK5(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                , intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setChoiceRNGF_D1_D1_IK4(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                , intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setChoiceRNGF_D1_D1_IK3(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                , intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setChoiceRNGF_D1_D1_IK2(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                , intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setChoiceRNGF_D1_D1_IK1(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                , intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setChoiceRNGF_D1_D1_LK5(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                , intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setChoiceRNGF_D1_D1_LK4(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                , intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setChoiceRNGF_D1_D1_LK3(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                , intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setChoiceRNGF_D1_D1_LK2(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                , intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setChoiceRNGF_D1_D1_LK1(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                , intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setChoiceRNGF_D1_D1_CK5(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        logical(LK)                 , intent(in)    , optional      :: unique
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                , intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setChoiceRNGF_D1_D1_CK4(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        logical(LK)                 , intent(in)    , optional      :: unique
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                , intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setChoiceRNGF_D1_D1_CK3(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        logical(LK)                 , intent(in)    , optional      :: unique
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                , intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setChoiceRNGF_D1_D1_CK2(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        logical(LK)                 , intent(in)    , optional      :: unique
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                , intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setChoiceRNGF_D1_D1_CK1(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        logical(LK)                 , intent(in)    , optional      :: unique
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                , intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setChoiceRNGF_D1_D1_RK5(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        logical(LK)                 , intent(in)    , optional      :: unique
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                   , intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setChoiceRNGF_D1_D1_RK4(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        logical(LK)                 , intent(in)    , optional      :: unique
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                   , intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setChoiceRNGF_D1_D1_RK3(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        logical(LK)                 , intent(in)    , optional      :: unique
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                   , intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setChoiceRNGF_D1_D1_RK2(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        logical(LK)                 , intent(in)    , optional      :: unique
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                   , intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setChoiceRNGF_D1_D1_RK1(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGF_D1_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        logical(LK)                 , intent(in)    , optional      :: unique
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                   , intent(out)   , contiguous    :: choice(:)
        type(rngf_type)             , intent(in)                    :: rng
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
    PURE module subroutine setChoiceRNGX_D0_D0_SK5(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D0_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)                    :: array
        character(*,SKC)            , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setChoiceRNGX_D0_D0_SK4(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D0_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)                    :: array
        character(*,SKC)            , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setChoiceRNGX_D0_D0_SK3(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D0_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)                    :: array
        character(*,SKC)            , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setChoiceRNGX_D0_D0_SK2(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D0_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)                    :: array
        character(*,SKC)            , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setChoiceRNGX_D0_D0_SK1(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D0_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)                    :: array
        character(*,SKC)            , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_SK5(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC), intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_SK4(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC), intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_SK3(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC), intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_SK2(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC), intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_SK1(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC), intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_IK5(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_IK4(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_IK3(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_IK2(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_IK1(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_LK5(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_LK4(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_LK3(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_LK2(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_LK1(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_CK5(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_CK4(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_CK3(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_CK2(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_CK1(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_RK5(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                   , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_RK4(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                   , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_RK3(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                   , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_RK2(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                   , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D0_RK1(rng, choice, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D0_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                   , intent(out)                   :: choice
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_SK5(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC), intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_SK4(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC), intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_SK3(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC), intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_SK2(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC), intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_SK1(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        logical(LK)                 , intent(in)    , optional      :: unique
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        character(len(array,IK),SKC), intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_IK5(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                , intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_IK4(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                , intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_IK3(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                , intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_IK2(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                , intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_IK1(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IKC)                , intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_LK5(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                , intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_LK4(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                , intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_LK3(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                , intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_LK2(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                , intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_LK1(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        logical(LKC)                , intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_CK5(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        logical(LK)                 , intent(in)    , optional      :: unique
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                , intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_CK4(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        logical(LK)                 , intent(in)    , optional      :: unique
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                , intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_CK3(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        logical(LK)                 , intent(in)    , optional      :: unique
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                , intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_CK2(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        logical(LK)                 , intent(in)    , optional      :: unique
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                , intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_CK1(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        logical(LK)                 , intent(in)    , optional      :: unique
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        complex(CKC)                , intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_RK5(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        logical(LK)                 , intent(in)    , optional      :: unique
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                   , intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_RK4(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        logical(LK)                 , intent(in)    , optional      :: unique
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                   , intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_RK3(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        logical(LK)                 , intent(in)    , optional      :: unique
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                   , intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_RK2(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        logical(LK)                 , intent(in)    , optional      :: unique
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                   , intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setChoiceRNGX_D1_D1_RK1(rng, choice, array, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoiceRNGX_D1_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        logical(LK)                 , intent(in)    , optional      :: unique
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        real(RKC)                   , intent(out)   , contiguous    :: choice(:)
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
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

end module pm_arrayChoice ! LCOV_EXCL_LINE