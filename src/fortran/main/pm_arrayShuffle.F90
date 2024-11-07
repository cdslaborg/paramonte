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
!>  This module contains procedures and generic interfaces for shuffling arrays of various types.
!>
!>  \benchmarks
!>
!>  \benchmark{getShuffled_vs_setShuffled, The runtime performance of [getShuffled](@ref pm_arrayShuffle::getShuffled) vs. [setShuffled](@ref pm_arrayShuffle::setShuffled)}
!>  \include{lineno} benchmark/pm_arrayShuffle/getShuffled_vs_setShuffled/main.F90
!>  \compilefb{getShuffled_vs_setShuffled}
!>  \postprocb{getShuffled_vs_setShuffled}
!>  \include{lineno} benchmark/pm_arrayShuffle/getShuffled_vs_setShuffled/main.py
!>  \visb{getShuffled_vs_setShuffled}
!>  \image html benchmark/pm_arrayShuffle/getShuffled_vs_setShuffled/benchmark.getShuffled_vs_setShuffled.runtime.png width=1000
!>  \image html benchmark/pm_arrayShuffle/getShuffled_vs_setShuffled/benchmark.getShuffled_vs_setShuffled.runtime.ratio.png width=1000
!>  \moralb{getShuffled_vs_setShuffled}
!>      -#  The procedures under the generic interface [getShuffled](@ref pm_arrayShuffle::getShuffled) are functions while
!>          the procedures under the generic interface [setShuffled](@ref pm_arrayShuffle::setShuffled) are subroutines.<br>
!>          The current implementation of the functional interface requires making a copy of the input array that is
!>          subsequently passed to the subroutine interface for random shuffling.<br>
!>          As such, the observed performance degradation is expected.<br>
!>          Note that an extra copy of the output array to the user-specified object is also needed,
!>          making the functional interface two-copies more expensive than the subroutine interface.<br>
!>
!>  \see
!>  [pm_arrayRemap](@ref pm_arrayRemap)<br>
!>  [pm_arrayChange](@ref pm_arrayChange)<br>
!>  [pm_arrayChoice](@ref pm_arrayChoice)<br>
!>  [pm_arrayShuffle](@ref pm_arrayShuffle)<br>
!>  [pm_distUnif](@ref pm_distUnif)<br>
!>
!>  \test
!>  [test_pm_arrayShuffle](@ref test_pm_arrayShuffle)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:20 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayShuffle

    use pm_kind, only: SK, IK
#if PDT_ENABLED
    use pm_container, only: css_pdt
#endif
    use pm_container, only: css_type
    use pm_distUnif, only: rngf_type, xoshiro256ssw_type

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_arrayShuffle"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Perform an unbiased random shuffling of the input array, known as the **Knuth** or
    !>  [**Fisher-Yates** shuffle](https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle),
    !>  and generate and return the new reshuffled array.
    !>
    !>  \param[in]   array  :   The input `contiguous` array of shape `(:)` of either <br>
    !>                          <ol>
    !>                              <li>    type `character` of kind \SKALL,<br>
    !>                              <li>    type `logical` of kind \LKALL,<br>
    !>                              <li>    type `integer` of kind \IKALL,<br>
    !>                              <li>    type `complex` of kind \CKALL,<br>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                              <li>    type [css_pdt](@ref pm_container::css_pdt), <br>
    !>                              <li>    type [css_type](@ref pm_container::css_type), <br>
    !>                          </ol>
    !>                          or,
    !>                          <ol>
    !>                              <li>    a scalar assumed-length `character` of kind \SKALL, <br>
    !>                          </ol>
    !>                          whose elements will be shuffled uniformly-randomly in the output `arrayShuffled` on return.<br>
    !>  \param[in]  count   :   The input positive scalar `integer` of default kind \IK, containing the number of
    !>                          elements of the unique uniformly-random draws (shuffled elements) from the input `array`.<br>
    !>                          The specified `count` must not be larger than the length of the input sequence `array`.<br>
    !>                          (**optional**, default = `len(array)` for scalar `character` input `array`, otherwise `size(array)`.)
    !>
    !>  \return
    !>  `arrayShuffled`     :   The output `allocatable` object of the same type, kind, rank, and shape as the input `array`
    !>                          of size `(1:count)` whose elements are set from the randomly-shuffled elements of `array`.<br>
    !>
    !>  \interface{getShuffled}
    !>  \code{.F90}
    !>
    !>      use pm_arrayShuffle, only: getShuffled
    !>
    !>      arrayShuffled = getShuffled(array, count = count)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= count` must hold for the corresponding input arguments.<br>
    !>  The condition `count <= lenArray` must hold for the corresponding input arguments where `lenArray` represents the length of the input sequence.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \remark
    !>  The procedures under this generic interface provide only a convenient functional interface to shuffling arrays.<br>
    !>  For high performance applications, use the subroutine interface [setShuffled](@ref pm_arrayShuffle::setShuffled).<br>
    !>
    !>  \see
    !>  [getChoice](@ref pm_arrayChoice::getChoice)<br>
    !>  [setChoice](@ref pm_arrayChoice::setChoice)<br>
    !>  [setShuffled](@ref pm_arrayShuffle::setShuffled)<br>
    !>  [getRemapped](@ref pm_arrayRemap::getRemapped)<br>
    !>  [setRemapped](@ref pm_arrayRemap::setRemapped)<br>
    !>  [getReversed](@ref pm_arrayReverse::getReversed)<br>
    !>  [setReversed](@ref pm_arrayReverse::setReversed)<br>
    !>
    !>  \example{getShuffled}
    !>  \include{lineno} example/pm_arrayShuffle/getShuffled/main.F90
    !>  \compilef{getShuffled}
    !>  \output{getShuffled}
    !>  \include{lineno} example/pm_arrayShuffle/getShuffled/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayShuffle](@ref test_pm_arrayShuffle)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to 2D input objects.<br>
    !>
    !>  \final{getShuffled}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getShuffled

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getShuffledRNGD_D0_SK5(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        character(:,SKG)            , allocatable                   :: arrayShuffled
        integer(IK)                 , intent(in)    , optional      :: count
    end function
#endif

#if SK4_ENABLED
    module function getShuffledRNGD_D0_SK4(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        character(:,SKG)            , allocatable                   :: arrayShuffled
        integer(IK)                 , intent(in)    , optional      :: count
    end function
#endif

#if SK3_ENABLED
    module function getShuffledRNGD_D0_SK3(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        character(:,SKG)            , allocatable                   :: arrayShuffled
        integer(IK)                 , intent(in)    , optional      :: count
    end function
#endif

#if SK2_ENABLED
    module function getShuffledRNGD_D0_SK2(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        character(:,SKG)            , allocatable                   :: arrayShuffled
        integer(IK)                 , intent(in)    , optional      :: count
    end function
#endif

#if SK1_ENABLED
    module function getShuffledRNGD_D0_SK1(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        character(:,SKG)            , allocatable                   :: arrayShuffled
        integer(IK)                 , intent(in)    , optional      :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getShuffledRNGD_D1_SK5(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        character(len(array,IK),SKG), allocatable                   :: arrayShuffled(:)
    end function
#endif

#if SK4_ENABLED
    module function getShuffledRNGD_D1_SK4(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        character(len(array,IK),SKG), allocatable                   :: arrayShuffled(:)
    end function
#endif

#if SK3_ENABLED
    module function getShuffledRNGD_D1_SK3(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        character(len(array,IK),SKG), allocatable                   :: arrayShuffled(:)
    end function
#endif

#if SK2_ENABLED
    module function getShuffledRNGD_D1_SK2(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        character(len(array,IK),SKG), allocatable                   :: arrayShuffled(:)
    end function
#endif

#if SK1_ENABLED
    module function getShuffledRNGD_D1_SK1(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        character(len(array,IK),SKG), allocatable                   :: arrayShuffled(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getShuffledRNGD_D1_IK5(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        integer(IKG)                , allocatable                   :: arrayShuffled(:)
    end function
#endif

#if IK4_ENABLED
    module function getShuffledRNGD_D1_IK4(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        integer(IKG)                , allocatable                   :: arrayShuffled(:)
    end function
#endif

#if IK3_ENABLED
    module function getShuffledRNGD_D1_IK3(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        integer(IKG)                , allocatable                   :: arrayShuffled(:)
    end function
#endif

#if IK2_ENABLED
    module function getShuffledRNGD_D1_IK2(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        integer(IKG)                , allocatable                   :: arrayShuffled(:)
    end function
#endif

#if IK1_ENABLED
    module function getShuffledRNGD_D1_IK1(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        integer(IKG)                , allocatable                   :: arrayShuffled(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getShuffledRNGD_D1_LK5(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        logical(LKG)                , allocatable                   :: arrayShuffled(:)
    end function
#endif

#if LK4_ENABLED
    module function getShuffledRNGD_D1_LK4(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        logical(LKG)                , allocatable                   :: arrayShuffled(:)
    end function
#endif

#if LK3_ENABLED
    module function getShuffledRNGD_D1_LK3(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        logical(LKG)                , allocatable                   :: arrayShuffled(:)
    end function
#endif

#if LK2_ENABLED
    module function getShuffledRNGD_D1_LK2(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        logical(LKG)                , allocatable                   :: arrayShuffled(:)
    end function
#endif

#if LK1_ENABLED
    module function getShuffledRNGD_D1_LK1(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        logical(LKG)                , allocatable                   :: arrayShuffled(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getShuffledRNGD_D1_CK5(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        complex(CKG)                , allocatable                   :: arrayShuffled(:)
    end function
#endif

#if CK4_ENABLED
    module function getShuffledRNGD_D1_CK4(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        complex(CKG)                , allocatable                   :: arrayShuffled(:)
    end function
#endif

#if CK3_ENABLED
    module function getShuffledRNGD_D1_CK3(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        complex(CKG)                , allocatable                   :: arrayShuffled(:)
    end function
#endif

#if CK2_ENABLED
    module function getShuffledRNGD_D1_CK2(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        complex(CKG)                , allocatable                   :: arrayShuffled(:)
    end function
#endif

#if CK1_ENABLED
    module function getShuffledRNGD_D1_CK1(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        complex(CKG)                , allocatable                   :: arrayShuffled(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getShuffledRNGD_D1_RK5(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        real(RKG)                   , allocatable                   :: arrayShuffled(:)
    end function
#endif

#if RK4_ENABLED
    module function getShuffledRNGD_D1_RK4(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        real(RKG)                   , allocatable                   :: arrayShuffled(:)
    end function
#endif

#if RK3_ENABLED
    module function getShuffledRNGD_D1_RK3(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        real(RKG)                   , allocatable                   :: arrayShuffled(:)
    end function
#endif

#if RK2_ENABLED
    module function getShuffledRNGD_D1_RK2(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        real(RKG)                   , allocatable                   :: arrayShuffled(:)
    end function
#endif

#if RK1_ENABLED
    module function getShuffledRNGD_D1_RK1(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        real(RKG)                   , allocatable                   :: arrayShuffled(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module function getShuffledRNGD_D1_PSSK5(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(css_pdt(SKG))          , allocatable                   :: arrayShuffled(:)
    end function
#endif

#if SK4_ENABLED
    module function getShuffledRNGD_D1_PSSK4(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(css_pdt(SKG))          , allocatable                   :: arrayShuffled(:)
    end function
#endif

#if SK3_ENABLED
    module function getShuffledRNGD_D1_PSSK3(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(css_pdt(SKG))          , allocatable                   :: arrayShuffled(:)
    end function
#endif

#if SK2_ENABLED
    module function getShuffledRNGD_D1_PSSK2(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(css_pdt(SKG))          , allocatable                   :: arrayShuffled(:)
    end function
#endif

#if SK1_ENABLED
    module function getShuffledRNGD_D1_PSSK1(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(css_pdt(SKG))          , allocatable                   :: arrayShuffled(:)
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function getShuffledRNGD_D1_BSSK(array, count) result(arrayShuffled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShuffledRNGD_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(css_type)              , allocatable                   :: arrayShuffled(:)
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Perform an unbiased random shuffling of the input array,
    !>  known as the **Knuth** or [**Fisher-Yates** shuffle](https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle).
    !>
    !>  \param[inout]   rng     :   The input/output scalar that can be an object of,
    !>                              <ol>
    !>                                  <li>    type [rngf_type](@ref pm_distUnif::rngf_type),
    !>                                          implying the use of intrinsic Fortran uniform RNG.<br>
    !>                                  <li>    type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type),
    !>                                          implying the use of [xoshiro256**](https://prng.di.unimi.it/) uniform RNG.<br>
    !>                              </ol>
    !>                              (**optional**, default = [rngf_type](@ref pm_distUnif::rngf_type), implying the use of the intrinsic Fortran URNG.)
    !>  \param[inout]   array   :   The input/output `contiguous` array of shape `(:)` of either <br>
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL,<br>
    !>                                  <li>    type `logical` of kind \LKALL,<br>
    !>                                  <li>    type `integer` of kind \IKALL,<br>
    !>                                  <li>    type `complex` of kind \CKALL,<br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ol>
    !>                              or,
    !>                              <ol>
    !>                                  <li>    a scalar assumed-length `character` of kind \SKALL, <br>
    !>                              </ol>
    !>                              whose first `(1:count)` elements will be shuffled uniformly-randomly on return.<br>
    !>                              The rest of the elements of `array(count + 1:)` will be shuffled but non-uniformly.
    !>  \param[in]      count   :   The input positive scalar `integer` of default kind \IK, containing the number of
    !>                              elements of the unique uniformly-random draws (shuffled elements) from the input `array`.<br>
    !>                              If specified, only the first `count` elements of the input/output sequence `array` are guaranteed to be randomly **uniformly** shuffled.<br>
    !>                              The specified `count` must not be larger than the length of the input sequence `array`.<br>
    !>                              (**optional**, default = `len(array)` for scalar `character` input `array`, otherwise `size(array)`.)
    !>
    !>  \interface{setShuffled}
    !>  \code{.F90}
    !>
    !>      use pm_arrayShuffle, only: setShuffled
    !>
    !>      call setShuffled(array, count = count)
    !>      call setShuffled(rng, array, count = count)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= count` must hold for the corresponding input arguments.<br>
    !>  The condition `count <= lenArray` must hold for the corresponding input arguments where `lenArray` represents the length of the input sequence.<br>
    !>  \vericons
    !>
    !>  \impure
    !>  The procedures of this generic interface become `pure` in release build mode when the input
    !>  argument `rng` is set to an object of type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type).<br>
    !>
    !>  \see
    !>  [getChoice](@ref pm_arrayChoice::getChoice)<br>
    !>  [setChoice](@ref pm_arrayChoice::setChoice)<br>
    !>  [getRemapped](@ref pm_arrayRemap::getRemapped)<br>
    !>  [setRemapped](@ref pm_arrayRemap::setRemapped)<br>
    !>  [getReversed](@ref pm_arrayReverse::getReversed)<br>
    !>  [setReversed](@ref pm_arrayReverse::setReversed)<br>
    !>
    !>  \example{setShuffled}
    !>  \include{lineno} example/pm_arrayShuffle/setShuffled/main.F90
    !>  \compilef{setShuffled}
    !>  \output{setShuffled}
    !>  \include{lineno} example/pm_arrayShuffle/setShuffled/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayShuffle](@ref test_pm_arrayShuffle)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to 2D input objects.<br>
    !>
    !>  \todo
    !>  \phigh
    !>  **Update 2021: This task is now resolved.**<br>
    !>  The current random integer generator uses a simple double precision `real` conversion to `integer` values.<br>
    !>  While this works fairly well for most use cases, it may biased for generating random `integer` of kind \IK4.<br>
    !>  A future remedy should use Bitmask with Rejection as described [here](https://www.pcg-random.org/posts/bounded-rands.html).<br>
    !>  As of 2021, the use of double precision (64-bit) vs. single-precision for random number generation increases the
    !>  computational cost of the algorithms by about three times.<br>
    !>
    !>  \final{setShuffled}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setShuffled

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setShuffledRNGD_D0_SK5(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setShuffledRNGD_D0_SK4(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setShuffledRNGD_D0_SK3(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setShuffledRNGD_D0_SK2(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setShuffledRNGD_D0_SK1(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setShuffledRNGD_D1_SK5(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setShuffledRNGD_D1_SK4(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setShuffledRNGD_D1_SK3(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setShuffledRNGD_D1_SK2(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setShuffledRNGD_D1_SK1(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setShuffledRNGD_D1_IK5(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setShuffledRNGD_D1_IK4(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setShuffledRNGD_D1_IK3(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setShuffledRNGD_D1_IK2(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setShuffledRNGD_D1_IK1(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setShuffledRNGD_D1_LK5(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setShuffledRNGD_D1_LK4(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setShuffledRNGD_D1_LK3(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setShuffledRNGD_D1_LK2(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setShuffledRNGD_D1_LK1(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setShuffledRNGD_D1_CK5(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setShuffledRNGD_D1_CK4(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setShuffledRNGD_D1_CK3(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setShuffledRNGD_D1_CK2(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setShuffledRNGD_D1_CK1(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setShuffledRNGD_D1_RK5(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setShuffledRNGD_D1_RK4(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setShuffledRNGD_D1_RK3(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setShuffledRNGD_D1_RK2(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setShuffledRNGD_D1_RK1(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setShuffledRNGD_D1_PSSK5(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setShuffledRNGD_D1_PSSK4(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setShuffledRNGD_D1_PSSK3(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setShuffledRNGD_D1_PSSK2(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setShuffledRNGD_D1_PSSK1(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setShuffledRNGD_D1_BSSK(array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGD_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        type(css_type)              , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setShuffledRNGF_D0_SK5(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setShuffledRNGF_D0_SK4(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setShuffledRNGF_D0_SK3(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setShuffledRNGF_D0_SK2(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setShuffledRNGF_D0_SK1(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setShuffledRNGF_D1_SK5(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setShuffledRNGF_D1_SK4(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setShuffledRNGF_D1_SK3(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setShuffledRNGF_D1_SK2(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setShuffledRNGF_D1_SK1(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setShuffledRNGF_D1_IK5(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setShuffledRNGF_D1_IK4(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setShuffledRNGF_D1_IK3(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setShuffledRNGF_D1_IK2(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setShuffledRNGF_D1_IK1(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setShuffledRNGF_D1_LK5(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setShuffledRNGF_D1_LK4(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setShuffledRNGF_D1_LK3(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setShuffledRNGF_D1_LK2(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setShuffledRNGF_D1_LK1(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setShuffledRNGF_D1_CK5(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setShuffledRNGF_D1_CK4(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setShuffledRNGF_D1_CK3(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setShuffledRNGF_D1_CK2(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setShuffledRNGF_D1_CK1(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setShuffledRNGF_D1_RK5(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setShuffledRNGF_D1_RK4(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setShuffledRNGF_D1_RK3(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setShuffledRNGF_D1_RK2(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setShuffledRNGF_D1_RK1(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setShuffledRNGF_D1_PSSK5(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setShuffledRNGF_D1_PSSK4(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setShuffledRNGF_D1_PSSK3(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setShuffledRNGF_D1_PSSK2(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setShuffledRNGF_D1_PSSK1(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setShuffledRNGF_D1_BSSK(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGF_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        type(css_type)              , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(rngf_type)             , intent(in)                    :: rng
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setShuffledRNGX_D0_SK5(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setShuffledRNGX_D0_SK4(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setShuffledRNGX_D0_SK3(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setShuffledRNGX_D0_SK2(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setShuffledRNGX_D0_SK1(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setShuffledRNGX_D1_SK5(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setShuffledRNGX_D1_SK4(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setShuffledRNGX_D1_SK3(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setShuffledRNGX_D1_SK2(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setShuffledRNGX_D1_SK1(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setShuffledRNGX_D1_IK5(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setShuffledRNGX_D1_IK4(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setShuffledRNGX_D1_IK3(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setShuffledRNGX_D1_IK2(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setShuffledRNGX_D1_IK1(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setShuffledRNGX_D1_LK5(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setShuffledRNGX_D1_LK4(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setShuffledRNGX_D1_LK3(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setShuffledRNGX_D1_LK2(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setShuffledRNGX_D1_LK1(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setShuffledRNGX_D1_CK5(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setShuffledRNGX_D1_CK4(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setShuffledRNGX_D1_CK3(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setShuffledRNGX_D1_CK2(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setShuffledRNGX_D1_CK1(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setShuffledRNGX_D1_RK5(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setShuffledRNGX_D1_RK4(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setShuffledRNGX_D1_RK3(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setShuffledRNGX_D1_RK2(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setShuffledRNGX_D1_RK1(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setShuffledRNGX_D1_PSSK5(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setShuffledRNGX_D1_PSSK4(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setShuffledRNGX_D1_PSSK3(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setShuffledRNGX_D1_PSSK2(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setShuffledRNGX_D1_PSSK1(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setShuffledRNGX_D1_BSSK(rng, array, count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShuffledRNGX_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        type(css_type)              , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , optional      :: count
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayShuffle ! LCOV_EXCL_LINE