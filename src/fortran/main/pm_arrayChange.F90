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
!>  This module contains procedures and generic interfaces for selecting uniformly-distributed
!>  random choices from a given `character` or `integer` range.<br>
!>
!>  \details
!>  The functionality of the procedures of this module are identical to those of [pm_arrayShuffle](@ref pm_arrayShuffle), except that
!>  <ol>
!>      <li>    the procedure within this module compute a shuffling of a (subset of a) range specified by lower and upper limits.<br>
!>      <li>    the procedure within this module additionally can compute a random sequence from a range **with replacement (non-unique)**.<br>
!>  </ol>
!>
!>  \remark
!>  The word **change** in naming of this generic interface and its output stands for <b>CH</b>oice in r<b>ANGE</b>.<br>
!>
!>  \note
!>  The case of length of the output sequence of **unique** elements being equal to the length of the specified range `[start, stop]`
!>  is equivalent to a [random shuffling](@ref pm_arrayShuffle) of all values within the specified range.<br>
!>
!>  \see
!>  [pm_arrayRemap](@ref pm_arrayRemap)<br>
!>  [pm_arrayChange](@ref pm_arrayChange)<br>
!>  [pm_arrayChoice](@ref pm_arrayChoice)<br>
!>  [pm_arrayShuffle](@ref pm_arrayShuffle)<br>
!>  [pm_distUnif](@ref pm_distUnif)<br>
!>
!>  \test
!>  [test_pm_arrayChange](@ref test_pm_arrayChange)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayChange

    use pm_kind, only: SK, IK, LK
    use pm_distUnif, only: rngf_type, xoshiro256ssw_type

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_arrayChange"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a randomly-uniform selected sequence of elements from
    !>  a scalar `character` or `integer` range specified via the input limits `[start, stop]`.<br>
    !>
    !>  \details
    !>  Should the sequence elements be unique, the length of the output sequence must be smaller than or equal to the specified range as `[start, stop]`.<br>
    !>
    !>  \param[in]  count   :   The input non-negative scalar `integer` of default kind \IK, representing the length
    !>                          of the output sequence of choices in the range specified by `[start, stop]`.<br>
    !>                          If `unique = .true.`, the condition `count <= len(getRange(start, stop, step))` or `count <= size(getRange(start, stop, step))`
    !>                          must hold to ensure uniqueness of the elements of the output sequence.<br>
    !>  \param[in]  start   :   See the documentation of the corresponding argument of [setChange(rng, change, start, stop, step, unique)](@ref pm_arrayChange::setChange).<br>
    !>  \param[in]  stop    :   See the documentation of the corresponding argument of [setChange(rng, change, start, stop, step, unique)](@ref pm_arrayChange::setChange).<br>
    !>  \param[in]  step    :   See the documentation of the corresponding argument of [setChange(rng, change, start, stop, step, unique)](@ref pm_arrayChange::setChange).<br>
    !>  \param[in]  unique  :   See the documentation of the corresponding argument of [setChange(rng, change, start, stop, step, unique)](@ref pm_arrayChange::setChange).<br>
    !>                          (**optional**, default = `.false._LK`)
    !>
    !>  \return
    !>  `change`            :   The output vector of size `count`,
    !>                          <ol>
    !>                              <li>    type `integer` of kind \IKALL,<br>
    !>                          </ol>
    !>                          or a scalar of,
    !>                          <ol>
    !>                              <li>    type `character` of kind \SKALL of length `count`, <br>
    !>                          </ol>
    !>                          whose elements are uniformly randomly selected from within the range [getRange(start, stop, step)](@ref pm_arrayRange::getRange).<br>
    !>
    !>  \interface{getChange}
    !>  \code{.F90}
    !>
    !>      use pm_arrayChange, only: getChange
    !>
    !>      change = getChange(count, start, stop, step, unique = unique)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  All warnings and conditions associated with [getRange(start, stop, step)](@ref pm_arrayRange::getRange)
    !>  and [setChoice(start, stop, step)](@ref pm_arrayChoice::setChoice) equally hold for the procedures of this generic interface.<br>
    !>  The length of the resulting sequence from the specified range `(start, stop, step)` must not be zero.<br>
    !>  \vericons
    !>
    !>  \impure
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
    !>  \example{getChange}
    !>  \include{lineno} example/pm_arrayChange/getChange/main.F90
    !>  \compilef{getChange}
    !>  \output{getChange}
    !>  \include{lineno} example/pm_arrayChange/getChange/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayChange](@ref test_pm_arrayChange)
    !>
    !>  \todo
    !>  The input argument `step` can be made `optional`.<br>
    !>  Currently, it is a required argument to avoid testing complexities in the case of input real-valued ranges.<br>
    !>
    !>  \final{getChange}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    interface getChange

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getChangeUnifRNGD_SK5(count, start, stop, step, unique) result(change)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChangeUnifRNGD_SK5
#endif
        use pm_kind, only: SKC => SK5
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: step
        integer(IK)                 , intent(in)                    :: count
        character(1,SKC)            , intent(in)                    :: start, stop
        character(count,SKC)                                        :: change
    end function
#endif

#if SK4_ENABLED
    module function getChangeUnifRNGD_SK4(count, start, stop, step, unique) result(change)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChangeUnifRNGD_SK4
#endif
        use pm_kind, only: SKC => SK4
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: step
        integer(IK)                 , intent(in)                    :: count
        character(1,SKC)            , intent(in)                    :: start, stop
        character(count,SKC)                                        :: change
    end function
#endif

#if SK3_ENABLED
    module function getChangeUnifRNGD_SK3(count, start, stop, step, unique) result(change)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChangeUnifRNGD_SK3
#endif
        use pm_kind, only: SKC => SK3
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: step
        integer(IK)                 , intent(in)                    :: count
        character(1,SKC)            , intent(in)                    :: start, stop
        character(count,SKC)                                        :: change
    end function
#endif

#if SK2_ENABLED
    module function getChangeUnifRNGD_SK2(count, start, stop, step, unique) result(change)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChangeUnifRNGD_SK2
#endif
        use pm_kind, only: SKC => SK2
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: step
        integer(IK)                 , intent(in)                    :: count
        character(1,SKC)            , intent(in)                    :: start, stop
        character(count,SKC)                                        :: change
    end function
#endif

#if SK1_ENABLED
    module function getChangeUnifRNGD_SK1(count, start, stop, step, unique) result(change)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChangeUnifRNGD_SK1
#endif
        use pm_kind, only: SKC => SK1
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: step
        integer(IK)                 , intent(in)                    :: count
        character(1,SKC)            , intent(in)                    :: start, stop
        character(count,SKC)                                        :: change
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getChangeUnifRNGD_IK5(count, start, stop, step, unique) result(change)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChangeUnifRNGD_IK5
#endif
        use pm_kind, only: IKC => IK5
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: count
        integer(IKC)                , intent(in)                    :: start, stop, step
        integer(IKC)                                                :: change(count)
    end function
#endif

#if IK4_ENABLED
    module function getChangeUnifRNGD_IK4(count, start, stop, step, unique) result(change)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChangeUnifRNGD_IK4
#endif
        use pm_kind, only: IKC => IK4
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: count
        integer(IKC)                , intent(in)                    :: start, stop, step
        integer(IKC)                                                :: change(count)
    end function
#endif

#if IK3_ENABLED
    module function getChangeUnifRNGD_IK3(count, start, stop, step, unique) result(change)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChangeUnifRNGD_IK3
#endif
        use pm_kind, only: IKC => IK3
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: count
        integer(IKC)                , intent(in)                    :: start, stop, step
        integer(IKC)                                                :: change(count)
    end function
#endif

#if IK2_ENABLED
    module function getChangeUnifRNGD_IK2(count, start, stop, step, unique) result(change)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChangeUnifRNGD_IK2
#endif
        use pm_kind, only: IKC => IK2
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: count
        integer(IKC)                , intent(in)                    :: start, stop, step
        integer(IKC)                                                :: change(count)
    end function
#endif

#if IK1_ENABLED
    module function getChangeUnifRNGD_IK1(count, start, stop, step, unique) result(change)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChangeUnifRNGD_IK1
#endif
        use pm_kind, only: IKC => IK1
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: count
        integer(IKC)                , intent(in)                    :: start, stop, step
        integer(IKC)                                                :: change(count)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getChangeUnifRNGD_RK5(count, start, stop, step, unique) result(change)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChangeUnifRNGD_RK5
#endif
        use pm_kind, only: RKC => RK5
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: count
        real(RKC)                   , intent(in)                    :: start, stop, step
        real(RKC)                                                   :: change(count)
    end function
#endif

#if RK4_ENABLED
    module function getChangeUnifRNGD_RK4(count, start, stop, step, unique) result(change)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChangeUnifRNGD_RK4
#endif
        use pm_kind, only: RKC => RK4
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: count
        real(RKC)                   , intent(in)                    :: start, stop, step
        real(RKC)                                                   :: change(count)
    end function
#endif

#if RK3_ENABLED
    module function getChangeUnifRNGD_RK3(count, start, stop, step, unique) result(change)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChangeUnifRNGD_RK3
#endif
        use pm_kind, only: RKC => RK3
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: count
        real(RKC)                   , intent(in)                    :: start, stop, step
        real(RKC)                                                   :: change(count)
    end function
#endif

#if RK2_ENABLED
    module function getChangeUnifRNGD_RK2(count, start, stop, step, unique) result(change)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChangeUnifRNGD_RK2
#endif
        use pm_kind, only: RKC => RK2
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: count
        real(RKC)                   , intent(in)                    :: start, stop, step
        real(RKC)                                                   :: change(count)
    end function
#endif

#if RK1_ENABLED
    module function getChangeUnifRNGD_RK1(count, start, stop, step, unique) result(change)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChangeUnifRNGD_RK1
#endif
        use pm_kind, only: RKC => RK1
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: count
        real(RKC)                   , intent(in)                    :: start, stop, step
        real(RKC)                                                   :: change(count)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a randomly-uniform uniquely-selected sequence of elements from
    !>  a scalar `character` or `integer` range specified via the input limits `[start, stop]`.<br>
    !>
    !>  \details
    !>  The length of the output sequence must be smaller than or equal to the specified range as `[start, stop]` to ensure uniqueness of elements.<br>
    !>  The case of length of the output sequence being equal to the length of the specified range `[start, stop]`
    !>  is equivalent to a random shuffling of all values within the specified range.<br>
    !>
    !>  \param[inout]   rng     :   The input/output scalar that can be an object of,
    !>                              <ol>
    !>                                  <li>    type [rngf_type](@ref pm_distUnif::rngf_type),
    !>                                          implying the use of intrinsic Fortran uniform RNG.<br>
    !>                                  <li>    type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type),
    !>                                          implying the use of [xoshiro256**](https://prng.di.unimi.it/) uniform RNG.<br>
    !>                              </ol>
    !>  \param[out]     change  :   The output vector of shape `(:)`,
    !>                              <ol>
    !>                                  <li>    type `integer` of kind \IKALL,<br>
    !>                              </ol>
    !>                              or a scalar of,
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL, <br>
    !>                              </ol>
    !>                              whose elements are uniformly randomly selected from within the specified range as `[start, stop]`.<br>
    !>                              The length/size of `change` must be smaller than or equal to the length of the range specified by `[start, stop]`.<br>
    !>  \param[in]      start   :   The input scalar of the same type and kind as the output `change` containing the lower limit of
    !>                              the range from which element(s) will be selected randomly and returned in the output `change`.<br>
    !>                              If `start` is of type `character`, then it must be of length type parameter `1`.<br>
    !>  \param[in]      step    :   The **non-zero** input scalar of type `integer` representing the size of the
    !>                              interval between adjacent values in the output `change`.<br>
    !>                              <ol>
    !>                                  <li>    If `change` is of type `character`, then `step` must have the default `integer` kind \IK.<br>
    !>                                  <li>    If `change` is of type `integer`, then `step` must have the same type kind parameter as the output `change`.<br>
    !>                              </ol>
    !>                              If `stop` is of type `character`, then it must be of length type parameter `1`.<br>
    !>  \param[in]      unique  :   The input scalar `logical` of default kind \LK.<br>
    !>                              If `.true.`, the elements of the output `change` will be **uniquely** selected from the input `array`.<br>
    !>                              If `unique = .true.`, the condition `count <= len(getRange(start, stop, step))` or `count <= size(getRange(start, stop, step))`
    !>                              must hold to ensure uniqueness of the elements of the output sequence.<br>
    !>                              (**optional**, default = `.false._LK`)
    !>
    !>
    !>  \interface{setChange}
    !>  \code{.F90}
    !>
    !>      use pm_arrayChange, only: setChange
    !>
    !>      change = setChange(rng, change, start, stop, step, unique = unique)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  All warnings and conditions associated with [getRange(start, stop, step)](@ref pm_arrayRange::getRange)
    !>  and [setChoice(start, stop, step)](@ref pm_arrayChoice::setChoice) equally hold for the procedures of this generic interface.<br>
    !>  The length of the resulting sequence from the specified range `(start, stop, step)` must not be zero.<br>
    !>  \vericons
    !>
    !>  \impure
    !>  The procedures of this generic interface are `pure` when the library is built with preprocessor macro `CHECK_ENABLED=0`
    !>  and the input argument `rng` is set to an object of type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type).<br>
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
    !>  \example{setChange}
    !>  \include{lineno} example/pm_arrayChange/setChange/main.F90
    !>  \compilef{setChange}
    !>  \output{setChange}
    !>  \include{lineno} example/pm_arrayChange/setChange/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayChange](@ref test_pm_arrayChange)
    !>
    !>  \final{setChange}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    interface setChange

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setChangeUnifRNGF_SK5(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGF_SK5
#endif
        use pm_kind, only: SKC => SK5
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: step
        type(rngf_type)             , intent(in)                    :: rng
        character(1,SKC)            , intent(in)                    :: start, stop
        character(*,SKC)            , intent(out)                   :: change
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setChangeUnifRNGF_SK4(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGF_SK4
#endif
        use pm_kind, only: SKC => SK4
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: step
        type(rngf_type)             , intent(in)                    :: rng
        character(1,SKC)            , intent(in)                    :: start, stop
        character(*,SKC)            , intent(out)                   :: change
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setChangeUnifRNGF_SK3(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGF_SK3
#endif
        use pm_kind, only: SKC => SK3
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: step
        type(rngf_type)             , intent(in)                    :: rng
        character(1,SKC)            , intent(in)                    :: start, stop
        character(*,SKC)            , intent(out)                   :: change
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setChangeUnifRNGF_SK2(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGF_SK2
#endif
        use pm_kind, only: SKC => SK2
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: step
        type(rngf_type)             , intent(in)                    :: rng
        character(1,SKC)            , intent(in)                    :: start, stop
        character(*,SKC)            , intent(out)                   :: change
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setChangeUnifRNGF_SK1(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGF_SK1
#endif
        use pm_kind, only: SKC => SK1
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: step
        type(rngf_type)             , intent(in)                    :: rng
        character(1,SKC)            , intent(in)                    :: start, stop
        character(*,SKC)            , intent(out)                   :: change
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setChangeUnifRNGF_IK5(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGF_IK5
#endif
        use pm_kind, only: IKC => IK5
        logical(LK)                 , intent(in)    , optional      :: unique
        type(rngf_type)             , intent(in)                    :: rng
        integer(IKC)                , intent(in)                    :: start, stop, step
        integer(IKC)                , intent(out)   , contiguous    :: change(:)
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setChangeUnifRNGF_IK4(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGF_IK4
#endif
        use pm_kind, only: IKC => IK4
        logical(LK)                 , intent(in)    , optional      :: unique
        type(rngf_type)             , intent(in)                    :: rng
        integer(IKC)                , intent(in)                    :: start, stop, step
        integer(IKC)                , intent(out)   , contiguous    :: change(:)
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setChangeUnifRNGF_IK3(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGF_IK3
#endif
        use pm_kind, only: IKC => IK3
        logical(LK)                 , intent(in)    , optional      :: unique
        type(rngf_type)             , intent(in)                    :: rng
        integer(IKC)                , intent(in)                    :: start, stop, step
        integer(IKC)                , intent(out)   , contiguous    :: change(:)
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setChangeUnifRNGF_IK2(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGF_IK2
#endif
        use pm_kind, only: IKC => IK2
        logical(LK)                 , intent(in)    , optional      :: unique
        type(rngf_type)             , intent(in)                    :: rng
        integer(IKC)                , intent(in)                    :: start, stop, step
        integer(IKC)                , intent(out)   , contiguous    :: change(:)
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setChangeUnifRNGF_IK1(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGF_IK1
#endif
        use pm_kind, only: IKC => IK1
        logical(LK)                 , intent(in)    , optional      :: unique
        type(rngf_type)             , intent(in)                    :: rng
        integer(IKC)                , intent(in)                    :: start, stop, step
        integer(IKC)                , intent(out)   , contiguous    :: change(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setChangeUnifRNGF_RK5(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGF_RK5
#endif
        use pm_kind, only: RKC => RK5
        logical(LK)                 , intent(in)    , optional      :: unique
        type(rngf_type)             , intent(in)                    :: rng
        real(RKC)                   , intent(in)                    :: start, stop, step
        real(RKC)                   , intent(out)   , contiguous    :: change(:)
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setChangeUnifRNGF_RK4(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGF_RK4
#endif
        use pm_kind, only: RKC => RK4
        logical(LK)                 , intent(in)    , optional      :: unique
        type(rngf_type)             , intent(in)                    :: rng
        real(RKC)                   , intent(in)                    :: start, stop, step
        real(RKC)                   , intent(out)   , contiguous    :: change(:)
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setChangeUnifRNGF_RK3(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGF_RK3
#endif
        use pm_kind, only: RKC => RK3
        logical(LK)                 , intent(in)    , optional      :: unique
        type(rngf_type)             , intent(in)                    :: rng
        real(RKC)                   , intent(in)                    :: start, stop, step
        real(RKC)                   , intent(out)   , contiguous    :: change(:)
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setChangeUnifRNGF_RK2(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGF_RK2
#endif
        use pm_kind, only: RKC => RK2
        logical(LK)                 , intent(in)    , optional      :: unique
        type(rngf_type)             , intent(in)                    :: rng
        real(RKC)                   , intent(in)                    :: start, stop, step
        real(RKC)                   , intent(out)   , contiguous    :: change(:)
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setChangeUnifRNGF_RK1(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGF_RK1
#endif
        use pm_kind, only: RKC => RK1
        logical(LK)                 , intent(in)    , optional      :: unique
        type(rngf_type)             , intent(in)                    :: rng
        real(RKC)                   , intent(in)                    :: start, stop, step
        real(RKC)                   , intent(out)   , contiguous    :: change(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setChangeUnifRNGX_SK5(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGX_SK5
#endif
        use pm_kind, only: SKC => SK5
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: step
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        character(1,SKC)            , intent(in)                    :: start, stop
        character(*,SKC)            , intent(out)                   :: change
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setChangeUnifRNGX_SK4(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGX_SK4
#endif
        use pm_kind, only: SKC => SK4
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: step
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        character(1,SKC)            , intent(in)                    :: start, stop
        character(*,SKC)            , intent(out)                   :: change
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setChangeUnifRNGX_SK3(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGX_SK3
#endif
        use pm_kind, only: SKC => SK3
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: step
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        character(1,SKC)            , intent(in)                    :: start, stop
        character(*,SKC)            , intent(out)                   :: change
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setChangeUnifRNGX_SK2(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGX_SK2
#endif
        use pm_kind, only: SKC => SK2
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: step
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        character(1,SKC)            , intent(in)                    :: start, stop
        character(*,SKC)            , intent(out)                   :: change
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setChangeUnifRNGX_SK1(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGX_SK1
#endif
        use pm_kind, only: SKC => SK1
        logical(LK)                 , intent(in)    , optional      :: unique
        integer(IK)                 , intent(in)                    :: step
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        character(1,SKC)            , intent(in)                    :: start, stop
        character(*,SKC)            , intent(out)                   :: change
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setChangeUnifRNGX_IK5(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGX_IK5
#endif
        use pm_kind, only: IKC => IK5
        logical(LK)                 , intent(in)    , optional      :: unique
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        integer(IKC)                , intent(in)                    :: start, stop, step
        integer(IKC)                , intent(out)   , contiguous    :: change(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setChangeUnifRNGX_IK4(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGX_IK4
#endif
        use pm_kind, only: IKC => IK4
        logical(LK)                 , intent(in)    , optional      :: unique
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        integer(IKC)                , intent(in)                    :: start, stop, step
        integer(IKC)                , intent(out)   , contiguous    :: change(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setChangeUnifRNGX_IK3(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGX_IK3
#endif
        use pm_kind, only: IKC => IK3
        logical(LK)                 , intent(in)    , optional      :: unique
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        integer(IKC)                , intent(in)                    :: start, stop, step
        integer(IKC)                , intent(out)   , contiguous    :: change(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setChangeUnifRNGX_IK2(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGX_IK2
#endif
        use pm_kind, only: IKC => IK2
        logical(LK)                 , intent(in)    , optional      :: unique
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        integer(IKC)                , intent(in)                    :: start, stop, step
        integer(IKC)                , intent(out)   , contiguous    :: change(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setChangeUnifRNGX_IK1(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGX_IK1
#endif
        use pm_kind, only: IKC => IK1
        logical(LK)                 , intent(in)    , optional      :: unique
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        integer(IKC)                , intent(in)                    :: start, stop, step
        integer(IKC)                , intent(out)   , contiguous    :: change(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setChangeUnifRNGX_RK5(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGX_RK5
#endif
        use pm_kind, only: RKC => RK5
        logical(LK)                 , intent(in)    , optional      :: unique
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKC)                   , intent(in)                    :: start, stop, step
        real(RKC)                   , intent(out)   , contiguous    :: change(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setChangeUnifRNGX_RK4(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGX_RK4
#endif
        use pm_kind, only: RKC => RK4
        logical(LK)                 , intent(in)    , optional      :: unique
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKC)                   , intent(in)                    :: start, stop, step
        real(RKC)                   , intent(out)   , contiguous    :: change(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setChangeUnifRNGX_RK3(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGX_RK3
#endif
        use pm_kind, only: RKC => RK3
        logical(LK)                 , intent(in)    , optional      :: unique
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKC)                   , intent(in)                    :: start, stop, step
        real(RKC)                   , intent(out)   , contiguous    :: change(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setChangeUnifRNGX_RK2(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGX_RK2
#endif
        use pm_kind, only: RKC => RK2
        logical(LK)                 , intent(in)    , optional      :: unique
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKC)                   , intent(in)                    :: start, stop, step
        real(RKC)                   , intent(out)   , contiguous    :: change(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setChangeUnifRNGX_RK1(rng, change, start, stop, step, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChangeUnifRNGX_RK1
#endif
        use pm_kind, only: RKC => RK1
        logical(LK)                 , intent(in)    , optional      :: unique
        type(xoshiro256ssw_type)    , intent(inout)                 :: rng
        real(RKC)                   , intent(in)                    :: start, stop, step
        real(RKC)                   , intent(out)   , contiguous    :: change(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayChange ! LCOV_EXCL_LINE