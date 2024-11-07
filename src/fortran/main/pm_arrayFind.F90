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
!>  This module contains procedures and generic interfaces for finding locations of a pattern in arrays of various types at the specified instances of occurrence of pattern.
!>
!>  \benchmarks
!>
!>  \benchmark{setLoc-scalarPattern_vs_vectorPattern, The runtime performance of [setLoc](@ref pm_arrayFind::setLoc) for scalar vs. vector input `pattern` argument.}
!>  \include{lineno} benchmark/pm_arrayFind/setLoc-scalarPattern_vs_vectorPattern/main.F90
!>  \compilefb{setLoc-scalarPattern_vs_vectorPattern}
!>  \postprocb{setLoc-scalarPattern_vs_vectorPattern}
!>  \include{lineno} benchmark/pm_arrayFind/setLoc-scalarPattern_vs_vectorPattern/main.py
!>  \visb{setLoc-scalarPattern_vs_vectorPattern}
!>  \image html benchmark/pm_arrayFind/setLoc-scalarPattern_vs_vectorPattern/benchmark.setLoc-scalarPattern_vs_vectorPattern.runtime.png width=1000
!>  \image html benchmark/pm_arrayFind/setLoc-scalarPattern_vs_vectorPattern/benchmark.setLoc-scalarPattern_vs_vectorPattern.runtime.ratio.png width=1000
!>  \moralb{setLoc-scalarPattern_vs_vectorPattern}
!>      -#  The procedures under the generic interface [setLoc](@ref pm_arrayFind::setLoc) take both scalar and vector `pattern` arguments.<br>
!>          As evidenced by the above benchmark, when the input `pattern` is vector of length `1`, it is much faster, possibly **up to 4X**,
!>          to pass `pattern` as a scalar instead of a whole array of length `1`.<br>
!>          Note that this benchmark is likely irrelevant to finding substrings in Fortran strings.<br>
!>
!>  \benchmark{getLoc_vs_setLoc, The runtime performance of [getLoc](@ref pm_arrayFind::getLoc) vs. [setLoc](@ref pm_arrayFind::setLoc)}
!>  \include{lineno} benchmark/pm_arrayFind/getLoc_vs_setLoc/main.F90
!>  \compilefb{getLoc_vs_setLoc}
!>  \postprocb{getLoc_vs_setLoc}
!>  \include{lineno} benchmark/pm_arrayFind/getLoc_vs_setLoc/main.py
!>  \visb{getLoc_vs_setLoc}
!>  \image html benchmark/pm_arrayFind/getLoc_vs_setLoc/benchmark.getLoc_vs_setLoc.runtime.png width=1000
!>  \image html benchmark/pm_arrayFind/getLoc_vs_setLoc/benchmark.getLoc_vs_setLoc.runtime.ratio.png width=1000
!>  \moralb{getLoc_vs_setLoc}
!>      -#  The procedures under the generic interface [getLoc](@ref pm_arrayFind::getLoc) are functions while
!>          the procedures under the generic interface [setLoc](@ref pm_arrayFind::setLoc) are subroutines.<br>
!>          From the benchmark results, it appears that the functional interface performs slightly less efficiently than
!>          the subroutine interface when the input `array` size is small.<br>
!>          Nevertheless, the difference appears to be marginal and insignificant in most practical situations.<br>
!>      -#  Note that this benchmark considers the worst-case scenario where all elements
!>          of the input `array` match the input `pattern` and must be therefore indexed.<br>
!>
!>  \benchmark{setLoc-sortedPositive_vs_defaultOptions, The runtime performance of [setLoc](@ref pm_arrayFind::setLoc) with and without `sorted = .true.` and `positive = .true.` optional input arguments.}
!>  \include{lineno} benchmark/pm_arrayFind/setLoc-sortedPositive_vs_defaultOptions/main.F90
!>  \compilefb{setLoc-sortedPositive_vs_defaultOptions}
!>  \postprocb{setLoc-sortedPositive_vs_defaultOptions}
!>  \include{lineno} benchmark/pm_arrayFind/setLoc-sortedPositive_vs_defaultOptions/main.py
!>  \visb{setLoc-sortedPositive_vs_defaultOptions}
!>  \image html benchmark/pm_arrayFind/setLoc-sortedPositive_vs_defaultOptions/benchmark.setLoc-sortedPositive_vs_defaultOptions.runtime.png width=1000
!>  \image html benchmark/pm_arrayFind/setLoc-sortedPositive_vs_defaultOptions/benchmark.setLoc-sortedPositive_vs_defaultOptions.runtime.ratio.png width=1000
!>  \moralb{setLoc-sortedPositive_vs_defaultOptions}
!>      -#  The procedures under the generic interface [setLoc](@ref pm_arrayFind::setLoc) take two optional `sorted` and `positive` input arguments
!>          when the input array of instances `instance` is also specified. Setting `sorted = .true._LK, positive = .true.` guarantees to the procedures
!>          that all elements of `instance` are sorted in **ascending-order** and that no negative or zero value is present in `instance`. Hence,
!>          the maximum value in `instance` corresponds to `instance(size(instance))`, saving the procedures the time to find the maximum value
!>          of the input `instance`.<br>
!>      -#  As evidenced by the above benchmark, when `sorted = .true._LK, positive = .true.` can be guaranteed by the user, the runtime performance
!>          of the procedures could possibly increase **2-3 times** compared to when the default values have to be assumed, corresponding to the
!>          worst case scenario (`sorted = .false._LK, positive = .false.`).<br>
!>      -#  The results of this benchmark equally hold for the functions under the generic interface [getLoc](@ref pm_arrayFind::getLoc).<br>
!>
!>  \benchmark{setLoc-iseqArg_vs_default, The runtime performance of [setLoc](@ref pm_arrayFind::setLoc) with and without `iseq` external optional input argument.}
!>  \include{lineno} benchmark/pm_arrayFind/setLoc-iseqArg_vs_default/main.F90
!>  \compilefb{setLoc-iseqArg_vs_default}
!>  \postprocb{setLoc-iseqArg_vs_default}
!>  \include{lineno} benchmark/pm_arrayFind/setLoc-iseqArg_vs_default/main.py
!>  \visb{setLoc-iseqArg_vs_default}
!>  \image html benchmark/pm_arrayFind/setLoc-iseqArg_vs_default/benchmark.setLoc-iseqArg_vs_default.runtime.png width=1000
!>  \image html benchmark/pm_arrayFind/setLoc-iseqArg_vs_default/benchmark.setLoc-iseqArg_vs_default.runtime.ratio.png width=1000
!>  \moralb{setLoc-iseqArg_vs_default}
!>      -#  The procedures under the generic interface [setLoc](@ref pm_arrayFind::setLoc) take an optional `iseq()` external-function input argument
!>          which can perform custom user-defined equivalence checks between the input `pattern` and `array` segments. By default, the equivalence check
!>          is equality (or equivalence for `logical` values).<br>
!>      -#  The presence of an external `iseq` argument does appear to deteriorate the runtime performance of [setLoc](@ref pm_arrayFind::setLoc).<br>
!>          However, the exact amount of this penalty appears to significantly depend on the computing platform and the ability of
!>          the compiler to inline the external function.<br>
!>      -#  The results of this benchmark equally hold for the functions under the generic interface [getLoc](@ref pm_arrayFind::getLoc).<br>
!>
!>  \test
!>  [test_pm_arrayFind](@ref test_pm_arrayFind)
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX
module pm_arrayFind

    use pm_kind, only: SK, IK, LK
    use pm_array, only: border_type, discrete, discrete_type

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_arrayFind"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the number of occurrences of the input `pattern` in the input `array` optionally subject to user-specified memory blindness.
    !>
    !>  \details
    !>  The instances of `pattern` are found via linear search.<br>
    !>  Therefore, the procedures under this generic interface have a worst-case complexity of `O(size(array))`.
    !>
    !>  \param[in]  array       :   The input `contiguous` array of rank `1` of either<br>
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL, or <br>
    !>                                  <li>    type `integer` of kind \IKALL, or <br>
    !>                                  <li>    type `logical` of kind \LKALL, or <br>
    !>                                  <li>    type `complex` of kind \CKALL, or <br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ol>
    !>                              or
    !>                              <ol>
    !>                                  <li>    a scalar assumed-length `character` of kind \SKALL,<br>
    !>                              </ol>
    !>                              within which the starting indices of the requested instances of `pattern` is to be found.<br>
    !>  \param[in]  pattern     :   The input object of the same or lower rank than the input `array`, and of the same type
    !>                              and kind as `array` containing the pattern that must be found within the input `array`.<br>
    !>  \param[in]  border      :   The input scalar constant object that can be any of the following,<br>
    !>                              <ol>
    !>                                  <li>    the constant [discrete](@ref pm_array::discrete) or an object of type [discrete_type](@ref pm_array::discrete_type),
    !>                                          implying that only non-overlapping and non-adjacent pattern locations must be identified.<br>
    !>                                          For example, if `pattern` is a blank character, and multiple adjacent blanks appear in `array`,
    !>                                          the only the location of the first blank is identified and the rest are ignored until a blank
    !>                                          reappears separately from the first group of contiguous blanks.<br>
    !>                                          Note that border adjacency and overlapping is measured by the input `blindness` and not by the length of `pattern`.<br>
    !>                                          This option is extremely useful for identifying and parsing separators that can be repeated in text file records,
    !>                                          for example, the white-space (blank) character in list-directed Fortran IO.<br>
    !>                              </ol>
    !>                              (**optional**. The default behavior does not recognize any borders for `pattern`.)
    !>  \param      iseq        :   The `external` user-specified function that takes two input **explicit-shape** arguments of the same type
    !>                              and kind as the input `array` and possibly, also the length of the arguments as the third argument, if the arguments are array-valued.<br>
    !>                              It returns a scalar `logical` of default kind \LK that is `.true.` if all elements of the two input arguments are equivalent (e.g., equal)
    !>                              according to the user-defined criterion, otherwise, it is `.false.`.<br>
    !>                              If `pattern` is an array of rank `1`, then the last argument to `iseq` is the length of the input `pattern`, preceded by a
    !>                              segment of `array` and `pattern` as the first and second arguments, whose lengths are given by the third argument `lenPattern`.<br>
    !>                              The following illustrates the generic interface of `iseq` where `pattern` is array-valued,
    !>                              \code{.F90}
    !>                                  function iseq(Segment, pattern, lenPattern) result(equivalent)
    !>                                      use pm_kind, only: IK, LK
    !>                                      integer(IK) , intent(in)    :: lenPattern
    !>                                      TYPE(KIND)  , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                      logical(LK)                 :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,<br>
    !>                              \code{.F90}
    !>                                  character(*, SK), intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  integer(IK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  logical(LK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  complex(CK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  real(RK)        , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              The following illustrates the generic interface of `iseq` where `pattern` is scalar-valued (**including Fortran scalar strings**),<br>
    !>                              \code{.F90}
    !>                                  function iseq(segment, pattern) result(equivalent)
    !>                                      use pm_kind, only: LK
    !>                                      TYPE(KIND)  , intent(in)    :: segment, pattern
    !>                                      logical(LK)                 :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,<br>
    !>                              \code{.F90}
    !>                                  use pm_kind, only: SK, IK, LK, CK, RK
    !>                                  character(*, SK), intent(in)    :: segment, pattern
    !>                                  integer(IK)     , intent(in)    :: segment, pattern
    !>                                  logical(LK)     , intent(in)    :: segment, pattern
    !>                                  complex(CK)     , intent(in)    :: segment, pattern
    !>                                  real(RK)        , intent(in)    :: segment, pattern
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              This user-defined equivalence check is extremely useful where a user-defined equivalence test other than exact equality
    !>                              or identity is needed, for example, when the array segments should match the input `pattern` only within a given threshold or,
    !>                              when the case-sensitivity in character comparisons do not matter.<br>
    !>                              In such cases, user can define a custom equivalence criterion within the user-defined external function `iseq` to achieve the goal.<br>
    !>                              (**optional**, the default equivalence operator is `.eqv.` if the input `array` is `logical`, otherwise `==`.)
    !>  \param[in]  blindness   :   The input **positive** `integer` of default kind \IK representing the length of the segment of `array` that should be ignored after finding an instance of the input `pattern` in the array.<br>
    !>                              Setting `blindness = len(pattern)` (for assumed-length character `pattern`) or `blindness = size(pattern)` (for other types of array-valued `pattern`)
    !>                              will lead to a search for exclusive non-overlapping instances of `pattern` in the input `array`.<br>
    !>                              See the examples below for more illustration of the utility of this input argument.<br>
    !>                              (**optional**, default = `1_IK`)
    !>
    !>  \return
    !>  `count`                 :   The output non-negative scalar of type `integer` of default kind \IK
    !>                              containing the number of times the input pattern appears in the input `array`.
    !>
    !>  \interface{getCountLoc}
    !>  \code{.F90}
    !>
    !>      use pm_arrayFind, only: getCountLoc
    !>
    !>      ! `array` and `pattern` are scalar strings.
    !>
    !>      count = getCountLoc(array, pattern, blindness = blindness)
    !>      count = getCountLoc(array, pattern, iseq, blindness = blindness)
    !>
    !>      ! `pattern` is a scalar of the same type and kind as `array`.
    !>
    !>      count = getCountLoc(array(:), pattern, blindness = blindness)
    !>      count = getCountLoc(array(:), pattern, iseq, blindness = blindness)
    !>
    !>      ! `array` and `pattern` are both vectors.
    !>
    !>      count = getCountLoc(array(:), pattern(:), blindness = blindness)
    !>      count = getCountLoc(array(:), pattern(:), iseq, blindness = blindness)
    !>
    !>      ! `array` and `pattern` are scalar strings.
    !>
    !>      count = getCountLoc(array, pattern, border, blindness = blindness)
    !>      count = getCountLoc(array, pattern, border, iseq, blindness = blindness)
    !>
    !>      ! `pattern` is a scalar of the same type and kind as `array`.
    !>
    !>      count = getCountLoc(array(:), pattern, border, blindness = blindness)
    !>      count = getCountLoc(array(:), pattern, border, iseq, blindness = blindness)
    !>
    !>      ! `array` and `pattern` are both vectors.
    !>
    !>      count = getCountLoc(array(:), pattern(:), border, blindness = blindness)
    !>      count = getCountLoc(array(:), pattern(:), border, iseq, blindness = blindness)
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < blindness` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warning
    !>  The procedures under this generic interface are `impure` when the user-specified `external` procedure `iseq` is specified as input argument.<br>
    !>
    !>  \warning
    !>  Note that in Fortran, trailing blanks are ignored in character comparison, that is, `"Fortran" == "Fortran "` yields `.true.`.<br>
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [setLoc](@ref pm_arrayFind::setLoc)<br>
    !>  [getBin](@ref pm_arraySearch::getBin)<br>
    !>  [setReplaced](@ref pm_arrayReplace::setReplaced)<br>
    !>  [getReplaced](@ref pm_arrayReplace::getReplaced)<br>
    !>  [setInserted](@ref pm_arrayInsert::setInserted)<br>
    !>  [setSplit](@ref pm_arraySplit::setSplit)<br>
    !>
    !>  \example{getCountLoc}
    !>  \include{lineno} example/pm_arrayFind/getCountLoc/main.F90
    !>  \compilef{getCountLoc}
    !>  \output{getCountLoc}
    !>  \include{lineno} example/pm_arrayFind/getCountLoc/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayFind](@ref test_pm_arrayFind)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.2.0}, \gfortran{10-12}
    !>  \desc
    !>  The Intel Fortran compiler Classic 2021.2.0 has a bug for the following interface definition
    !>  \code{.F90}
    !>      character(len(array),SK), allocatable   :: count(:)
    !>  \endcode
    !>  leading to an **internal compiler error**.<br>
    !>  For now, the remedy seems to be to redefine the interface as,<br>
    !>  \code{.F90}
    !>      character(:, SK), allocatable   :: count(:)
    !>  \endcode
    !>  and changing the allocation method accordingly in the implementation to,
    !>  \code{.F90}
    !>      allocate(character(len(array, kind = IK)) :: count(lenLoc))
    !>  \endcode
    !>  However, this introduces `internal compiler error: Segmentation fault` with gfortran versions 10 and 11.<br>
    !>  Here is a code snippet to regenerate the bug in Intel ifort (uncomment the commented line to reproduce the gfortran bug),<br>
    !>  \code{.F90}
    !>
    !>      module pm_explicitLenResult
    !>          implicit none
    !>          interface
    !>              pure module function bug(array) result(count)
    !>                  character(*, SK), intent(in), contiguous    :: array(:)
    !>                  character(len(array),SK)    , allocatable   :: count(:) ! catastrophic internal error with ifort 2021.2. Fine with gfortran 10.3
    !>                 !character(:, SK)            , allocatable   :: count(:) ! catastrophic internal error with gfortran 10.3. Fine with ifort 2021.2
    !>              end function
    !>          end interface
    !>      end module pm_explicitLenResult
    !>
    !>      submodule (pm_explicitLenResult) routines
    !>          implicit none
    !>      contains
    !>          module procedure bug
    !>             allocate(count, source = array)
    !>          end procedure
    !>      end submodule routines
    !>
    !>      program main
    !>          use pm_explicitLenResult, only: bug
    !>          character(2) :: array(3) = ["AA", "BB", "CC"]
    !>          character(2), allocatable :: count(:)
    !>          count = bug(array)
    !>      end program main
    !>
    !>  \endcode
    !>  It turns out that both gfortran and Intel do not tolerate the separation of interface from implementation in the above code snippet.<br>
    !>  \remedy
    !>  If one duplicates the interface in the implementation submodule, then both compilers compile and run the code with no errors.<br>
    !>  This is the remedy that is currently used in this [getCountLoc](@ref pm_arrayFind::getCountLoc) generic interface (interface duplication where the bug exists).<br>
    !>  Here is a workaround example for the bug in the above code snippet,<br>
    !>  \code{.F90}
    !>
    !>      module pm_explicitLenResult
    !>          implicit none
    !>          interface
    !>              pure module function bug(array) result(count)
    !>                  character(*, SK), intent(in), contiguous    :: array(:)
    !>                  character(len(array),SK), allocatable   :: count(:) ! catastrophic internal error with ifort 2021.2. Fine with gfortran 10.3
    !>              end function
    !>          end interface
    !>      end module pm_explicitLenResult
    !>
    !>      submodule (pm_explicitLenResult) routines
    !>          implicit none
    !>      contains
    !>          module procedure bug
    !>             allocate(count, source = array)
    !>          end procedure
    !>      end submodule routines
    !>
    !>      program main
    !>          use pm_explicitLenResult, only: bug
    !>          character(2) :: array(3) = ["AA", "BB", "CC"]
    !>          character(2), allocatable :: count(:)
    !>          count = bug(array)
    !>      end program main
    !>
    !>  \endcode
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to higher-dimensional input arrays.<br>
    !>
    !>  \todo
    !>  \pvhigh
    !>  Currently, the value of `blindness` is checked for being non-zero in the implementation.<br>
    !>  However, the documentation of `blindness` requires it to be positive.<br>
    !>  This conflict between the implementation and documentation must be resolved.<br>
    !>
    !>  \todo
    !>  \pmed
    !>  The functionality of this generic interface can be extended with an optional `border` argument as in [getCountLoc](@ref pm_arrayFind::getCountLoc).<br>
    !>
    !>  \final{getCountLoc}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! border default D0_D0

    interface getCountLoc

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getCountLocDefBorDefCom_D0_D0_SK5(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK4_ENABLED
    PURE module function getCountLocDefBorDefCom_D0_D0_SK4(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK3_ENABLED
    PURE module function getCountLocDefBorDefCom_D0_D0_SK3(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK2_ENABLED
    PURE module function getCountLocDefBorDefCom_D0_D0_SK2(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK1_ENABLED
    PURE module function getCountLocDefBorDefCom_D0_D0_SK1(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getCountLocDefBorCusCom_D0_D0_SK5(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK4_ENABLED
    module function getCountLocDefBorCusCom_D0_D0_SK4(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK3_ENABLED
    module function getCountLocDefBorCusCom_D0_D0_SK3(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK2_ENABLED
    module function getCountLocDefBorCusCom_D0_D0_SK2(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK1_ENABLED
    module function getCountLocDefBorCusCom_D0_D0_SK1(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! border default D1_D0

    interface getCountLoc

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_SK5(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK4_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_SK4(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK3_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_SK3(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK2_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_SK2(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK1_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_SK1(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_IK5(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK4_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_IK4(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK3_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_IK3(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK2_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_IK2(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK1_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_IK1(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_LK5(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK4_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_LK4(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK3_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_LK3(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK2_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_LK2(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK1_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_LK1(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_CK5(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK4_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_CK4(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK3_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_CK3(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK2_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_CK2(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK1_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_CK1(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_RK5(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK4_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_RK4(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK3_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_RK3(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK2_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_RK2(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK1_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D0_RK1(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_SK5(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK4_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_SK4(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK3_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_SK3(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK2_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_SK2(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK1_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_SK1(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_IK5(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK4_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_IK4(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK3_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_IK3(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK2_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_IK2(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK1_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_IK1(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_LK5(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK4_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_LK4(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK3_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_LK3(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK2_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_LK2(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK1_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_LK1(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_CK5(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK4_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_CK4(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK3_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_CK3(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK2_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_CK2(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK1_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_CK1(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_RK5(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK4_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_RK4(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK3_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_RK3(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK2_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_RK2(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK1_ENABLED
    module function getCountLocDefBorCusCom_D1_D0_RK1(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! border default D1_D1

    interface getCountLoc

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_SK5(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK4_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_SK4(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK3_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_SK3(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK2_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_SK2(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK1_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_SK1(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_IK5(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK4_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_IK4(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK3_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_IK3(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK2_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_IK2(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK1_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_IK1(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_LK5(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK4_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_LK4(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK3_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_LK3(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK2_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_LK2(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK1_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_LK1(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_CK5(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK4_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_CK4(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK3_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_CK3(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK2_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_CK2(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK1_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_CK1(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_RK5(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK4_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_RK4(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK3_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_RK3(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK2_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_RK2(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK1_ENABLED
    PURE module function getCountLocDefBorDefCom_D1_D1_RK1(array, pattern, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorDefCom_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_SK5(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK4_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_SK4(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK3_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_SK3(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK2_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_SK2(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK1_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_SK1(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_IK5(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK4_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_IK4(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK3_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_IK3(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK2_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_IK2(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK1_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_IK1(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_LK5(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK4_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_LK4(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK3_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_LK3(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK2_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_LK2(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK1_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_LK1(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_CK5(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK4_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_CK4(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK3_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_CK3(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK2_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_CK2(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK1_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_CK1(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_RK5(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK4_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_RK4(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK3_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_RK3(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK2_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_RK2(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK1_ENABLED
    module function getCountLocDefBorCusCom_D1_D1_RK1(array, pattern, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDefBorCusCom_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! border discrete D0_D0

    interface getCountLoc

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getCountLocDisBorDefCom_D0_D0_SK5(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK4_ENABLED
    PURE module function getCountLocDisBorDefCom_D0_D0_SK4(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK3_ENABLED
    PURE module function getCountLocDisBorDefCom_D0_D0_SK3(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK2_ENABLED
    PURE module function getCountLocDisBorDefCom_D0_D0_SK2(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK1_ENABLED
    PURE module function getCountLocDisBorDefCom_D0_D0_SK1(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getCountLocDisBorCusCom_D0_D0_SK5(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK4_ENABLED
    module function getCountLocDisBorCusCom_D0_D0_SK4(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK3_ENABLED
    module function getCountLocDisBorCusCom_D0_D0_SK3(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK2_ENABLED
    module function getCountLocDisBorCusCom_D0_D0_SK2(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK1_ENABLED
    module function getCountLocDisBorCusCom_D0_D0_SK1(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! border discrete D1_D0

    interface getCountLoc

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_SK5(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK4_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_SK4(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK3_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_SK3(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK2_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_SK2(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK1_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_SK1(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_IK5(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK4_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_IK4(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK3_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_IK3(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK2_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_IK2(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK1_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_IK1(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_LK5(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK4_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_LK4(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK3_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_LK3(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK2_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_LK2(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK1_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_LK1(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_CK5(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK4_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_CK4(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK3_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_CK3(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK2_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_CK2(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK1_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_CK1(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_RK5(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK4_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_RK4(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK3_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_RK3(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK2_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_RK2(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK1_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D0_RK1(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_SK5(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK4_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_SK4(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK3_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_SK3(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK2_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_SK2(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK1_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_SK1(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_IK5(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK4_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_IK4(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK3_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_IK3(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK2_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_IK2(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK1_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_IK1(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_LK5(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK4_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_LK4(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK3_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_LK3(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK2_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_LK2(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK1_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_LK1(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_CK5(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK4_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_CK4(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK3_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_CK3(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK2_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_CK2(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK1_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_CK1(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_RK5(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK4_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_RK4(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK3_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_RK3(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK2_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_RK2(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK1_ENABLED
    module function getCountLocDisBorCusCom_D1_D0_RK1(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! border discrete D1_D1

    interface getCountLoc

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_SK5(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK4_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_SK4(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK3_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_SK3(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK2_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_SK2(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK1_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_SK1(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_IK5(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK4_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_IK4(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK3_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_IK3(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK2_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_IK2(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK1_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_IK1(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_LK5(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK4_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_LK4(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK3_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_LK3(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK2_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_LK2(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK1_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_LK1(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_CK5(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK4_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_CK4(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK3_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_CK3(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK2_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_CK2(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK1_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_CK1(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_RK5(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK4_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_RK4(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK3_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_RK3(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK2_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_RK2(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK1_ENABLED
    PURE module function getCountLocDisBorDefCom_D1_D1_RK1(array, pattern, border, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorDefCom_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_SK5(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK4_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_SK4(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK3_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_SK3(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK2_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_SK2(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if SK1_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_SK1(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_IK5(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK4_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_IK4(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK3_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_IK3(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK2_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_IK2(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if IK1_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_IK1(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_LK5(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK4_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_LK4(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK3_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_LK3(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK2_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_LK2(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if LK1_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_LK1(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_CK5(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK4_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_CK4(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK3_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_CK3(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK2_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_CK2(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if CK1_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_CK1(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_RK5(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK4_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_RK4(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK3_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_RK3(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK2_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_RK2(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

#if RK1_ENABLED
    module function getCountLocDisBorCusCom_D1_D1_RK1(array, pattern, border, iseq, blindness) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountLocDisBorCusCom_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        type(discrete_type)     , intent(in)                :: border
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                                         :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return an `allocatable` array containing the indices of the locations within the input array where the input `pattern` exists.
    !>
    !>  \details
    !>  If an input vector of `instance` is specified, containing the indices of the specific instances of `pattern` that must be found,
    !>  then the indices of only those specific instances will be returned.<br>
    !>  The instances of `pattern` are found via linear search.<br>
    !>  Therefore, the procedures under this generic interface have a worst-case complexity of `O(size(array))`.
    !>
    !>  \param[in]  array       :   The input `contiguous` array of rank `1` of either<br>
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL, or <br>
    !>                                  <li>    type `integer` of kind \IKALL, or <br>
    !>                                  <li>    type `logical` of kind \LKALL, or <br>
    !>                                  <li>    type `complex` of kind \CKALL, or <br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ol>
    !>                              or
    !>                              <ol>
    !>                                  <li>    a scalar assumed-length `character` of kind \SKALL,<br>
    !>                              </ol>
    !>                              within which the starting indices of the requested instances of `pattern` is to be found.<br>
    !>  \param[in]  pattern     :   The input object of the same or lower rank than the input `array`, and of the same type
    !>                              and kind as `array` containing the pattern that must be found within the input `array`.<br>
    !>  \param      iseq        :   The `external` user-specified function that takes two input **explicit-shape** arguments of the same type
    !>                              and kind as the input `array` and possibly, also the length of the arguments as the third argument, if the arguments are array-valued.<br>
    !>                              It returns a scalar `logical` of default kind \LK that is `.true.` if all elements of the two input arguments are equivalent (e.g., equal)
    !>                              according to the user-defined criterion, otherwise, it is `.false.`.<br>
    !>                              If `pattern` is an array of rank `1`, then the last argument to `iseq` is the length of the input `pattern`, preceded by a
    !>                              segment of `array` and `pattern` as the first and second arguments, whose lengths are given by the third argument `lenPattern`.<br>
    !>                              The following illustrates the generic interface of `iseq` where `pattern` is array-valued,
    !>                              \code{.F90}
    !>                                  function iseq(Segment, pattern, lenPattern) result(equivalent)
    !>                                      use pm_kind, only: IK, LK
    !>                                      integer(IK) , intent(in)    :: lenPattern
    !>                                      TYPE(KIND)  , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                      logical(LK)                 :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,<br>
    !>                              \code{.F90}
    !>                                  character(*, SK), intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  integer(IK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  logical(LK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  complex(CK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  real(RK)        , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              The following illustrates the generic interface of `iseq` where `pattern` is scalar-valued (**including Fortran scalar strings**),<br>
    !>                              \code{.F90}
    !>                                  function iseq(segment, pattern) result(equivalent)
    !>                                      use pm_kind, only: LK
    !>                                      TYPE(KIND)  , intent(in)    :: segment, pattern
    !>                                      logical(LK)                 :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,<br>
    !>                              \code{.F90}
    !>                                  use pm_kind, only: SK, IK, LK, CK, RK
    !>                                  character(*, SK), intent(in)    :: segment, pattern
    !>                                  integer(IK)     , intent(in)    :: segment, pattern
    !>                                  logical(LK)     , intent(in)    :: segment, pattern
    !>                                  complex(CK)     , intent(in)    :: segment, pattern
    !>                                  real(RK)        , intent(in)    :: segment, pattern
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              This user-defined equivalence check is extremely useful where a user-defined equivalence test other than exact equality
    !>                              or identity is needed, for example, when the array segments should match the input `pattern` only within a given threshold or,
    !>                              when the case-sensitivity in character comparisons do not matter.<br>
    !>                              In such cases, user can define a custom equivalence criterion within the user-defined external function `iseq` to achieve the goal.<br>
    !>                              (**optional**, the default equivalence operator is `.eqv.` if the input `array` is `logical`, otherwise `==`.)
    !>  \param[in]  instance    :   The input `contiguous` array of rank `1` of type `integer` of default kind \IK,
    !>                              containing the specific instances of the input `pattern` in the input `array` that should be found.<br>
    !>                              Any element of `instance` that points to an out-of-scope instance of `pattern` in the input `array` will be ignored.<br>
    !>                              Any element of `instance` that is negatively valued will be counted from end of the input `array`.<br>
    !>                              Any element of `instance` that is duplicated will yield duplicated output `array` indices if `pattern` is found, otherwise will be ignored.<br>
    !>                              For example, `instance = [2,-1]` requests finding the indices of second instance of `pattern` in `array` from the beginning and
    !>                              the first instance of `pattern` starting from the end of `array`.<br>
    !>                              (**optional**, the default value corresponds to finding all instances of `pattern` in `array`.)
    !>  \param[in]  sorted      :   The input `logical` of default kind \LK indicating whether the elements of the specified input `instance` are all in **ascending-order**.<br>
    !>                              This includes the negative elements of `instance` **after** they are translated to the corresponding **positive** instances from the beginning of the input `array`.<br>
    !>                              Setting `sorted = .true.` will lead to faster runtime of the procedure. However, the onus will be
    !>                              strictly on the user to ensure all elements of `instance` are in **ascending-order**.<br>
    !>                              This is generally not an easy guarantee to make if there are negative elements in `instance`.<br>
    !>                              Therefore, set `sorted = .true.` **only if** you can guarantee the validity of the condition.<br>
    !>                              See also the relevant benchmark at [pm_arrayFind](@ref pm_arrayFind).<br>
    !>                              (**optional**, default = `.false.`. It can be present as input argument **only if** the input argument `instance` is present.)
    !>  \param[in]  positive    :   The input `logical` of default kind \LK indicating whether the elements of the specified input `instance` are all positive.<br>
    !>                              Setting `positive = .true.` will lead to faster runtime of the procedure.<br>
    !>                              However, the onus will be strictly on the user to ensure all elements of `instance` are positive.<br>
    !>                              This is generally not an easy guarantee to make if there are negative elements in `instance`.<br>
    !>                              Therefore, set `positive = .true.` **only if** you can guarantee the validity of the condition.<br>
    !>                              See also the relevant benchmark at [pm_arrayFind](@ref pm_arrayFind).<br>
    !>                              (**optional**, default = `.false.`. It can be present as input argument **only if** the input argument `instance` is present.)
    !>  \param[in]  blindness   :   The input **positive** `integer` of default kind \IK representing the length of the segment of `array` that should be ignored after finding an instance of the input `pattern` in the array.<br>
    !>                              Setting `blindness = len(pattern)` (for assumed-length character `pattern`) or `blindness = size(pattern)` (for other types of array-valued `pattern`)
    !>                              will lead to a search for exclusive non-overlapping instances of `pattern` in the input `array`.<br>
    !>                              See the examples below for more illustration of the utility of this input argument.<br>
    !>                              (**optional**, default = `1_IK`)
    !>
    !>  \return
    !>  `loc`                   :   The output `allocatable` array of shape `(1:*)` of type `integer` of default kind \IK,
    !>                              containing the indices of the input `array` where the requested instances of `pattern` start.<br>
    !>                              The output `loc` is a vector of size `0` if there are no instances of the input `pattern` in the `array`.<br>
    !>
    !>  \interface{getLoc}
    !>  \code{.F90}
    !>
    !>      use pm_arrayFind, only: getLoc
    !>
    !>      ! `array` and `pattern` are scalar strings.
    !>
    !>      loc = getLoc(array, pattern, blindness = blindness)
    !>      loc = getLoc(array, pattern, iseq, blindness = blindness)
    !>      loc = getLoc(array, pattern, instance(:), sorted = sorted, positive = positive, blindness = blindness)
    !>      loc = getLoc(array, pattern, iseq, instance(:), sorted = sorted, positive = positive, blindness = blindness)
    !>
    !>      ! `pattern` is a scalar of the same type and kind as `array`.
    !>
    !>      loc = getLoc(array(:), pattern, blindness = blindness)
    !>      loc = getLoc(array(:), pattern, iseq, blindness = blindness)
    !>      loc = getLoc(array(:), pattern, instance(:), sorted = sorted, positive = positive, blindness = blindness)
    !>      loc = getLoc(array(:), pattern, iseq, instance(:), sorted = sorted, positive = positive, blindness = blindness)
    !>
    !>      ! `array` and `pattern` are both vectors.
    !>
    !>      loc = getLoc(array(:), pattern(:), blindness = blindness)
    !>      loc = getLoc(array(:), pattern(:), iseq, blindness = blindness)
    !>      loc = getLoc(array(:), pattern(:), instance(:), sorted = sorted, positive = positive, blindness = blindness)
    !>      loc = getLoc(array(:), pattern(:), iseq, instance(:), sorted = sorted, positive = positive, blindness = blindness)
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < blindness` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warning
    !>  The procedures under this generic interface are `impure` when the user-specified `external` procedure `iseq` is specified as input argument.<br>
    !>
    !>  \warning
    !>  Note that in Fortran, trailing blanks are ignored in character comparison, that is, `"Fortran" == "Fortran "` yields `.true.`.<br>
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  The functions under this generic interface might be slightly slower than [setLoc](@ref pm_arrayFind::setLoc)
    !>  subroutine implementations due to the extra copy action required upon return from the function.<br>
    !>  See [pm_arrayFind](@ref pm_arrayFind) for the relevant benchmarks.<br>
    !>
    !>  \note
    !>  If the input `array` does not contain any (of the requested) instances of the input `pattern`,
    !>  the output `loc` will be an `allocated` array of size `0`.<br>
    !>
    !>  \see
    !>  [setLoc](@ref pm_arrayFind::setLoc)<br>
    !>  [getBin](@ref pm_arraySearch::getBin)<br>
    !>  [setReplaced](@ref pm_arrayReplace::setReplaced)<br>
    !>  [getReplaced](@ref pm_arrayReplace::getReplaced)<br>
    !>  [setInserted](@ref pm_arrayInsert::setInserted)<br>
    !>  [setSplit](@ref pm_arraySplit::setSplit)<br>
    !>
    !>  \example{getLoc}
    !>  \include{lineno} example/pm_arrayFind/getLoc/main.F90
    !>  \compilef{getLoc}
    !>  \output{getLoc}
    !>  \include{lineno} example/pm_arrayFind/getLoc/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayFind](@ref test_pm_arrayFind)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.2.0}, \gfortran{10-12}
    !>  \desc
    !>  The Intel Fortran compiler Classic 2021.2.0 has a bug for the following interface definition
    !>  \code{.F90}
    !>      character(len(array),SK), allocatable   :: loc(:)
    !>  \endcode
    !>  leading to an **internal compiler error**.<br>
    !>  For now, the remedy seems to be to redefine the interface as,<br>
    !>  \code{.F90}
    !>      character(:, SK), allocatable   :: loc(:)
    !>  \endcode
    !>  and changing the allocation method accordingly in the implementation to,
    !>  \code{.F90}
    !>      allocate(character(len(array, kind = IK)) :: loc(lenLoc))
    !>  \endcode
    !>  However, this introduces `internal compiler error: Segmentation fault` with gfortran versions 10 and 11.<br>
    !>  Here is a code snippet to regenerate the bug in Intel ifort (uncomment the commented line to reproduce the gfortran bug),<br>
    !>  \code{.F90}
    !>
    !>      module pm_explicitLenResult
    !>          implicit none
    !>          interface
    !>              pure module function bug(array) result(loc)
    !>                  character(*, SK), intent(in), contiguous    :: array(:)
    !>                  character(len(array),SK)    , allocatable   :: loc(:) ! catastrophic internal error with ifort 2021.2. Fine with gfortran 10.3
    !>                 !character(:, SK)            , allocatable   :: loc(:) ! catastrophic internal error with gfortran 10.3. Fine with ifort 2021.2
    !>              end function
    !>          end interface
    !>      end module pm_explicitLenResult
    !>
    !>      submodule (pm_explicitLenResult) routines
    !>          implicit none
    !>      contains
    !>          module procedure bug
    !>             allocate(loc, source = array)
    !>          end procedure
    !>      end submodule routines
    !>
    !>      program main
    !>          use pm_explicitLenResult, only: bug
    !>          character(2) :: array(3) = ["AA", "BB", "CC"]
    !>          character(2), allocatable :: loc(:)
    !>          loc = bug(array)
    !>      end program main
    !>
    !>  \endcode
    !>  It turns out that both gfortran and Intel do not tolerate the separation of interface from implementation in the above code snippet.<br>
    !>  \remedy
    !>  If one duplicates the interface in the implementation submodule, then both compilers compile and run the code with no errors.<br>
    !>  This is the remedy that is currently used in this [getLoc](@ref pm_arrayFind::getLoc) generic interface (interface duplication where the bug exists).<br>
    !>  Here is a workaround example for the bug in the above code snippet,<br>
    !>  \code{.F90}
    !>
    !>      module pm_explicitLenResult
    !>          implicit none
    !>          interface
    !>              pure module function bug(array) result(loc)
    !>                  character(*, SK), intent(in), contiguous    :: array(:)
    !>                  character(len(array),SK), allocatable   :: loc(:) ! catastrophic internal error with ifort 2021.2. Fine with gfortran 10.3
    !>              end function
    !>          end interface
    !>      end module pm_explicitLenResult
    !>
    !>      submodule (pm_explicitLenResult) routines
    !>          implicit none
    !>      contains
    !>          module procedure bug
    !>             allocate(loc, source = array)
    !>          end procedure
    !>      end submodule routines
    !>
    !>      program main
    !>          use pm_explicitLenResult, only: bug
    !>          character(2) :: array(3) = ["AA", "BB", "CC"]
    !>          character(2), allocatable :: loc(:)
    !>          loc = bug(array)
    !>      end program main
    !>
    !>  \endcode
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to higher-dimensional input arrays.<br>
    !>
    !>  \todo
    !>  \pvhigh
    !>  Currently, the value of `blindness` is checked for being non-zero in the implementation.<br>
    !>  However, the documentation of `blindness` requires it to be positive.<br>
    !>  This conflict between the implementation and documentation must be resolved.<br>
    !>
    !>  \todo
    !>  \pmed
    !>  The functionality of this generic interface can be extended with an optional `border` argument as in [getCountLoc](@ref pm_arrayFind::getCountLoc).<br>
    !>
    !>  \final{getLoc}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! D0_D0

    interface getLoc

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getLocDefComDefIns_D0_D0_SK5(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getLocDefComDefIns_D0_D0_SK4(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getLocDefComDefIns_D0_D0_SK3(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getLocDefComDefIns_D0_D0_SK2(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK1_ENABLED
    PURE module function getLocDefComDefIns_D0_D0_SK1(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getLocCusComDefIns_D0_D0_SK5(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK4_ENABLED
    module function getLocCusComDefIns_D0_D0_SK4(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK3_ENABLED
    module function getLocCusComDefIns_D0_D0_SK3(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK2_ENABLED
    module function getLocCusComDefIns_D0_D0_SK2(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK1_ENABLED
    module function getLocCusComDefIns_D0_D0_SK1(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getLocDefComCusIns_D0_D0_SK5(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getLocDefComCusIns_D0_D0_SK4(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getLocDefComCusIns_D0_D0_SK3(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getLocDefComCusIns_D0_D0_SK2(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK1_ENABLED
    PURE module function getLocDefComCusIns_D0_D0_SK1(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getLocCusComCusIns_D0_D0_SK5(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK4_ENABLED
    module function getLocCusComCusIns_D0_D0_SK4(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK3_ENABLED
    module function getLocCusComCusIns_D0_D0_SK3(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK2_ENABLED
    module function getLocCusComCusIns_D0_D0_SK2(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK1_ENABLED
    module function getLocCusComCusIns_D0_D0_SK1(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D1_D0

    interface getLoc

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_SK5(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_SK4(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_SK3(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_SK2(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK1_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_SK1(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_IK5(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_IK4(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_IK3(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_IK2(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK1_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_IK1(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_LK5(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK4_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_LK4(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK3_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_LK3(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK2_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_LK2(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK1_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_LK1(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_CK5(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_CK4(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_CK3(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_CK2(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK1_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_CK1(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_RK5(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_RK4(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_RK3(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_RK2(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK1_ENABLED
    PURE module function getLocDefComDefIns_D1_D0_RK1(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getLocCusComDefIns_D1_D0_SK5(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK4_ENABLED
    module function getLocCusComDefIns_D1_D0_SK4(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK3_ENABLED
    module function getLocCusComDefIns_D1_D0_SK3(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK2_ENABLED
    module function getLocCusComDefIns_D1_D0_SK2(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK1_ENABLED
    module function getLocCusComDefIns_D1_D0_SK1(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getLocCusComDefIns_D1_D0_IK5(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK4_ENABLED
    module function getLocCusComDefIns_D1_D0_IK4(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK3_ENABLED
    module function getLocCusComDefIns_D1_D0_IK3(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK2_ENABLED
    module function getLocCusComDefIns_D1_D0_IK2(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK1_ENABLED
    module function getLocCusComDefIns_D1_D0_IK1(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getLocCusComDefIns_D1_D0_LK5(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK4_ENABLED
    module function getLocCusComDefIns_D1_D0_LK4(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK3_ENABLED
    module function getLocCusComDefIns_D1_D0_LK3(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK2_ENABLED
    module function getLocCusComDefIns_D1_D0_LK2(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK1_ENABLED
    module function getLocCusComDefIns_D1_D0_LK1(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getLocCusComDefIns_D1_D0_CK5(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK4_ENABLED
    module function getLocCusComDefIns_D1_D0_CK4(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK3_ENABLED
    module function getLocCusComDefIns_D1_D0_CK3(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK2_ENABLED
    module function getLocCusComDefIns_D1_D0_CK2(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK1_ENABLED
    module function getLocCusComDefIns_D1_D0_CK1(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getLocCusComDefIns_D1_D0_RK5(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK4_ENABLED
    module function getLocCusComDefIns_D1_D0_RK4(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK3_ENABLED
    module function getLocCusComDefIns_D1_D0_RK3(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK2_ENABLED
    module function getLocCusComDefIns_D1_D0_RK2(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK1_ENABLED
    module function getLocCusComDefIns_D1_D0_RK1(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_SK5(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_SK4(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_SK3(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_SK2(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK1_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_SK1(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_IK5(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_IK4(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_IK3(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_IK2(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK1_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_IK1(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_LK5(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK4_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_LK4(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK3_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_LK3(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK2_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_LK2(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK1_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_LK1(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_CK5(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_CK4(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_CK3(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_CK2(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK1_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_CK1(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_RK5(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_RK4(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_RK3(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_RK2(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK1_ENABLED
    PURE module function getLocDefComCusIns_D1_D0_RK1(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getLocCusComCusIns_D1_D0_SK5(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK4_ENABLED
    module function getLocCusComCusIns_D1_D0_SK4(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK3_ENABLED
    module function getLocCusComCusIns_D1_D0_SK3(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK2_ENABLED
    module function getLocCusComCusIns_D1_D0_SK2(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK1_ENABLED
    module function getLocCusComCusIns_D1_D0_SK1(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getLocCusComCusIns_D1_D0_IK5(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK4_ENABLED
    module function getLocCusComCusIns_D1_D0_IK4(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK3_ENABLED
    module function getLocCusComCusIns_D1_D0_IK3(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK2_ENABLED
    module function getLocCusComCusIns_D1_D0_IK2(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK1_ENABLED
    module function getLocCusComCusIns_D1_D0_IK1(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getLocCusComCusIns_D1_D0_LK5(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK4_ENABLED
    module function getLocCusComCusIns_D1_D0_LK4(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK3_ENABLED
    module function getLocCusComCusIns_D1_D0_LK3(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK2_ENABLED
    module function getLocCusComCusIns_D1_D0_LK2(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK1_ENABLED
    module function getLocCusComCusIns_D1_D0_LK1(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getLocCusComCusIns_D1_D0_CK5(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK4_ENABLED
    module function getLocCusComCusIns_D1_D0_CK4(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK3_ENABLED
    module function getLocCusComCusIns_D1_D0_CK3(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK2_ENABLED
    module function getLocCusComCusIns_D1_D0_CK2(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK1_ENABLED
    module function getLocCusComCusIns_D1_D0_CK1(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getLocCusComCusIns_D1_D0_RK5(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK4_ENABLED
    module function getLocCusComCusIns_D1_D0_RK4(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK3_ENABLED
    module function getLocCusComCusIns_D1_D0_RK3(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK2_ENABLED
    module function getLocCusComCusIns_D1_D0_RK2(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK1_ENABLED
    module function getLocCusComCusIns_D1_D0_RK1(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: pattern
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D1_D1

    interface getLoc

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_SK5(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_SK4(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_SK3(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_SK2(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK1_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_SK1(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_IK5(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_IK4(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_IK3(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_IK2(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK1_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_IK1(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_LK5(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK4_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_LK4(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK3_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_LK3(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK2_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_LK2(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK1_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_LK1(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_CK5(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_CK4(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_CK3(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_CK2(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK1_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_CK1(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_RK5(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_RK4(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_RK3(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_RK2(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK1_ENABLED
    PURE module function getLocDefComDefIns_D1_D1_RK1(array, pattern, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComDefIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getLocCusComDefIns_D1_D1_SK5(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK4_ENABLED
    module function getLocCusComDefIns_D1_D1_SK4(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK3_ENABLED
    module function getLocCusComDefIns_D1_D1_SK3(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK2_ENABLED
    module function getLocCusComDefIns_D1_D1_SK2(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK1_ENABLED
    module function getLocCusComDefIns_D1_D1_SK1(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getLocCusComDefIns_D1_D1_IK5(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK4_ENABLED
    module function getLocCusComDefIns_D1_D1_IK4(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK3_ENABLED
    module function getLocCusComDefIns_D1_D1_IK3(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK2_ENABLED
    module function getLocCusComDefIns_D1_D1_IK2(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK1_ENABLED
    module function getLocCusComDefIns_D1_D1_IK1(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getLocCusComDefIns_D1_D1_LK5(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK4_ENABLED
    module function getLocCusComDefIns_D1_D1_LK4(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK3_ENABLED
    module function getLocCusComDefIns_D1_D1_LK3(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK2_ENABLED
    module function getLocCusComDefIns_D1_D1_LK2(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK1_ENABLED
    module function getLocCusComDefIns_D1_D1_LK1(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getLocCusComDefIns_D1_D1_CK5(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK4_ENABLED
    module function getLocCusComDefIns_D1_D1_CK4(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK3_ENABLED
    module function getLocCusComDefIns_D1_D1_CK3(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK2_ENABLED
    module function getLocCusComDefIns_D1_D1_CK2(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK1_ENABLED
    module function getLocCusComDefIns_D1_D1_CK1(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getLocCusComDefIns_D1_D1_RK5(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK4_ENABLED
    module function getLocCusComDefIns_D1_D1_RK4(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK3_ENABLED
    module function getLocCusComDefIns_D1_D1_RK3(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK2_ENABLED
    module function getLocCusComDefIns_D1_D1_RK2(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK1_ENABLED
    module function getLocCusComDefIns_D1_D1_RK1(array, pattern, iseq, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComDefIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_SK5(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_SK4(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_SK3(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_SK2(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK1_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_SK1(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_IK5(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_IK4(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_IK3(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_IK2(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK1_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_IK1(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_LK5(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK4_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_LK4(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK3_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_LK3(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK2_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_LK2(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK1_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_LK1(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_CK5(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_CK4(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_CK3(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_CK2(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK1_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_CK1(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_RK5(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_RK4(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_RK3(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_RK2(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK1_ENABLED
    PURE module function getLocDefComCusIns_D1_D1_RK1(array, pattern, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocDefComCusIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getLocCusComCusIns_D1_D1_SK5(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK4_ENABLED
    module function getLocCusComCusIns_D1_D1_SK4(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK3_ENABLED
    module function getLocCusComCusIns_D1_D1_SK3(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK2_ENABLED
    module function getLocCusComCusIns_D1_D1_SK2(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if SK1_ENABLED
    module function getLocCusComCusIns_D1_D1_SK1(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getLocCusComCusIns_D1_D1_IK5(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK4_ENABLED
    module function getLocCusComCusIns_D1_D1_IK4(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK3_ENABLED
    module function getLocCusComCusIns_D1_D1_IK3(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK2_ENABLED
    module function getLocCusComCusIns_D1_D1_IK2(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if IK1_ENABLED
    module function getLocCusComCusIns_D1_D1_IK1(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getLocCusComCusIns_D1_D1_LK5(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK4_ENABLED
    module function getLocCusComCusIns_D1_D1_LK4(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK3_ENABLED
    module function getLocCusComCusIns_D1_D1_LK3(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK2_ENABLED
    module function getLocCusComCusIns_D1_D1_LK2(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if LK1_ENABLED
    module function getLocCusComCusIns_D1_D1_LK1(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in), contiguous    :: array(:)
        logical(LKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getLocCusComCusIns_D1_D1_CK5(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK4_ENABLED
    module function getLocCusComCusIns_D1_D1_CK4(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK3_ENABLED
    module function getLocCusComCusIns_D1_D1_CK3(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK2_ENABLED
    module function getLocCusComCusIns_D1_D1_CK2(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if CK1_ENABLED
    module function getLocCusComCusIns_D1_D1_CK1(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getLocCusComCusIns_D1_D1_RK5(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK4_ENABLED
    module function getLocCusComCusIns_D1_D1_RK4(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK3_ENABLED
    module function getLocCusComCusIns_D1_D1_RK3(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK2_ENABLED
    module function getLocCusComCusIns_D1_D1_RK2(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

#if RK1_ENABLED
    module function getLocCusComCusIns_D1_D1_RK1(array, pattern, iseq, instance, sorted, positive, blindness) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocCusComCusIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in), contiguous    :: pattern(:)
        procedure(logical(LK))                              :: iseq
        integer(IK)             , intent(in), contiguous    :: instance(:)
        logical(LK)             , intent(in), optional      :: sorted
        logical(LK)             , intent(in), optional      :: positive
        integer(IK)             , intent(in), optional      :: blindness
        integer(IK)                         , allocatable   :: loc(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return an `allocatable` array containing the indices of the locations within the input array where the input `pattern` exists.<br>
    !>
    !>  \details
    !>  If an input vector of `instance` is specified, containing the indices of the specific instances of `pattern` that must be found,
    !>  then, the indices of only those specific instances will be returned. The instances of `pattern` are found via linear search.<br>
    !>  Therefore, the procedures under this generic interface have a worst-case complexity of `O(size(array))`.<br>
    !>
    !>  \param[inout]   loc         :   The input/output `allocatable` array of shape `(:)` of type `integer` of default kind \IK
    !>                                  containing the indices of the input `array` where the requested instances of `pattern` start.
    !>                                  On input, `loc` must be allocated to the best-guess number of `pattern` instances that are to be found in the input `array`.<br>
    !>                                  If necessary, `loc` will be [resized](@ref pm_arrayResize::setResized) to fit all identified instances in the vector `loc`.<br>
    !>                                  On output, only the first `nloc` elements of `loc` are filled with the starting indices of `pattern` instances in `array`.<br>
    !>                                  This approach avoids a potentially unnecessary allocation at the beginning of search and a final redundant re-allocation,
    !>                                  thus improving the performance of the algorithm for repeated calls to this generic interface while ensuring correctness
    !>                                  at the cost of a small performance penalty for runtime resizing checks at iteration of the loop.<br>
    !>  \param[out]     nloc        :   The output scalar `integer` of default kind \IK, containing the number of pattern instances identified in `array`
    !>                                  whose starting position indices within `array` are stored and returned in `loc`.<br>
    !>                                  If there no instances of `pattern` in `array`, the return value for `nloc` is `0`.<br>
    !>  \param[in]      array       :   The input `contiguous` array of rank `1` of either<br>
    !>                                  <ol>
    !>                                      <li>    type `character` of kind \SKALL, or <br>
    !>                                      <li>    type `integer` of kind \IKALL, or <br>
    !>                                      <li>    type `logical` of kind \LKALL, or <br>
    !>                                      <li>    type `complex` of kind \CKALL, or <br>
    !>                                      <li>    type `real` of kind \RKALL, <br>
    !>                                  </ol>
    !>                                  or
    !>                                  <ol>
    !>                                      <li>    a scalar assumed-length `character` of kind \SKALL,<br>
    !>                                  </ol>
    !>                                  within which the starting indices of the requested instances of `pattern` is to be found.<br>
    !>  \param[in]      pattern     :   The input `contiguous` array of rank `1` or scalar of the same type and kind as the input `array`
    !>                                  containing the pattern that must be found within the input `array`.<br>
    !>  \param          iseq        :   The `external` user-specified function that takes two input **explicit-shape** arguments of the same type
    !>                                  and kind as the input `array` and possibly, also the length of the arguments as the third argument, if the arguments are array-valued.<br>
    !>                                  It returns a scalar `logical` of default kind \LK that is `.true.` if all elements of the two input arguments are equivalent (e.g., equal)
    !>                                  according to the user-defined criterion, otherwise, it is `.false.`.<br>
    !>                                  If `pattern` is an array of rank `1`, then the last argument
    !>                                  to `iseq` is the length of the input `pattern`, preceded by a segment of `array` and `pattern` as the first and
    !>                                  second arguments, whose lengths are given by the third argument `lenPattern`.<br>
    !>                                  The following illustrates the generic interface of `iseq` where `pattern` is array-valued,<br>
    !>                                  \code{.F90}
    !>                                      function iseq(Segment, pattern, lenPattern) result(equivalent)
    !>                                          use pm_kind, only: IK, LK
    !>                                          integer(IK) , intent(in)    :: lenPattern
    !>                                          TYPE(KIND)  , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                          logical(LK)                 :: equivalent
    !>                                      end function
    !>                                  \endcode
    !>                                  where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
    !>                                  \code{.F90}
    !>                                      use pm_kind, only: SK, IK, LK, CK, RK
    !>                                      character(*, SK), intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                      integer(IK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                      logical(LK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                      complex(CK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                      real(RK)        , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  \endcode
    !>                                  where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                                  The following illustrates the generic interface of `iseq` where `pattern` is scalar-valued (**including Fortran scalar strings**),
    !>                                  \code{.F90}
    !>                                      function iseq(segment, pattern) result(equivalent)
    !>                                          use pm_kind, only: LK
    !>                                          TYPE(KIND)  , intent(in)    :: segment, pattern
    !>                                          logical(LK)                 :: equivalent
    !>                                      end function
    !>                                  \endcode
    !>                                  where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
    !>                                  \code{.F90}
    !>                                      use pm_kind, only: SK, IK, LK, CK, RK
    !>                                      character(*, SK), intent(in)    :: segment, pattern
    !>                                      integer(IK)     , intent(in)    :: segment, pattern
    !>                                      logical(LK)     , intent(in)    :: segment, pattern
    !>                                      complex(CK)     , intent(in)    :: segment, pattern
    !>                                      real(RK)        , intent(in)    :: segment, pattern
    !>                                  \endcode
    !>                                  where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                                  This user-defined equivalence check is extremely useful where a user-defined equivalence test other than exact equality
    !>                                  or identity is needed, for example, when the array segments should match the input `pattern` only within a given threshold or,
    !>                                  when the case-sensitivity in character comparisons do not matter.<br>
    !>                                  In such cases, user can define a custom equivalence criterion within the user-defined external function `iseq` to achieve the goal.<br>
    !>                                  (**optional**, the default equivalence operator is `.eqv.` if the input `array` is `logical`, otherwise `==`.)
    !>  \param[in]      instance    :   The input `contiguous` array of rank `1` of non-zero values of type `integer` of default kind \IK,
    !>                                  containing the specific instances of the input `pattern` in the input `array` that should be found.<br>
    !>                                  Any element of `instance` that points to an out-of-scope instance of `pattern` in the input `array` will be ignored.<br>
    !>                                  Any element of `instance` that is negatively valued will be counted from end of the input `array`.<br>
    !>                                  Any element of `instance` that is duplicated will yield duplicated output `array` indices if `pattern` is found, otherwise will be ignored.<br>
    !>                                  For example, `instance = [2,-1]` requests finding the indices of second instance of `pattern` in `array` from the beginning and
    !>                                  the first instance of `pattern` starting from the end of `array`.<br>
    !>                                  (**optional**, the default value corresponds to finding all instances of `pattern` in `array`. It must be present **if and only if** the input arguments `sorted` and `positive` are also present.)
    !>  \param[in]      sorted      :   The input `logical` of default kind \LK indicating whether the elements of the specified input `instance` are all in **ascending-order**.<br>
    !>                                  This includes the negative elements of `instance` **after** they are translated to the corresponding **positive** instances from the beginning of the input `array`.<br>
    !>                                  Setting `sorted = .true.` will lead to faster runtime of the procedure.<br>
    !>                                  However, the onus will be strictly on the user to ensure all elements of `instance` are in **ascending-order**.<br>
    !>                                  This is generally not an easy guarantee to make if there are negative elements in `instance`.<br>
    !>                                  Therefore, set `sorted = .true.` **only if** you can guarantee the validity of the condition.<br>
    !>                                  See also the relevant benchmark at [pm_arrayFind](@ref pm_arrayFind).<br>
    !>                                  (**optional**. It must be present **if and only if** the input arguments `instance` and `positive` are also present.)
    !>  \param[in]      positive    :   The input `logical` of default kind \LK indicating whether the elements of the specified input `instance` are all positive.<br>
    !>                                  Setting `positive = .true.` will lead to faster runtime of the procedure.<br>
    !>                                  However, the onus will be strictly on the user to ensure all elements of `instance` are positive.<br>
    !>                                  This is generally not an easy guarantee to make if there are negative elements in `instance`.<br>
    !>                                  Therefore, set `positive = .true.` **only if** you can guarantee the validity of the condition.<br>
    !>                                  See also the relevant benchmark at [pm_arrayFind](@ref pm_arrayFind).<br>
    !>                                  (**optional**. It must be present **if and only if** the input arguments `sorted` and `positive` are also present.)
    !>  \param[in]      blindness   :   The input **positive** `integer` of default kind \IK representing the length of the segment of `array` that should be ignored after finding an instance of the input `pattern` in the array.<br>
    !>                                  Setting `blindness = len(pattern)` (for assumed-length character `pattern`) or `blindness = size(pattern)` (for other types of array-valued `pattern`)
    !>                                  will lead to a search for exclusive non-overlapping instances of `pattern` in the input `array`.<br>
    !>                                  See the examples below for more illustration of the utility of this input argument.<br>
    !>
    !>  \interface{setLoc}
    !>  \code{.F90}
    !>
    !>      use pm_arrayFind, only: setLoc
    !>
    !>      ! `array` and `pattern` are scalar strings.
    !>
    !>      call setLoc(loc(:), loc, array, pattern, blindness)
    !>      call setLoc(loc(:), loc, array, pattern, iseq, blindness)
    !>      call setLoc(loc(:), loc, array, pattern, instance(:), sorted, positive, blindness)
    !>      call setLoc(loc(:), loc, array, pattern, iseq, instance(:), sorted, positive, blindness)
    !>
    !>      ! `pattern` is a scalar.
    !>
    !>      call setLoc(loc(:), loc, array(:), pattern, blindness)
    !>      call setLoc(loc(:), loc, array(:), pattern, iseq, blindness)
    !>      call setLoc(loc(:), loc, array(:), pattern, instance(:), sorted, positive, blindness)
    !>      call setLoc(loc(:), loc, array(:), pattern, iseq, instance(:), sorted, positive, blindness)
    !>
    !>      ! `array` and `pattern` are both vectors.
    !>
    !>      call setLoc(loc(:), loc, array(:), pattern(:), blindness)
    !>      call setLoc(loc(:), loc, array(:), pattern(:), iseq, blindness)
    !>      call setLoc(loc(:), loc, array(:), pattern(:), instance(:), sorted, positive, blindness)
    !>      call setLoc(loc(:), loc, array(:), pattern(:), iseq, instance(:), sorted, positive, blindness)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < blindness` must hold for the corresponding input arguments.<br>
    !>  The condition `allocated(loc)` must hold for the corresponding input arguments.<br>
    !>  The condition `positive .and. all(0 < instance)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warning
    !>  The procedures under this generic interface are `impure` when the user-specified `external` procedure `iseq` is specified as input argument.<br>
    !>
    !>  \warning
    !>  Note that in Fortran, trailing blanks are ignored in character comparison, that is, `"Fortran" == "Fortran "` yields `.true.`.<br>
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  The subroutines under this generic interface is slightly faster than [getLoc](@ref pm_arrayFind::getLoc)
    !>  function implementations due to avoiding multiple allocations and copy actions.<br>
    !>  See [pm_arrayFind](@ref pm_arrayFind) for the relevant benchmarks.<br>
    !>
    !>  \note
    !>  This algorithm has been carefully designed to take into account the possibility of a lower bound for `array` other than `1`.<br>
    !>
    !>  \see
    !>  [getLoc](@ref pm_arrayFind::getLoc)<br>
    !>  [getBin](@ref pm_arraySearch::getBin)<br>
    !>  [setReplaced](@ref pm_arrayReplace::setReplaced)<br>
    !>  [getReplaced](@ref pm_arrayReplace::getReplaced)<br>
    !>  [setInserted](@ref pm_arrayInsert::setInserted)<br>
    !>  [setSplit](@ref pm_arraySplit::setSplit)<br>
    !>
    !>  \example{setLoc}
    !>  \include{lineno} example/pm_arrayFind/setLoc/main.F90
    !>  \compilef{setLoc}
    !>  \output{setLoc}
    !>  \include{lineno} example/pm_arrayFind/setLoc/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayFind](@ref test_pm_arrayFind)
    !>
    !>  \todo
    !>  \plow
    !>  This generic interface can be extended to higher-dimensional input arrays.<br>
    !>
    !>  \todo
    !>  \pvhigh
    !>  Currently, the value of `blindness` is checked for being non-zero in the implementation.<br>
    !>  However, the documentation of `blindness` requires it to be positive.<br>
    !>  This conflict between the implementation and documentation must be resolved.<br>
    !>
    !>  \final{setLoc}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! D0_D0

    interface setLoc

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setLocDefComDefIns_D0_D0_SK5(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setLocDefComDefIns_D0_D0_SK4(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setLocDefComDefIns_D0_D0_SK3(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setLocDefComDefIns_D0_D0_SK2(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setLocDefComDefIns_D0_D0_SK1(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setLocCusComDefIns_D0_D0_SK5(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setLocCusComDefIns_D0_D0_SK4(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setLocCusComDefIns_D0_D0_SK3(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setLocCusComDefIns_D0_D0_SK2(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setLocCusComDefIns_D0_D0_SK1(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setLocDefComCusIns_D0_D0_SK5(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setLocDefComCusIns_D0_D0_SK4(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setLocDefComCusIns_D0_D0_SK3(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setLocDefComCusIns_D0_D0_SK2(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setLocDefComCusIns_D0_D0_SK1(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setLocCusComCusIns_D0_D0_SK5(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setLocCusComCusIns_D0_D0_SK4(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setLocCusComCusIns_D0_D0_SK3(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setLocCusComCusIns_D0_D0_SK2(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setLocCusComCusIns_D0_D0_SK1(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D1_D0

    interface setLoc

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_SK5(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_SK4(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_SK3(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_SK2(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_SK1(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_IK5(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_IK4(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_IK3(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_IK2(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_IK1(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_LK5(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_LK4(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_LK3(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_LK2(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_LK1(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_CK5(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_CK4(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_CK3(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_CK2(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_CK1(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_RK5(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_RK4(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_RK3(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_RK2(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D0_RK1(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_SK5(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_SK4(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_SK3(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_SK2(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_SK1(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_IK5(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_IK4(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_IK3(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_IK2(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_IK1(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_LK5(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_LK4(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_LK3(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_LK2(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_LK1(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_CK5(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_CK4(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_CK3(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_CK2(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_CK1(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_RK5(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_RK4(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_RK3(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_RK2(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setLocCusComDefIns_D1_D0_RK1(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_SK5(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_SK4(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_SK3(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_SK2(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_SK1(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_IK5(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_IK4(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_IK3(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_IK2(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_IK1(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_LK5(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_LK4(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_LK3(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_LK2(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_LK1(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_CK5(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_CK4(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_CK3(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_CK2(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_CK1(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_RK5(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_RK4(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_RK3(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_RK2(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D0_RK1(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_SK5(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_SK4(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_SK3(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_SK2(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_SK1(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_IK5(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_IK4(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_IK3(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_IK2(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_IK1(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_LK5(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_LK4(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_LK3(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_LK2(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_LK1(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_CK5(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_CK4(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_CK3(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_CK2(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_CK1(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_RK5(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_RK4(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_RK3(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_RK2(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setLocCusComCusIns_D1_D0_RK1(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D1_D1

    interface setLoc

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_SK5(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_SK4(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_SK3(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_SK2(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_SK1(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_IK5(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_IK4(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_IK3(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_IK2(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_IK1(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_LK5(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_LK4(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_LK3(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_LK2(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_LK1(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_CK5(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_CK4(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_CK3(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_CK2(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_CK1(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_RK5(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_RK4(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_RK3(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_RK2(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setLocDefComDefIns_D1_D1_RK1(loc, nloc, array, pattern, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComDefIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_SK5(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_SK4(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_SK3(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_SK2(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_SK1(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_IK5(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_IK4(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_IK3(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_IK2(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_IK1(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_LK5(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_LK4(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_LK3(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_LK2(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_LK1(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_CK5(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_CK4(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_CK3(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_CK2(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_CK1(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_RK5(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_RK4(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_RK3(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_RK2(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setLocCusComDefIns_D1_D1_RK1(loc, nloc, array, pattern, iseq, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComDefIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , optional      :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_SK5(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_SK4(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_SK3(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_SK2(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_SK1(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_IK5(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_IK4(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_IK3(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_IK2(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_IK1(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_LK5(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_LK4(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_LK3(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_LK2(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_LK1(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_CK5(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_CK4(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_CK3(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_CK2(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_CK1(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_RK5(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_RK4(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_RK3(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_RK2(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setLocDefComCusIns_D1_D1_RK1(loc, nloc, array, pattern, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocDefComCusIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_SK5(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_SK4(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_SK3(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_SK2(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_SK1(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_IK5(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_IK4(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_IK3(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_IK2(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_IK1(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_LK5(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_LK4(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_LK3(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_LK2(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_LK1(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_CK5(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_CK4(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_CK3(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_CK2(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_CK1(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_RK5(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_RK4(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_RK3(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_RK2(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setLocCusComCusIns_D1_D1_RK1(loc, nloc, array, pattern, iseq, instance, sorted, positive, blindness)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLocCusComCusIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)                    :: sorted
        logical(LK)             , intent(in)                    :: positive
        integer(IK)             , intent(in)                    :: blindness
        integer(IK)             , intent(inout) , allocatable   :: loc(:)
        integer(IK)             , intent(out)                   :: nloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayFind ! LCOV_EXCL_LINE