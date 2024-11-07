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
!>  This module contains procedures and generic interfaces for removing a pattern from arrays of various types at the specified instances of occurrence of pattern.
!>
!>  \benchmarks
!>
!>  \benchmark{setRemoved-scalarPattern_vs_vectorPattern, The runtime performance of [setRemoved](@ref pm_arrayRemove::setRemoved) for scalar vs. vector input `pattern` argument.}
!>  \include{lineno} benchmark/pm_arrayRemove/setRemoved-scalarPattern_vs_vectorPattern/main.F90
!>  \compilefb{setRemoved-scalarPattern_vs_vectorPattern}
!>  \postprocb{setRemoved-scalarPattern_vs_vectorPattern}
!>  \include{lineno} benchmark/pm_arrayRemove/setRemoved-scalarPattern_vs_vectorPattern/main.py
!>  \visb{setRemoved-scalarPattern_vs_vectorPattern}
!>  \image html benchmark/pm_arrayRemove/setRemoved-scalarPattern_vs_vectorPattern/benchmark.setRemoved-scalarPattern_vs_vectorPattern.runtime.png width=1000
!>  \image html benchmark/pm_arrayRemove/setRemoved-scalarPattern_vs_vectorPattern/benchmark.setRemoved-scalarPattern_vs_vectorPattern.runtime.ratio.png width=1000
!>  \moralb{setRemoved-scalarPattern_vs_vectorPattern}
!>  -#  The procedures under the generic interface [setRemoved](@ref pm_arrayRemove::setRemoved) take both scalar and vector `pattern` arguments.<br>
!>      As evidenced by the above benchmark, when the input `pattern` is vector of length `1`, it is much faster, **up to 4X**,
!>      to pass `pattern` as a scalar instead of a whole array of length `1`.<br>
!>      Note that this benchmark is likely irrelevant to removing substrings from Fortran strings.<br>
!>
!>  \benchmark{getRemoved_vs_setRemoved, The runtime performance of [getRemoved](@ref pm_arrayRemove::getRemoved) vs. [setRemoved](@ref pm_arrayRemove::setRemoved)}
!>  \include{lineno} benchmark/pm_arrayRemove/getRemoved_vs_setRemoved/main.F90
!>  \compilefb{getRemoved_vs_setRemoved}
!>  \postprocb{getRemoved_vs_setRemoved}
!>  \include{lineno} benchmark/pm_arrayRemove/getRemoved_vs_setRemoved/main.py
!>  \visb{getRemoved_vs_setRemoved}
!>  \image html benchmark/pm_arrayRemove/getRemoved_vs_setRemoved/benchmark.getRemoved_vs_setRemoved.runtime.png width=1000
!>  \image html benchmark/pm_arrayRemove/getRemoved_vs_setRemoved/benchmark.getRemoved_vs_setRemoved.runtime.ratio.png width=1000
!>  \moralb{getRemoved_vs_setRemoved}
!>  -#  The procedures under the generic interface [getRemoved](@ref pm_arrayRemove::getRemoved) are functions while
!>      the procedures under the generic interface [setRemoved](@ref pm_arrayRemove::setRemoved) are subroutines.<br>
!>      From the benchmark results, it appears that the functional interface performs slightly less efficiently than the subroutine interface,
!>      despite the two algorithms having the same implementation.
!>  -#  Note that this benchmark does not even include the cost of repeated reallcations, that is, the allocation of `Removed` happen only once in all tests.
!>  -#  Note that this benchmark considers the worst-case scenario where all elements of the input `array` match the
!>      input `pattern` and must be therefore, removed.
!>
!>  \test
!>  [test_pm_arrayRemove](@ref test_pm_arrayRemove)
!>
!>  \bug
!>  \status \unresolved
!>  \source \ifort{2021.2.0}
!>  \desc
!>  The \ifort{2021.2.0} has a bug for `use`ing the following two modules
!>  simultaneously in the implementation of the procedures in this module,
!>  \code{.F90}
!>      use pm_arrayUnique, only: getUnique
!>      use pm_arraySort, only: setSorted
!>  \endcode
!>  The following is example error message from the \ifort,
!>  \verbatim
!>      pm_arrayRemove@routines@setRemoved_D1.inc.F90(196): error #6405: The same named entity from different modules and/or program units cannot be referenced.   [UNIQUE]
!>                      if (present(unique)) unique_def = unique
!>      ----------------------------^
!>  \endverbatim
!>  Searching this error on the web points to the possibility that the internal representation of entities by \ifort has a naming conflict.
!>  \remedy
!>  For now, the remedy was to isolate the use of one of the modules to exactly where it is needed, like,
!>  \code{.F90}
!>      block
!>          use pm_arrayUnique, only: getUnique
!>          InstanceNew = getUnique(InstanceNew(1:lenInstanceNew))
!>      end block
!>  \endcode
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayRemove

!>  \cond excluded
!   \bug
!   The following bypasses the bug reported below that creates a conflict between Intel and gfortran.
#if     __INTEL_COMPILER
#define LEN_ARRAY :
#else
#define LEN_ARRAY len(array,IK)
#endif
!>  \endcond excluded

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_arrayRemove"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return an `allocatable` array containing the remaining parts of the input array as a sequence after removing
    !>  the input `pattern` at the requested occurrences.
    !>
    !>  \details
    !>  If an input vector of `instance` is specified, containing the indices of the specific instances of `pattern` that must be removed,
    !>  then only those instances will be removed from the array.
    !>
    !>  \param[in]  array   :       The input `contiguous` array of rank `1` of either <br>
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
    !>                              within which the requested instances of `pattern` is to be removed.<br>
    !>  \param[in]  pattern :       The input `contiguous` array of rank `1` or scalar of the same type and kind as the input `array`
    !>                              containing the pattern that must be removed from the input `array`.
    !>  \param      iseq    :       The `external` user-specified function that takes two input **explicit-shape** arguments of the same type
    !>                              and kind as the input `array` and possibly, also the length of the arguments as the third argument.<br>
    !>                              It returns a scalar `logical` of default kind \LK that is `.true.` if all elements of the two input arguments are equivalent (e.g., equal) according to the
    !>                              user-defined criterion, otherwise, it is `.false.`.<br>
    !>                              If `pattern` is an array of rank `1`, then the last argument to `iseq`
    !>                              is the length of the input `pattern`, preceded by a segment of `array` and `pattern` as the first and second arguments
    !>                              whose lengths are given by the third argument `lenPattern`.<br>
    !>                              The following illustrates the generic interface of `iseq` where `pattern` is array-valued,
    !>                              \code{.F90}
    !>                                  function iseq(Segment, pattern, lenPattern) result(equivalent)
    !>                                      use pm_kind, only: IK, LK
    !>                                      integer(IK) , intent(in)    :: lenPattern
    !>                                      TYPE(KIND)  , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                      logical(LK)                 :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
    !>                              \code{.F90}
    !>                                  use pm_kind, only: SK, IK, LK, CK, RK
    !>                                  character(*, SK), intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  integer(IK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  logical(LK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  complex(CK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  real(RK)        , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              The following illustrates the generic interface of `iseq` where `pattern` is scalar-valued (**including Fortran scalar strings**),
    !>                              \code{.F90}
    !>                                  function iseq(segment, pattern) result(equivalent)
    !>                                      use pm_kind, only: LK
    !>                                      TYPE(KIND)  , intent(in)    :: segment, pattern
    !>                                      logical(LK)                 :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
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
    !>                              (**optional**, the default equivalence operator is `.eqv.` if the input `array` is `logical`, otherwise `==`)
    !>  \param[in]  instance    :   The input `contiguous` array of rank `1` of type `integer` of default kind \IK,
    !>                              containing the instances of the input `pattern` in the input `array` that should be removed.<br>
    !>                              Any element of `instance` that points to an out-of-scope instance of `pattern` in the input `array` will be ignored.<br>
    !>                              Any element of `instance` that is negatively valued will be counted from end of the input `array`.<br>
    !>                              For example, `instance = [2,-1]` requests removing the second instance of `pattern` in `array` from the beginning and
    !>                              removing the first instance of `pattern` starting from the end of `array`.<br>
    !>                              (**optional**, the default value corresponds to removing all instances of `pattern` in `array`)
    !>  \param[in]  sorted      :   The input `logical` of default kind \LK indicating whether the elements of the specified input `instance` are all in ascending-order.<br>
    !>                              This includes the negative elements of `instance` **after** they are translated to the corresponding **positive** instances from
    !>                              the beginning of the input `array`. Setting `sorted = .true.` will lead to faster runtime of the procedure.<br>
    !>                              However, the onus will be strictly on the user to ensure all elements of `instance` are in ascending-order.<br>
    !>                              This is generally not an easy guarantee to make if there are negative elements in `instance`.<br>
    !>                              Therefore, set `sorted = .true.` **only if** you can guarantee the validity of the condition.<br>
    !>                              (**optional**, default = `.false.`. It can be present as input argument **only if** the input argument `instance` is present.)
    !>  \param[in]  unique      :   The input `logical` of default kind \LK indicating whether the elements of the specified input `instance` are all unique.<br>
    !>                              This includes the negative elements of `instance` **after** they are translated to the corresponding **positive** instances from the beginning of the input `array`.<br>
    !>                              Setting `unique = .true.` will lead to faster runtime of the procedure.<br>
    !>                              However, the onus will be strictly on the user to ensure all elements of `instance` are unique.<br>
    !>                              This is generally not an easy guarantee to make if there are negative elements in `instance`.<br>
    !>                              Therefore, set `unique = .true.` **only if** you can guarantee the validity of the condition.<br>
    !>                              (**optional**, default = `.false.`. It can be present as input argument **only if** the input argument `instance` is present.)
    !>
    !>  \return
    !>  `ArrayRemoved`          :   The output allocatable of the same type, kind, and rank as the input `array` argument containing
    !>                              the remaining parts of the array in sequence after removing the requested (or all) instances of `pattern`.
    !>
    !>  \interface{getRemoved}
    !>  \code{.F90}
    !>
    !>      use pm_arrayRemove, only: getRemoved
    !>
    !>      ArrayRemoved = getRemoved(array, pattern)
    !>      ArrayRemoved = getRemoved(array, pattern, iseq)
    !>      ArrayRemoved = getRemoved(array, pattern, instance, sorted = sorted, unique = unique)
    !>      ArrayRemoved = getRemoved(array, pattern, iseq, instance, sorted = sorted, unique = unique)
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
    !>  The functions under this generic interface could be slightly slower than the [setRemoved](@ref pm_arrayRemove::setRemoved)
    !>  subroutine implementations due to the extra copy action required upon return from the function.<br>
    !>  See [pm_arrayRemove](@ref pm_arrayRemove) for the relevant benchmarks.<br>
    !>
    !>  \see
    !>  [setRemoved](@ref pm_arrayRemove::setRemoved)<br>
    !>  [setReplaced](@ref pm_arrayReplace::setReplaced)<br>
    !>  [getReplaced](@ref pm_arrayReplace::getReplaced)<br>
    !>  [setInserted](@ref pm_arrayInsert::setInserted)<br>
    !>  [setSplit](@ref pm_arraySplit::setSplit)<br>
    !>
    !>  \example{getRemoved}
    !>  \include{lineno} example/pm_arrayRemove/getRemoved/main.F90
    !>  \compilef{getRemoved}
    !>  \output{getRemoved}
    !>  \include{lineno} example/pm_arrayRemove/getRemoved/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayRemove](@ref test_pm_arrayRemove)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.2.0}, \gfortran{10, 11}
    !>  \desc
    !>  The \ifort{2021.2.0} has a bug for the following interface definition
    !>  \code{.F90}
    !>      character(len(array,IK),SK), allocatable   :: ArrayRemoved(:)
    !>  \endcode
    !>  leading to an **internal compiler error**. For now, the remedy seems to be to redefine the interface as,
    !>  \code{.F90}
    !>      character(:, SK), allocatable   :: ArrayRemoved(:)
    !>  \endcode
    !>  and changing the allocation method accordingly in the implementation to,
    !>  \code{.F90}
    !>      allocate(character(len(array, kind = IK)) :: ArrayRemoved(lenArrayRemoved))
    !>  \endcode
    !>  However, this introduces `internal compiler error: Segmentation fault` with \gfortran versions 10 and 11.
    !>  Here is a code snippet to regenerate the bug in \ifort (uncomment the commented line to reproduce the gfortran bug),
    !>  \code{.F90}
    !>
    !>      module pm_explicitLenResult
    !>          implicit none
    !>          interface
    !>              pure module function bug(array) result(ArrayRemoved)
    !>                  character(*, SK), intent(in), contiguous    :: array(:)
    !>                  character(len(array),SK)    , allocatable   :: ArrayRemoved(:) ! catastrophic internal error with ifort 2021.2. Fine with gfortran 10.3
    !>                 !character(:, SK)            , allocatable   :: ArrayRemoved(:) ! catastrophic internal error with gfortran 10.3. Fine with ifort 2021.2
    !>              end function
    !>          end interface
    !>      end module pm_explicitLenResult
    !>
    !>      submodule (pm_explicitLenResult) routines
    !>          implicit none
    !>      contains
    !>          module procedure bug
    !>             allocate(ArrayRemoved, source = array)
    !>          end procedure
    !>      end submodule routines
    !>
    !>      program main
    !>          use pm_explicitLenResult, only: bug
    !>          character(2) :: array(3) = ["AA", "BB", "CC"]
    !>          character(2), allocatable :: ArrayRemoved(:)
    !>          ArrayRemoved = bug(array)
    !>      end program main
    !>
    !>  \endcode
    !>  \remedy
    !>  It turns out that both gfortran and Intel do not tolerate the separation of interface from implementation in the above code snippet.<br>
    !>  If one duplicates the interface in the implementation submodule, then both compilers compile and run the code with no errors.<br>
    !>  This is the remedy that is currently used in this [getRemoved](@ref pm_arrayRemove::getRemoved) generic interface
    !>  (interface duplication where the bug exists).<br>
    !>  Here is a workaround example for the bug in the above code snippet,<br>
    !>  \code{.F90}
    !>
    !>      module pm_explicitLenResult
    !>          implicit none
    !>          interface
    !>              pure module function bug(array) result(ArrayRemoved)
    !>                  character(*, SK), intent(in), contiguous    :: array(:)
    !>                  character(len(array),SK), allocatable   :: ArrayRemoved(:) ! catastrophic internal error with ifort 2021.2. Fine with gfortran 10.3
    !>              end function
    !>          end interface
    !>      end module pm_explicitLenResult
    !>
    !>      submodule (pm_explicitLenResult) routines
    !>          implicit none
    !>      contains
    !>          module procedure bug
    !>             allocate(ArrayRemoved, source = array)
    !>          end procedure
    !>      end submodule routines
    !>
    !>      program main
    !>          use pm_explicitLenResult, only: bug
    !>          character(2) :: array(3) = ["AA", "BB", "CC"]
    !>          character(2), allocatable :: ArrayRemoved(:)
    !>          ArrayRemoved = bug(array)
    !>      end program main
    !>
    !>  \endcode
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to 2D input objects.<br>
    !>
    !>  \final{getRemoved}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getRemoved

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRemovedDefComDefIns_D0_D0_SK5(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(:,SKG)                        , allocatable   :: ArrayRemoved
    end function
#endif

#if SK4_ENABLED
    PURE module function getRemovedDefComDefIns_D0_D0_SK4(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(:,SKG)                        , allocatable   :: ArrayRemoved
    end function
#endif

#if SK3_ENABLED
    PURE module function getRemovedDefComDefIns_D0_D0_SK3(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(:,SKG)                        , allocatable   :: ArrayRemoved
    end function
#endif

#if SK2_ENABLED
    PURE module function getRemovedDefComDefIns_D0_D0_SK2(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(:,SKG)                        , allocatable   :: ArrayRemoved
    end function
#endif

#if SK1_ENABLED
    PURE module function getRemovedDefComDefIns_D0_D0_SK1(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(:,SKG)                        , allocatable   :: ArrayRemoved
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRemovedCusComDefIns_D0_D0_SK5(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        character(:,SKG)                        , allocatable   :: ArrayRemoved
    end function
#endif

#if SK4_ENABLED
    module function getRemovedCusComDefIns_D0_D0_SK4(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        character(:,SKG)                        , allocatable   :: ArrayRemoved
    end function
#endif

#if SK3_ENABLED
    module function getRemovedCusComDefIns_D0_D0_SK3(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        character(:,SKG)                        , allocatable   :: ArrayRemoved
    end function
#endif

#if SK2_ENABLED
    module function getRemovedCusComDefIns_D0_D0_SK2(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        character(:,SKG)                        , allocatable   :: ArrayRemoved
    end function
#endif

#if SK1_ENABLED
    module function getRemovedCusComDefIns_D0_D0_SK1(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        character(:,SKG)                        , allocatable   :: ArrayRemoved
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRemovedDefComCusIns_D0_D0_SK5(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(:,SKG)                        , allocatable   :: ArrayRemoved
    end function
#endif

#if SK4_ENABLED
    PURE module function getRemovedDefComCusIns_D0_D0_SK4(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(:,SKG)                        , allocatable   :: ArrayRemoved
    end function
#endif

#if SK3_ENABLED
    PURE module function getRemovedDefComCusIns_D0_D0_SK3(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(:,SKG)                        , allocatable   :: ArrayRemoved
    end function
#endif

#if SK2_ENABLED
    PURE module function getRemovedDefComCusIns_D0_D0_SK2(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(:,SKG)                        , allocatable   :: ArrayRemoved
    end function
#endif

#if SK1_ENABLED
    PURE module function getRemovedDefComCusIns_D0_D0_SK1(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(:,SKG)                        , allocatable   :: ArrayRemoved
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRemovedCusComCusIns_D0_D0_SK5(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(:,SKG)                        , allocatable   :: ArrayRemoved
    end function
#endif

#if SK4_ENABLED
    module function getRemovedCusComCusIns_D0_D0_SK4(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(:,SKG)                        , allocatable   :: ArrayRemoved
    end function
#endif

#if SK3_ENABLED
    module function getRemovedCusComCusIns_D0_D0_SK3(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(:,SKG)                        , allocatable   :: ArrayRemoved
    end function
#endif

#if SK2_ENABLED
    module function getRemovedCusComCusIns_D0_D0_SK2(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(:,SKG)                        , allocatable   :: ArrayRemoved
    end function
#endif

#if SK1_ENABLED
    module function getRemovedCusComCusIns_D0_D0_SK1(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(:,SKG)                        , allocatable   :: ArrayRemoved
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_SK5(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_SK4(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_SK3(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_SK2(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK1_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_SK1(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_IK5(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_IK4(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_IK3(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_IK2(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK1_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_IK1(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_LK5(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK4_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_LK4(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK3_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_LK3(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK2_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_LK2(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK1_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_LK1(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_CK5(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_CK4(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_CK3(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_CK2(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK1_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_CK1(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_RK5(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_RK4(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_RK3(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_RK2(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK1_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D0_RK1(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRemovedCusComDefIns_D1_D0_SK5(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK4_ENABLED
    module function getRemovedCusComDefIns_D1_D0_SK4(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK3_ENABLED
    module function getRemovedCusComDefIns_D1_D0_SK3(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK2_ENABLED
    module function getRemovedCusComDefIns_D1_D0_SK2(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK1_ENABLED
    module function getRemovedCusComDefIns_D1_D0_SK1(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getRemovedCusComDefIns_D1_D0_IK5(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK4_ENABLED
    module function getRemovedCusComDefIns_D1_D0_IK4(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK3_ENABLED
    module function getRemovedCusComDefIns_D1_D0_IK3(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK2_ENABLED
    module function getRemovedCusComDefIns_D1_D0_IK2(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK1_ENABLED
    module function getRemovedCusComDefIns_D1_D0_IK1(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getRemovedCusComDefIns_D1_D0_LK5(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK4_ENABLED
    module function getRemovedCusComDefIns_D1_D0_LK4(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK3_ENABLED
    module function getRemovedCusComDefIns_D1_D0_LK3(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK2_ENABLED
    module function getRemovedCusComDefIns_D1_D0_LK2(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK1_ENABLED
    module function getRemovedCusComDefIns_D1_D0_LK1(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getRemovedCusComDefIns_D1_D0_CK5(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK4_ENABLED
    module function getRemovedCusComDefIns_D1_D0_CK4(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK3_ENABLED
    module function getRemovedCusComDefIns_D1_D0_CK3(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK2_ENABLED
    module function getRemovedCusComDefIns_D1_D0_CK2(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK1_ENABLED
    module function getRemovedCusComDefIns_D1_D0_CK1(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getRemovedCusComDefIns_D1_D0_RK5(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK4_ENABLED
    module function getRemovedCusComDefIns_D1_D0_RK4(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK3_ENABLED
    module function getRemovedCusComDefIns_D1_D0_RK3(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK2_ENABLED
    module function getRemovedCusComDefIns_D1_D0_RK2(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK1_ENABLED
    module function getRemovedCusComDefIns_D1_D0_RK1(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_SK5(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_SK4(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_SK3(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_SK2(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK1_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_SK1(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_IK5(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_IK4(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_IK3(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_IK2(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK1_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_IK1(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_LK5(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK4_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_LK4(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK3_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_LK3(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK2_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_LK2(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK1_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_LK1(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_CK5(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_CK4(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_CK3(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_CK2(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK1_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_CK1(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_RK5(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_RK4(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_RK3(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_RK2(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK1_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D0_RK1(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRemovedCusComCusIns_D1_D0_SK5(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK4_ENABLED
    module function getRemovedCusComCusIns_D1_D0_SK4(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK3_ENABLED
    module function getRemovedCusComCusIns_D1_D0_SK3(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK2_ENABLED
    module function getRemovedCusComCusIns_D1_D0_SK2(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK1_ENABLED
    module function getRemovedCusComCusIns_D1_D0_SK1(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getRemovedCusComCusIns_D1_D0_IK5(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK4_ENABLED
    module function getRemovedCusComCusIns_D1_D0_IK4(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK3_ENABLED
    module function getRemovedCusComCusIns_D1_D0_IK3(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK2_ENABLED
    module function getRemovedCusComCusIns_D1_D0_IK2(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK1_ENABLED
    module function getRemovedCusComCusIns_D1_D0_IK1(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getRemovedCusComCusIns_D1_D0_LK5(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK4_ENABLED
    module function getRemovedCusComCusIns_D1_D0_LK4(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK3_ENABLED
    module function getRemovedCusComCusIns_D1_D0_LK3(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK2_ENABLED
    module function getRemovedCusComCusIns_D1_D0_LK2(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK1_ENABLED
    module function getRemovedCusComCusIns_D1_D0_LK1(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getRemovedCusComCusIns_D1_D0_CK5(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK4_ENABLED
    module function getRemovedCusComCusIns_D1_D0_CK4(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK3_ENABLED
    module function getRemovedCusComCusIns_D1_D0_CK3(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK2_ENABLED
    module function getRemovedCusComCusIns_D1_D0_CK2(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK1_ENABLED
    module function getRemovedCusComCusIns_D1_D0_CK1(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getRemovedCusComCusIns_D1_D0_RK5(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK4_ENABLED
    module function getRemovedCusComCusIns_D1_D0_RK4(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK3_ENABLED
    module function getRemovedCusComCusIns_D1_D0_RK3(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK2_ENABLED
    module function getRemovedCusComCusIns_D1_D0_RK2(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK1_ENABLED
    module function getRemovedCusComCusIns_D1_D0_RK1(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_SK5(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_SK4(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_SK3(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_SK2(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK1_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_SK1(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_IK5(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_IK4(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_IK3(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_IK2(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK1_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_IK1(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_LK5(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK4_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_LK4(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK3_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_LK3(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK2_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_LK2(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK1_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_LK1(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_CK5(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_CK4(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_CK3(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_CK2(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK1_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_CK1(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_RK5(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_RK4(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_RK3(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_RK2(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK1_ENABLED
    PURE module function getRemovedDefComDefIns_D1_D1_RK1(array, pattern) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComDefIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRemovedCusComDefIns_D1_D1_SK5(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK4_ENABLED
    module function getRemovedCusComDefIns_D1_D1_SK4(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK3_ENABLED
    module function getRemovedCusComDefIns_D1_D1_SK3(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK2_ENABLED
    module function getRemovedCusComDefIns_D1_D1_SK2(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK1_ENABLED
    module function getRemovedCusComDefIns_D1_D1_SK1(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getRemovedCusComDefIns_D1_D1_IK5(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK4_ENABLED
    module function getRemovedCusComDefIns_D1_D1_IK4(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK3_ENABLED
    module function getRemovedCusComDefIns_D1_D1_IK3(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK2_ENABLED
    module function getRemovedCusComDefIns_D1_D1_IK2(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK1_ENABLED
    module function getRemovedCusComDefIns_D1_D1_IK1(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getRemovedCusComDefIns_D1_D1_LK5(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK4_ENABLED
    module function getRemovedCusComDefIns_D1_D1_LK4(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK3_ENABLED
    module function getRemovedCusComDefIns_D1_D1_LK3(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK2_ENABLED
    module function getRemovedCusComDefIns_D1_D1_LK2(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK1_ENABLED
    module function getRemovedCusComDefIns_D1_D1_LK1(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getRemovedCusComDefIns_D1_D1_CK5(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK4_ENABLED
    module function getRemovedCusComDefIns_D1_D1_CK4(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK3_ENABLED
    module function getRemovedCusComDefIns_D1_D1_CK3(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK2_ENABLED
    module function getRemovedCusComDefIns_D1_D1_CK2(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK1_ENABLED
    module function getRemovedCusComDefIns_D1_D1_CK1(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getRemovedCusComDefIns_D1_D1_RK5(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK4_ENABLED
    module function getRemovedCusComDefIns_D1_D1_RK4(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK3_ENABLED
    module function getRemovedCusComDefIns_D1_D1_RK3(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK2_ENABLED
    module function getRemovedCusComDefIns_D1_D1_RK2(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK1_ENABLED
    module function getRemovedCusComDefIns_D1_D1_RK1(array, pattern, iseq) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComDefIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_SK5(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_SK4(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_SK3(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_SK2(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK1_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_SK1(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_IK5(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_IK4(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_IK3(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_IK2(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK1_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_IK1(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_LK5(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK4_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_LK4(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK3_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_LK3(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK2_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_LK2(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK1_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_LK1(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_CK5(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_CK4(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_CK3(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_CK2(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK1_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_CK1(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_RK5(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_RK4(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_RK3(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_RK2(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK1_ENABLED
    PURE module function getRemovedDefComCusIns_D1_D1_RK1(array, pattern, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedDefComCusIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRemovedCusComCusIns_D1_D1_SK5(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK4_ENABLED
    module function getRemovedCusComCusIns_D1_D1_SK4(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK3_ENABLED
    module function getRemovedCusComCusIns_D1_D1_SK3(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK2_ENABLED
    module function getRemovedCusComCusIns_D1_D1_SK2(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if SK1_ENABLED
    module function getRemovedCusComCusIns_D1_D1_SK1(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY, SKG)               , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getRemovedCusComCusIns_D1_D1_IK5(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK4_ENABLED
    module function getRemovedCusComCusIns_D1_D1_IK4(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK3_ENABLED
    module function getRemovedCusComCusIns_D1_D1_IK3(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK2_ENABLED
    module function getRemovedCusComCusIns_D1_D1_IK2(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if IK1_ENABLED
    module function getRemovedCusComCusIns_D1_D1_IK1(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getRemovedCusComCusIns_D1_D1_LK5(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK4_ENABLED
    module function getRemovedCusComCusIns_D1_D1_LK4(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK3_ENABLED
    module function getRemovedCusComCusIns_D1_D1_LK3(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK2_ENABLED
    module function getRemovedCusComCusIns_D1_D1_LK2(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if LK1_ENABLED
    module function getRemovedCusComCusIns_D1_D1_LK1(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getRemovedCusComCusIns_D1_D1_CK5(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK4_ENABLED
    module function getRemovedCusComCusIns_D1_D1_CK4(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK3_ENABLED
    module function getRemovedCusComCusIns_D1_D1_CK3(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK2_ENABLED
    module function getRemovedCusComCusIns_D1_D1_CK2(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if CK1_ENABLED
    module function getRemovedCusComCusIns_D1_D1_CK1(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getRemovedCusComCusIns_D1_D1_RK5(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK4_ENABLED
    module function getRemovedCusComCusIns_D1_D1_RK4(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK3_ENABLED
    module function getRemovedCusComCusIns_D1_D1_RK3(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK2_ENABLED
    module function getRemovedCusComCusIns_D1_D1_RK2(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

#if RK1_ENABLED
    module function getRemovedCusComCusIns_D1_D1_RK1(array, pattern, iseq, instance, sorted, unique) result(ArrayRemoved)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRemovedCusComCusIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: ArrayRemoved(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the remaining parts of the input array as a sequence after removing the input `pattern` at the requested occurrences.
    !>
    !>  \details
    !>  If an input vector of `instance` is specified, containing the indices of the specific instances of `pattern` that must be removed,
    !>  then only those instances will be removed from the array.
    !>
    !>  \param[inout]   array   :   The input/output `allocatable` array of rank `1` of either <br>
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
    !>                              within which the requested instances of `pattern` is to be removed.<br>
    !>                              On output, `array` is reallocated to contain the input `array` with the requested instances of `pattern` removed from it.
    !>  \param[in]      pattern :   The input `contiguous` array of rank `1` or scalar of the same type and kind as the input `array`
    !>                              containing the pattern that must be removed from the input `array`.
    !>  \param          iseq    :   The `external` user-specified function that takes two input **explicit-shape** arguments of the same type
    !>                              and kind as the input `array` and possibly, also the length of the arguments as the third argument.<br>
    !>                              It returns a scalar `logical` of default kind \LK that is `.true.` if all elements of the two input arguments are equivalent (e.g., equal) according to the
    !>                              user-defined criterion, otherwise, it is `.false.`.<br>
    !>                              If `pattern` is an array of rank `1`, then the last argument to `iseq`
    !>                              is the length of the input `pattern`, preceded by a segment of `array` and `pattern` as the first and second arguments
    !>                              whose lengths are given by the third argument `lenPattern`.<br>
    !>                              The following illustrates the generic interface of `iseq` where `pattern` is array-valued,
    !>                              \code{.F90}
    !>                                  function iseq(Segment, pattern, lenPattern) result(equivalent)
    !>                                      use pm_kind, only: IK, LK
    !>                                      integer(IK) , intent(in)    :: lenPattern
    !>                                      TYPE(KIND)  , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                      logical(LK)                 :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
    !>                              \code{.F90}
    !>                                  use pm_kind, only: SK, IK, LK, CK, RK
    !>                                  character(*, SK), intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  integer(IK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  logical(LK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  complex(CK)     , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                                  real(RK)        , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              The following illustrates the generic interface of `iseq` where `pattern` is scalar-valued (**including Fortran scalar strings**),
    !>                              \code{.F90}
    !>                                  function iseq(segment, pattern) result(equivalent)
    !>                                      use pm_kind, only: LK
    !>                                      TYPE(KIND)  , intent(in)    :: segment, pattern
    !>                                      logical(LK)                 :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
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
    !>                              (**optional**, the default equivalence operator is `.eqv.` if the input `array` is `logical`, otherwise `==`)
    !>  \param[in]  instance    :   The input `contiguous` array of rank `1` of type `integer` of default kind \IK,
    !>                              containing the instances of the input `pattern` in the input `array` that should be removed.<br>
    !>                              Any element of `instance` that points to an out-of-scope instance of `pattern` in the input `array` will be ignored.<br>
    !>                              Any element of `instance` that is negatively valued will be counted from end of the input `array`.<br>
    !>                              For example, `instance = [2,-1]` requests removing the second instance of `pattern` in `array` from the beginning and
    !>                              removing the first instance of `pattern` starting from the end of `array`.<br>
    !>                              (**optional**, the default value corresponds to removing all instances of `pattern` in `array`)
    !>  \param[in]  sorted      :   The input `logical` of default kind \LK indicating whether the elements of the specified input `instance` are all in ascending-order.<br>
    !>                              This includes the negative elements of `instance` **after** they are translated to the corresponding **positive** instances from
    !>                              the beginning of the input `array`.<br>Setting `sorted = .true.` will lead to faster runtime of the procedure.<br>
    !>                              However, the onus will be strictly on the user to ensure all elements of `instance` are in ascending-order.<br>
    !>                              This is generally not an easy guarantee to make if there are negative elements in `instance`.<br>
    !>                              Therefore, set `sorted = .true.` **only if** you can guarantee the validity of the condition.<br>
    !>                              (**optional**, default = `.false.`. It can be present as input argument **only if** the input argument `instance` is present.)
    !>  \param[in]  unique      :   The input `logical` of default kind \LK indicating whether the elements of the specified input `instance` are all unique.<br>
    !>                              This includes the negative elements of `instance` **after** they are translated to the corresponding **positive** instances from
    !>                              the beginning of the input `array`.<br>
    !>                              Setting `unique = .true.` will lead to faster runtime of the procedure.<br>
    !>                              However, the onus will be strictly on the user to ensure all elements of `instance` are unique.<br>
    !>                              This is generally not an easy guarantee to make if there are negative elements in `instance`.<br>
    !>                              Therefore, set `unique = .true.` **only if** you can guarantee the validity of the condition.<br>
    !>                              (**optional**, default = `.false.`. It can be present as input argument **only if** the input argument `instance` is present.)
    !>
    !>  \interface{setRemoved}
    !>  \code{.F90}
    !>
    !>      use pm_arrayRemove, only: setRemoved
    !>
    !>      call setRemoved(array, pattern)
    !>      call setRemoved(array, pattern, iseq)
    !>      call setRemoved(array, pattern, instance, sorted = sorted, unique = unique)
    !>      call setRemoved(array, pattern, iseq, instance, sorted = sorted, unique = unique)
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
    !>  The functions under this generic interface are slightly slower than the [setRemoved](@ref pm_arrayRemove::setRemoved) subroutine implementations.<br>
    !>  See [pm_arrayRemove](@ref pm_arrayRemove) for the relevant benchmarks.<br>
    !>
    !>  \note
    !>  Upon return, the output allocatable `array` is guaranteed to have the same lower bound as before, but its upper bound will likely be different.<br>
    !>
    !>  \see
    !>  [getRemoved](@ref pm_arrayRemove::getRemoved)<br>
    !>  [getReplaced](@ref pm_arrayReplace::getReplaced)<br>
    !>  [setReplaced](@ref pm_arrayReplace::setReplaced)<br>
    !>  [setInserted](@ref pm_arrayInsert::setInserted)<br>
    !>  [setSplit](@ref pm_arraySplit::setSplit)<br>
    !>
    !>  \example{setRemoved}
    !>  \include{lineno} example/pm_arrayRemove/setRemoved/main.F90
    !>  \compilef{setRemoved}
    !>  \output{setRemoved}
    !>  \include{lineno} example/pm_arrayRemove/setRemoved/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayRemove](@ref test_pm_arrayRemove)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to 2D input objects.<br>
    !>
    !>  \final{setRemoved}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setRemoved

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D0_D0_SK5(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D0_D0_SK4(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D0_D0_SK3(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D0_D0_SK2(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D0_D0_SK1(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setRemovedCusComDefIns_D0_D0_SK5(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setRemovedCusComDefIns_D0_D0_SK4(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setRemovedCusComDefIns_D0_D0_SK3(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setRemovedCusComDefIns_D0_D0_SK2(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setRemovedCusComDefIns_D0_D0_SK1(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D0_D0_SK5(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D0_D0_SK4(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D0_D0_SK3(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D0_D0_SK2(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D0_D0_SK1(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setRemovedCusComCusIns_D0_D0_SK5(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setRemovedCusComCusIns_D0_D0_SK4(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setRemovedCusComCusIns_D0_D0_SK3(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setRemovedCusComCusIns_D0_D0_SK2(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setRemovedCusComCusIns_D0_D0_SK1(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
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
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_SK5(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_SK4(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_SK3(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_SK2(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_SK1(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_IK5(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_IK4(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_IK3(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_IK2(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_IK1(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_LK5(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_LK4(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_LK3(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_LK2(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_LK1(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_CK5(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_CK4(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_CK3(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_CK2(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_CK1(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_RK5(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_RK4(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_RK3(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_RK2(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D0_RK1(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_SK5(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_SK4(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_SK3(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_SK2(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_SK1(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_IK5(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_IK4(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_IK3(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_IK2(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_IK1(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_LK5(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_LK4(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_LK3(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_LK2(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_LK1(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_CK5(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_CK4(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_CK3(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_CK2(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_CK1(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_RK5(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_RK4(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_RK3(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_RK2(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D0_RK1(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_SK5(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_SK4(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_SK3(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_SK2(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_SK1(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_IK5(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_IK4(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_IK3(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_IK2(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_IK1(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_LK5(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_LK4(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_LK3(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_LK2(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_LK1(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_CK5(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_CK4(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_CK3(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_CK2(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_CK1(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_RK5(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_RK4(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_RK3(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_RK2(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D0_RK1(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_SK5(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_SK4(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_SK3(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_SK2(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_SK1(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_IK5(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_IK4(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_IK3(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_IK2(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_IK1(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_LK5(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_LK4(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_LK3(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_LK2(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_LK1(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_CK5(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_CK4(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_CK3(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_CK2(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_CK1(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_RK5(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_RK4(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_RK3(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_RK2(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D0_RK1(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
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
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_SK5(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_SK4(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_SK3(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_SK2(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_SK1(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_IK5(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_IK4(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_IK3(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_IK2(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_IK1(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_LK5(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_LK4(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_LK3(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_LK2(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_LK1(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_CK5(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_CK4(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_CK3(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_CK2(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_CK1(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_RK5(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_RK4(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_RK3(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_RK2(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRemovedDefComDefIns_D1_D1_RK1(array, pattern)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_SK5(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_SK4(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_SK3(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_SK2(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_SK1(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_IK5(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_IK4(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_IK3(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_IK2(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_IK1(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_LK5(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_LK4(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_LK3(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_LK2(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_LK1(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_CK5(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_CK4(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_CK3(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_CK2(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_CK1(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_RK5(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_RK4(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_RK3(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_RK2(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setRemovedCusComDefIns_D1_D1_RK1(array, pattern, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_SK5(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_SK4(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_SK3(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_SK2(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_SK1(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_IK5(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_IK4(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_IK3(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_IK2(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_IK1(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_LK5(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_LK4(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_LK3(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_LK2(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_LK1(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_CK5(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_CK4(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_CK3(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_CK2(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_CK1(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_RK5(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_RK4(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_RK3(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_RK2(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRemovedDefComCusIns_D1_D1_RK1(array, pattern, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_SK5(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_SK4(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_SK3(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_SK2(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_SK1(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_IK5(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_IK4(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_IK3(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_IK2(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_IK1(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_LK5(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_LK4(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_LK3(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_LK2(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_LK1(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_CK5(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_CK4(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_CK3(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_CK2(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_CK1(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_RK5(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_RK4(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_RK3(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_RK2(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setRemovedCusComCusIns_D1_D1_RK1(array, pattern, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
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
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! The following interfaces return a new array and can be activated in the future if needed.

!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D0_D0_SK5(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D0_D0_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        character(:,SKG)        , intent(inout) , allocatable   :: array
!        character(*,SKG)        , intent(in)                    :: pattern
!        character(:,SKG)        , intent(out)   , allocatable   :: ArrayRemoved
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D0_D0_SK4(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D0_D0_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        character(:,SKG)        , intent(inout) , allocatable   :: array
!        character(*,SKG)        , intent(in)                    :: pattern
!        character(:,SKG)        , intent(out)   , allocatable   :: ArrayRemoved
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D0_D0_SK3(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D0_D0_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        character(:,SKG)        , intent(inout) , allocatable   :: array
!        character(*,SKG)        , intent(in)                    :: pattern
!        character(:,SKG)        , intent(out)   , allocatable   :: ArrayRemoved
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D0_D0_SK2(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D0_D0_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        character(:,SKG)        , intent(inout) , allocatable   :: array
!        character(*,SKG)        , intent(in)                    :: pattern
!        character(:,SKG)        , intent(out)   , allocatable   :: ArrayRemoved
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D0_D0_SK1(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D0_D0_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        character(:,SKG)        , intent(inout) , allocatable   :: array
!        character(*,SKG)        , intent(in)                    :: pattern
!        character(:,SKG)        , intent(out)   , allocatable   :: ArrayRemoved
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
!    module subroutine setRemovedCusComDefIns_D0_D0_SK5(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D0_D0_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        character(:,SKG)        , intent(inout) , allocatable   :: array
!        character(*,SKG)        , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        character(:,SKG)        , intent(out)   , allocatable   :: ArrayRemoved
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    module subroutine setRemovedCusComDefIns_D0_D0_SK4(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D0_D0_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        character(:,SKG)        , intent(inout) , allocatable   :: array
!        character(*,SKG)        , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        character(:,SKG)        , intent(out)   , allocatable   :: ArrayRemoved
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    module subroutine setRemovedCusComDefIns_D0_D0_SK3(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D0_D0_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        character(:,SKG)        , intent(inout) , allocatable   :: array
!        character(*,SKG)        , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        character(:,SKG)        , intent(out)   , allocatable   :: ArrayRemoved
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    module subroutine setRemovedCusComDefIns_D0_D0_SK2(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D0_D0_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        character(:,SKG)        , intent(inout) , allocatable   :: array
!        character(*,SKG)        , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        character(:,SKG)        , intent(out)   , allocatable   :: ArrayRemoved
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    module subroutine setRemovedCusComDefIns_D0_D0_SK1(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D0_D0_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        character(:,SKG)        , intent(inout) , allocatable   :: array
!        character(*,SKG)        , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        character(:,SKG)        , intent(out)   , allocatable   :: ArrayRemoved
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
!    PURE module subroutine setRemovedDefComCusIns_D0_D0_SK5(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D0_D0_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        character(:,SKG)        , intent(inout) , allocatable   :: array
!        character(*,SKG)        , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(:,SKG)        , intent(out)   , allocatable   :: ArrayRemoved
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D0_D0_SK4(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D0_D0_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        character(:,SKG)        , intent(inout) , allocatable   :: array
!        character(*,SKG)        , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(:,SKG)        , intent(out)   , allocatable   :: ArrayRemoved
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D0_D0_SK3(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D0_D0_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        character(:,SKG)        , intent(inout) , allocatable   :: array
!        character(*,SKG)        , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(:,SKG)        , intent(out)   , allocatable   :: ArrayRemoved
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D0_D0_SK2(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D0_D0_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        character(:,SKG)        , intent(inout) , allocatable   :: array
!        character(*,SKG)        , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(:,SKG)        , intent(out)   , allocatable   :: ArrayRemoved
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D0_D0_SK1(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D0_D0_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        character(:,SKG)        , intent(inout) , allocatable   :: array
!        character(*,SKG)        , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(:,SKG)        , intent(out)   , allocatable   :: ArrayRemoved
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
!    module subroutine setRemovedCusComCusIns_D0_D0_SK5(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D0_D0_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        character(:,SKG)        , intent(inout) , allocatable   :: array
!        character(*,SKG)        , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(:,SKG)        , intent(out)   , allocatable   :: ArrayRemoved
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    module subroutine setRemovedCusComCusIns_D0_D0_SK4(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D0_D0_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        character(:,SKG)        , intent(inout) , allocatable   :: array
!        character(*,SKG)        , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(:,SKG)        , intent(out)   , allocatable   :: ArrayRemoved
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    module subroutine setRemovedCusComCusIns_D0_D0_SK3(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D0_D0_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        character(:,SKG)        , intent(inout) , allocatable   :: array
!        character(*,SKG)        , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(:,SKG)        , intent(out)   , allocatable   :: ArrayRemoved
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    module subroutine setRemovedCusComCusIns_D0_D0_SK2(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D0_D0_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        character(:,SKG)        , intent(inout) , allocatable   :: array
!        character(*,SKG)        , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(:,SKG)        , intent(out)   , allocatable   :: ArrayRemoved
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    module subroutine setRemovedCusComCusIns_D0_D0_SK1(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D0_D0_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        character(:,SKG)        , intent(inout) , allocatable   :: array
!        character(*,SKG)        , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(:,SKG)        , intent(out)   , allocatable   :: ArrayRemoved
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_SK5(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)                    :: pattern
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_SK4(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)                    :: pattern
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_SK3(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)                    :: pattern
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_SK2(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)                    :: pattern
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_SK1(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)                    :: pattern
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_IK5(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_IK5
!#endif
!        use pm_kind, only: IKG => IK5
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)                    :: pattern
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_IK4(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_IK4
!#endif
!        use pm_kind, only: IKG => IK4
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)                    :: pattern
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_IK3(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_IK3
!#endif
!        use pm_kind, only: IKG => IK3
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)                    :: pattern
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_IK2(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_IK2
!#endif
!        use pm_kind, only: IKG => IK2
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)                    :: pattern
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_IK1(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_IK1
!#endif
!        use pm_kind, only: IKG => IK1
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)                    :: pattern
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if LK5_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_LK5(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_LK5
!#endif
!        use pm_kind, only: LKG => LK5
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)                    :: pattern
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK4_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_LK4(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_LK4
!#endif
!        use pm_kind, only: LKG => LK4
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)                    :: pattern
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK3_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_LK3(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_LK3
!#endif
!        use pm_kind, only: LKG => LK3
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)                    :: pattern
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK2_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_LK2(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_LK2
!#endif
!        use pm_kind, only: LKG => LK2
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)                    :: pattern
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK1_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_LK1(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_LK1
!#endif
!        use pm_kind, only: LKG => LK1
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)                    :: pattern
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if CK5_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_CK5(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_CK5
!#endif
!        use pm_kind, only: CKG => CK5
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)                    :: pattern
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK4_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_CK4(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_CK4
!#endif
!        use pm_kind, only: CKG => CK4
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)                    :: pattern
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK3_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_CK3(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_CK3
!#endif
!        use pm_kind, only: CKG => CK3
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)                    :: pattern
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK2_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_CK2(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_CK2
!#endif
!        use pm_kind, only: CKG => CK2
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)                    :: pattern
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK1_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_CK1(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_CK1
!#endif
!        use pm_kind, only: CKG => CK1
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)                    :: pattern
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_RK5(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)                    :: pattern
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_RK4(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)                    :: pattern
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_RK3(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)                    :: pattern
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_RK2(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)                    :: pattern
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D0_RK1(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D0_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)                    :: pattern
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
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
!    module subroutine setRemovedCusComDefIns_D1_D0_SK5(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_SK4(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_SK3(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_SK2(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_SK1(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_IK5(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_IK5
!#endif
!        use pm_kind, only: IKG => IK5
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_IK4(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_IK4
!#endif
!        use pm_kind, only: IKG => IK4
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_IK3(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_IK3
!#endif
!        use pm_kind, only: IKG => IK3
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_IK2(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_IK2
!#endif
!        use pm_kind, only: IKG => IK2
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_IK1(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_IK1
!#endif
!        use pm_kind, only: IKG => IK1
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if LK5_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_LK5(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_LK5
!#endif
!        use pm_kind, only: LKG => LK5
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK4_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_LK4(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_LK4
!#endif
!        use pm_kind, only: LKG => LK4
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK3_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_LK3(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_LK3
!#endif
!        use pm_kind, only: LKG => LK3
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK2_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_LK2(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_LK2
!#endif
!        use pm_kind, only: LKG => LK2
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK1_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_LK1(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_LK1
!#endif
!        use pm_kind, only: LKG => LK1
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if CK5_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_CK5(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_CK5
!#endif
!        use pm_kind, only: CKG => CK5
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK4_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_CK4(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_CK4
!#endif
!        use pm_kind, only: CKG => CK4
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK3_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_CK3(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_CK3
!#endif
!        use pm_kind, only: CKG => CK3
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK2_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_CK2(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_CK2
!#endif
!        use pm_kind, only: CKG => CK2
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK1_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_CK1(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_CK1
!#endif
!        use pm_kind, only: CKG => CK1
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_RK5(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_RK4(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_RK3(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_RK2(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D0_RK1(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D0_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
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
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_SK5(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_SK4(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_SK3(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_SK2(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_SK1(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_IK5(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_IK5
!#endif
!        use pm_kind, only: IKG => IK5
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_IK4(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_IK4
!#endif
!        use pm_kind, only: IKG => IK4
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_IK3(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_IK3
!#endif
!        use pm_kind, only: IKG => IK3
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_IK2(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_IK2
!#endif
!        use pm_kind, only: IKG => IK2
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_IK1(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_IK1
!#endif
!        use pm_kind, only: IKG => IK1
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if LK5_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_LK5(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_LK5
!#endif
!        use pm_kind, only: LKG => LK5
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK4_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_LK4(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_LK4
!#endif
!        use pm_kind, only: LKG => LK4
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK3_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_LK3(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_LK3
!#endif
!        use pm_kind, only: LKG => LK3
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK2_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_LK2(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_LK2
!#endif
!        use pm_kind, only: LKG => LK2
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK1_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_LK1(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_LK1
!#endif
!        use pm_kind, only: LKG => LK1
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if CK5_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_CK5(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_CK5
!#endif
!        use pm_kind, only: CKG => CK5
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK4_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_CK4(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_CK4
!#endif
!        use pm_kind, only: CKG => CK4
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK3_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_CK3(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_CK3
!#endif
!        use pm_kind, only: CKG => CK3
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK2_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_CK2(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_CK2
!#endif
!        use pm_kind, only: CKG => CK2
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK1_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_CK1(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_CK1
!#endif
!        use pm_kind, only: CKG => CK1
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_RK5(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_RK4(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_RK3(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_RK2(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D0_RK1(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D0_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)                    :: pattern
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
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
!    module subroutine setRemovedCusComCusIns_D1_D0_SK5(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_SK4(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_SK3(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_SK2(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_SK1(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_IK5(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_IK5
!#endif
!        use pm_kind, only: IKG => IK5
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_IK4(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_IK4
!#endif
!        use pm_kind, only: IKG => IK4
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_IK3(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_IK3
!#endif
!        use pm_kind, only: IKG => IK3
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_IK2(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_IK2
!#endif
!        use pm_kind, only: IKG => IK2
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_IK1(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_IK1
!#endif
!        use pm_kind, only: IKG => IK1
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if LK5_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_LK5(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_LK5
!#endif
!        use pm_kind, only: LKG => LK5
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK4_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_LK4(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_LK4
!#endif
!        use pm_kind, only: LKG => LK4
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK3_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_LK3(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_LK3
!#endif
!        use pm_kind, only: LKG => LK3
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK2_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_LK2(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_LK2
!#endif
!        use pm_kind, only: LKG => LK2
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK1_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_LK1(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_LK1
!#endif
!        use pm_kind, only: LKG => LK1
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if CK5_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_CK5(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_CK5
!#endif
!        use pm_kind, only: CKG => CK5
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK4_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_CK4(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_CK4
!#endif
!        use pm_kind, only: CKG => CK4
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK3_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_CK3(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_CK3
!#endif
!        use pm_kind, only: CKG => CK3
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK2_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_CK2(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_CK2
!#endif
!        use pm_kind, only: CKG => CK2
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK1_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_CK1(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_CK1
!#endif
!        use pm_kind, only: CKG => CK1
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_RK5(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_RK4(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_RK3(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_RK2(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D0_RK1(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D0_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)                    :: pattern
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_SK5(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_SK4(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_SK3(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_SK2(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_SK1(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_IK5(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_IK5
!#endif
!        use pm_kind, only: IKG => IK5
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_IK4(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_IK4
!#endif
!        use pm_kind, only: IKG => IK4
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_IK3(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_IK3
!#endif
!        use pm_kind, only: IKG => IK3
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_IK2(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_IK2
!#endif
!        use pm_kind, only: IKG => IK2
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_IK1(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_IK1
!#endif
!        use pm_kind, only: IKG => IK1
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if LK5_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_LK5(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_LK5
!#endif
!        use pm_kind, only: LKG => LK5
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK4_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_LK4(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_LK4
!#endif
!        use pm_kind, only: LKG => LK4
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK3_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_LK3(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_LK3
!#endif
!        use pm_kind, only: LKG => LK3
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK2_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_LK2(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_LK2
!#endif
!        use pm_kind, only: LKG => LK2
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK1_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_LK1(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_LK1
!#endif
!        use pm_kind, only: LKG => LK1
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if CK5_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_CK5(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_CK5
!#endif
!        use pm_kind, only: CKG => CK5
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK4_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_CK4(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_CK4
!#endif
!        use pm_kind, only: CKG => CK4
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK3_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_CK3(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_CK3
!#endif
!        use pm_kind, only: CKG => CK3
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK2_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_CK2(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_CK2
!#endif
!        use pm_kind, only: CKG => CK2
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK1_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_CK1(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_CK1
!#endif
!        use pm_kind, only: CKG => CK1
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_RK5(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_RK4(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_RK3(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_RK2(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE module subroutine setRemovedDefComDefIns_D1_D1_RK1(ArrayRemoved, array, pattern)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComDefIns_D1_D1_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
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
!    module subroutine setRemovedCusComDefIns_D1_D1_SK5(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_SK4(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_SK3(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_SK2(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_SK1(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_IK5(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_IK5
!#endif
!        use pm_kind, only: IKG => IK5
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_IK4(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_IK4
!#endif
!        use pm_kind, only: IKG => IK4
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_IK3(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_IK3
!#endif
!        use pm_kind, only: IKG => IK3
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_IK2(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_IK2
!#endif
!        use pm_kind, only: IKG => IK2
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_IK1(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_IK1
!#endif
!        use pm_kind, only: IKG => IK1
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if LK5_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_LK5(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_LK5
!#endif
!        use pm_kind, only: LKG => LK5
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK4_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_LK4(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_LK4
!#endif
!        use pm_kind, only: LKG => LK4
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK3_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_LK3(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_LK3
!#endif
!        use pm_kind, only: LKG => LK3
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK2_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_LK2(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_LK2
!#endif
!        use pm_kind, only: LKG => LK2
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK1_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_LK1(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_LK1
!#endif
!        use pm_kind, only: LKG => LK1
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if CK5_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_CK5(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_CK5
!#endif
!        use pm_kind, only: CKG => CK5
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK4_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_CK4(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_CK4
!#endif
!        use pm_kind, only: CKG => CK4
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK3_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_CK3(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_CK3
!#endif
!        use pm_kind, only: CKG => CK3
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK2_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_CK2(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_CK2
!#endif
!        use pm_kind, only: CKG => CK2
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK1_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_CK1(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_CK1
!#endif
!        use pm_kind, only: CKG => CK1
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_RK5(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_RK4(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_RK3(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_RK2(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    module subroutine setRemovedCusComDefIns_D1_D1_RK1(ArrayRemoved, array, pattern, iseq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComDefIns_D1_D1_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
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
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_SK5(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_SK4(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_SK3(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_SK2(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_SK1(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_IK5(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_IK5
!#endif
!        use pm_kind, only: IKG => IK5
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_IK4(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_IK4
!#endif
!        use pm_kind, only: IKG => IK4
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_IK3(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_IK3
!#endif
!        use pm_kind, only: IKG => IK3
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_IK2(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_IK2
!#endif
!        use pm_kind, only: IKG => IK2
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_IK1(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_IK1
!#endif
!        use pm_kind, only: IKG => IK1
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if LK5_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_LK5(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_LK5
!#endif
!        use pm_kind, only: LKG => LK5
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK4_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_LK4(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_LK4
!#endif
!        use pm_kind, only: LKG => LK4
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK3_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_LK3(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_LK3
!#endif
!        use pm_kind, only: LKG => LK3
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK2_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_LK2(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_LK2
!#endif
!        use pm_kind, only: LKG => LK2
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK1_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_LK1(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_LK1
!#endif
!        use pm_kind, only: LKG => LK1
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if CK5_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_CK5(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_CK5
!#endif
!        use pm_kind, only: CKG => CK5
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK4_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_CK4(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_CK4
!#endif
!        use pm_kind, only: CKG => CK4
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK3_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_CK3(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_CK3
!#endif
!        use pm_kind, only: CKG => CK3
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK2_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_CK2(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_CK2
!#endif
!        use pm_kind, only: CKG => CK2
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK1_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_CK1(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_CK1
!#endif
!        use pm_kind, only: CKG => CK1
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_RK5(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_RK4(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_RK3(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_RK2(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE module subroutine setRemovedDefComCusIns_D1_D1_RK1(ArrayRemoved, array, pattern, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedDefComCusIns_D1_D1_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
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
!    module subroutine setRemovedCusComCusIns_D1_D1_SK5(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_SK4(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_SK3(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_SK2(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_SK1(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        character(:,SKG)        , intent(inout) , allocatable   :: array(:)
!        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        character(*,SKG)                        , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_IK5(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_IK5
!#endif
!        use pm_kind, only: IKG => IK5
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_IK4(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_IK4
!#endif
!        use pm_kind, only: IKG => IK4
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_IK3(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_IK3
!#endif
!        use pm_kind, only: IKG => IK3
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_IK2(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_IK2
!#endif
!        use pm_kind, only: IKG => IK2
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_IK1(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_IK1
!#endif
!        use pm_kind, only: IKG => IK1
!        integer(IKG)            , intent(inout) , allocatable   :: array(:)
!        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        integer(IKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if LK5_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_LK5(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_LK5
!#endif
!        use pm_kind, only: LKG => LK5
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK4_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_LK4(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_LK4
!#endif
!        use pm_kind, only: LKG => LK4
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK3_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_LK3(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_LK3
!#endif
!        use pm_kind, only: LKG => LK3
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK2_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_LK2(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_LK2
!#endif
!        use pm_kind, only: LKG => LK2
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if LK1_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_LK1(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_LK1
!#endif
!        use pm_kind, only: LKG => LK1
!        logical(LKG)            , intent(inout) , allocatable   :: array(:)
!        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        logical(LKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if CK5_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_CK5(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_CK5
!#endif
!        use pm_kind, only: CKG => CK5
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK4_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_CK4(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_CK4
!#endif
!        use pm_kind, only: CKG => CK4
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK3_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_CK3(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_CK3
!#endif
!        use pm_kind, only: CKG => CK3
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK2_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_CK2(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_CK2
!#endif
!        use pm_kind, only: CKG => CK2
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if CK1_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_CK1(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_CK1
!#endif
!        use pm_kind, only: CKG => CK1
!        complex(CKG)            , intent(inout) , allocatable   :: array(:)
!        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        complex(CKG)            , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_RK5(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_RK4(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_RK3(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_RK2(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    module subroutine setRemovedCusComCusIns_D1_D1_RK1(ArrayRemoved, array, pattern, iseq, instance, sorted, unique)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setRemovedCusComCusIns_D1_D1_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)               , intent(inout) , allocatable   :: array(:)
!        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
!        procedure(logical(LK))                                  :: iseq
!        integer(IK)             , intent(in)    , contiguous    :: instance(:)
!        logical(LK)             , intent(in)    , optional      :: sorted
!        logical(LK)             , intent(in)    , optional      :: unique
!        real(RKG)               , intent(out)   , allocatable   :: ArrayRemoved(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayRemove ! LCOV_EXCL_LINE